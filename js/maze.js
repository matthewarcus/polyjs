// Take a polyhedron (should be connected with at most 2 faces meeting
// at an edge) described in OFF format, compute graph topology and use
// thingies algorithm to generate a set of edges defining a spanning
// tree for the dual graph; ie. a maze.

PolyContext.prototype.mazify = function(off,options) {
    // off format is very simple:
    // vertices: a list of points, specified as 3 coordinates.
    // faces: a list of faces, specified as a list of vertices,
    // as indexes into the vertices list.

    if (off.mazified) return off;
    console.log("mazifying")
    var vertices = off.vertices
    var faces = off.faces
    var newfaces = []
    // Walls are pairs of vertices, which separate regions
    var wmap = new Map() // Build up map of wall -> face here
    for (var i = 0; i < faces.length; i++) {
        var face = faces[i]
        var vlist = face.vlist
        if (vlist.length > 2) {
            var faceindex = newfaces.length
            newfaces.push(face)
            for (var j = 0; j < vlist.length; j++) {
                var v1 = vlist[j]
                var v2 = vlist[(j+1)%vlist.length]
                var key;
                // Code walls as integers
                if (v1 < v2) key = (v1 << 15) + v2;
                else key = (v2 << 15) + v1;
                if (!wmap[key]) wmap[key] = [];
                wmap[key].push(faceindex);
            }
        }
    }
    // Build up a list of walls:
    // [ cell1, cell2, vertex1, vertex2 ]
    var walls = []
    for (var k in wmap) {
        var v1 = k >> 15
        var v2 = k & ((1<<15)-1)
        var clist = wmap[k]
        if (clist.length == 1) {
            //  A border edge, so just add to result
            newfaces.push({ vlist: [v1,v2] });
        } else if (clist.length == 2) {
            walls.push([clist[0],clist[1],v1,v2])
        } else {
            console.assert(false,"Wall separates 3 or more regions");
        }
    }

    // Union find
    var regions = []
    for (var i = 0; i < newfaces.length; i++) {
        regions[i] = i
    }

    // Find the region for face n, and update
    // the region array.
    function region(n) {
        if (n == regions[n]) return n;
        var r = region(regions[n])
        regions[n] = r
        return r
    }

    function random(n) {
        var r = Math.random()
        return Math.floor(r * n)
    }

    // Construct maze: a straightforward
    // implementation of thingies algorithm
    for (var i = 0; i < walls.length; i++) {
        var j = i + random(walls.length - i)
        // Select a random wall to consider.
        var w = walls[j]
        walls[j] = walls[i]
        walls[i] = w
        var r0 = region(w[0])
        var r1 = region(w[1])
        if (r0 == r1) {
            // Include this wall in the result, the
            // two sides are already connected.
            newfaces.push({ vlist: [w[2],w[3]] });
        } else {
            // Sides not connected, so connect them
            // and coalesce regions
            regions[r1] = r0;
        }
    }

    off.faces = newfaces
    off.mazified = true
    return off;
}
