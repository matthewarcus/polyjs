"use strict";

// TODO:
// *Snub duals
// *Store separate zrot value
// *Reset zrot (and other values)
// *7 point region drawing
// *Suppress face display
// *Add gradation from centre color mode
// *Sierpinski triangles
// *Reuse geometry when drawing compounds
// *Add check that we have a valid Schwarz triangle
// *Invert though midsphere (minimum of distance to midedge points)
// *Superimpose dual (needs midsphere inversion)
// *Exploding faces
// *Animated inversion
// *Better animated inversion
// *Step to next position
// *Explicitly specify trilinear coords
// *Get rid of global variables
// *Don't draw degenerate sectors of faces
// *Modularize
// *Use proper options object rather than globals
// Reuse point geometry between compounds
// Fix snubcentre for degenerate snub triangles
// Table driven options
// Better separation of generation and display
// Separate logical face structure from display structure
// Hollow faces
// Neo filling
// zrotation should be function of time
// Invertible snub duals 
// Snub stellations
// Dual of stellations
// Stellation of duals
// More generic dual generation
// Stellation parameters
// Make URL for current scene
// More compact URL scheme
// Better colors
// Symmetric colors
// Use homogeneous coordinates
// The Fourth Dimension

// Our state object.
function PolyContext(options) {
    // Options, modifiable during runtime
    this.dostellate = false;
    this.dosnubify = false;
    this.drawtype = 0;
    this.invertinc = -0.01;
    this.ifact = 0;
    this.donormalize = false;
    this.docompound = false;
    this.dorotate = false;
    this.dozrot = false;
    this.dorotonly = false;
    this.colorstyle = 0;
    this.hideface = [false,false,false,false]
    this.zrotation = 0;
    this.tridepth = 0;
    this.regionstyle = 0;
    this.explode = 0;
    this.dohemi = false;
    this.animstep = 10; // Seconds per animation step

    // Take a walk around the Schwarz triangle
    this.tours = [
        [], // Can be set from query params
        [ [1,0,0],[1,1,0],[0,1,0],[0,1,1],[0,0,1],[1,0,1],[1,0,0],
          [1,1,1],[1,1,0],[0,1,0],[1,1,1],[0,1,1],[0,0,1],[1,1,1],[1,0,1] ],
        [ [2,-1,0], [1,0,0], [0,1,0], [-1,2,0],
          [0,2,-1], [0,1,0], [0,0,1], [0,-1,2],
          [-1,0,2], [0,0,1], [1,0,0], [2,0,-1] ],
        [ [2,-1,0],[0,-1,2],[-1,0,2], [-1,2,0],[0,2,-1],[2,0,-1] ]
    ];
    this.tournum = 1;

    this.initrunning = false;
    this.initz = 3;
    this.initsym = [2,3,5];
    this.angles  = [2,3,5];
    this.initt = 0;

    this.processoptions(options);
}

PolyContext.Stopwatch = function(init,step,running) {
    // When running, reported time is difference between Date.now()
    // and stored time. When stopped, reported time is just the
    // stored time.
    this.running = running;
    this.scale = step*1000;  // Return time in seconds
    this.setTime(init || 0); // Different semantics in running and non-running states
}

PolyContext.Stopwatch.prototype.toggle = function() { 
    this.running = !this.running;
    this.time = Date.now()-this.time;
}
PolyContext.Stopwatch.prototype.setTime = function(t) { 
    if (this.running) this.time = Date.now()-t*this.scale;
    else this.time = t*this.scale;
}
PolyContext.Stopwatch.prototype.incTime = function(t) { 
    if (this.running) {
        this.time -= t*this.scale;
    } else {
        // Apply increment and step to next position
        this.time += t*this.scale;
        this.time = this.scale*(Math.floor(this.time/this.scale));
    }
}

PolyContext.Stopwatch.prototype.getTime = function() {
    if (this.running) return (Date.now()-this.time)/this.scale;
    else return this.time/this.scale;
}

PolyContext.prototype.makeinfostring = function() {
    // TBD: put more here
    var s = [];
    s.push("tri = [" + this.tri + "]");
    s.push("bary = [" + this.bary + "]");
    if (this.docompound) {
        s.push("z = " + this.zrotation);
    }
    return s.join("; ");
}

// Given a time t, interpolate between the list of points
PolyContext.prototype.interpolatepoint = function(t,points) {
    t = t%points.length;
    if (t < 0) t += points.length;
    var s = t%1;
    t = Math.floor(t);
    var start = points[t];
    var end = points[(t+1)%points.length];
    var res = [];
    for (var i = 0; i < 3; i++) {
        res[i] = (1-s)*start[i] + s*end[i];
    }
    return res;
}

PolyContext.prototype.makeangles = function(a,b,c,d,e,f) {
    b = b || 1; d = d || 1; f = f || 1;
    a = Number(a); b = Number(b);
    c = Number(c); d = Number(d);
    e = Number(e); f = Number(f);
    // Check b/a + d/c + f/e > 1
    // ie. bce + dea + fac > ace
    if (a < b || c < d || e < f ||
        b*c*e + d*e*a + f*a*c <= a*c*e) {
        alert("Not a Schwarz Triangle");
    }
    return [a/b,c/d,e/f];
}

PolyContext.prototype.processoptions = function(options) {
    // TBD: A more compact parameter format would be good.
    var args = options.split('&');
    var context = this; // Need to see this in the forEach function.
    args.forEach(function (arg) {
        //console.log("Doing parameter '" + arg + "'");
        var matches;
        if (matches = arg.match(/^args=([\d]+)(?:\/([\d]+))?:([\d]+)(?:\/([\d]+))?:([\d]+)(?:\/([\d]+))?$/)) {
            context.angles = context.makeangles(matches[1],matches[2],matches[3],matches[4],matches[5],matches[6]);
        } else if (matches = arg.match(/^sym=([\d]+):([\d]+):([\d]+)$/)) {
            context.initsym = context.makeangles(matches[1],1,matches[2],1,matches[3],1);
        } else if (matches = arg.match(/^tri=([\d.]+):([\d.]+):([\d.]+)$/)) {
            context.tours[0].push([Number(matches[1]),Number(matches[2]),Number(matches[3])]);
            context.tournum = 0;
        } else if (matches = arg.match(/^zrot=([\d]+)$/)) {
            context.zrotation = Math.PI/Number(matches[1]);
        } else if (matches = arg.match(/^zrot=([\d]+.[\d]+)$/)) {
            context.zrotation = Number(matches[1]);
        } else if (matches = arg.match(/^colorstyle=([\d]+)$/)) {
            context.colorstyle = Number(matches[1]);
        } else if (matches = arg.match(/^hide=([\d]+)$/)) {
            context.hideface[Number(matches[1])-1] = true;
        } else if (matches = arg.match(/^t=([\d.]+)$/)) {
            context.initt = Number(matches[1]);
        } else if (matches = arg.match(/^z=([\d.]+)$/)) {
            context.initz = Number(matches[1]);
        } else if (matches = arg.match(/^tridepth=([\d]+)$/)) {
            context.tridepth = Number(matches[1]);
        } else if (matches = arg.match(/^tour=([\d]+)$/)) {
            context.tournum = Number(matches[1]);
        } else if (matches = arg.match(/^compound$/)) {
            context.docompound = true;
        } else if (matches = arg.match(/^snub$/)) {
            context.dosnubify = true;
        } else if (matches = arg.match(/^stellate$/)) {
            context.dostellate = true;
        } else if (matches = arg.match(/^both$/)) {
            context.drawtype = 1;
        } else if (matches = arg.match(/^dual$/)) {
            context.drawtype = 2;
        } else if (matches = arg.match(/^invert$/)) {
            context.invertinc = -context.invertinc;
            context.ifact = 1;
        } else if (matches = arg.match(/^dozrot$/)) {
            context.dozrot = true;
        } else if (matches = arg.match(/^dorotonly$/)) {
            context.dorotonly = true;
        } else if (matches = arg.match(/^dorotate$/)) {
            context.dorotate = true;
        } else if (matches = arg.match(/^hemi$/)) {
            context.dohemi = true;
        } else if (matches = arg.match(/^run$/)) {
            context.initrunning = true;
        } else {
            console.log("Ignoring parameter '" + arg + "'");
        }
    });
    if (this.tours[0].length == 0) {
        this.tours[0].push([1,1,1]);
    }
}

PolyContext.prototype.handleKey = function(key) {
    var handled = true;
    switch(key) {
    case ' ':
        this.stopwatch.toggle();
        break;
    case 'z':
        this.dozrot = !this.dozrot;
        break;
    case 'a':
        this.zrotation = 0;
        break;
    case '[':
        this.stopwatch.incTime(-1);
        break;
    case ']':
        this.stopwatch.incTime(1);
        break;
    case 'p':
        this.stopwatch.setTime(0);
        break;
    case 'r':
        this.dorotate = !this.dorotate;
        break;
    case 'h':
        this.dohemi = !this.dohemi;
        break;
    case 's':
        this.dosnubify = !this.dosnubify;
        break;
    case 'q':
        this.dostellate = !this.dostellate;
        break;
    case 'e':
        this.explode += 1/16;
        break;
    case 'w':
        this.explode -= 1/16;
        break;
    case 'd':
        this.drawtype = (this.drawtype+1)%3
        break;
    case 'n':
        this.donormalize = !this.donormalize;
        break;
    case 'i':
        this.invertinc = -this.invertinc;
        break;
    case 'f':
        this.colorstyle = (this.colorstyle+1) % this.numcolorstyles;
        this.geometry.colorsNeedUpdate = true;
        break;
    case 't':
        this.tournum = (this.tournum+1)%this.tours.length;
        this.stopwatch.setTime(0);
        break;
    case 'u':
        this.rotatecolors();
        this.geometry.colorsNeedUpdate = true;
        break;
    case 'c':
        this.docompound = !this.docompound;
        break;
    case 'x':
        this.dorotonly = !this.dorotonly;
        break;
    case '=':
        this.tridepth = this.tridepth+1;
        break;
    case '-':
        if (this.tridepth > 0) this.tridepth = this.tridepth-1;
        break;
    case 'y':
        this.regionstyle = (this.regionstyle+1)%2;
        break;
    case '?':
        alert(this.makeinfostring());
        break;
    case '1':
        this.hideface[0] = !this.hideface[0];
        break;
    case '2':
        this.hideface[1] = !this.hideface[1];
        break;
    case '3':
        this.hideface[2] = !this.hideface[2];
        break;
    case '4':
        this.hideface[3] = !this.hideface[3];
        break;
    default:
        handled = false;
    }
    return handled;
}

// Finally, draw a point (via. the context.drawpoint0 function).
// This is sort of time-critical so inline a lot of stuff.
// Currently, we call this function several times for each actual point
// (unless we are doing an exploded view).
// so there is some significant optimizations possible here.
// Return the index in the geometries list of points.
PolyContext.prototype.drawpoint = function(p,offset) {
    var Vector = Geometry.Vector;
    var w = p[3] || 1; // Homogeneous coords
    // For fun, do the explode before the inversion
    if (this.donormalize) w = Vector.length(p);
    var ifact = this.ifact;
    if (ifact != 0) {
        //var q = Vector.invert(p,this.midsphere);
        var k = Vector.dot(p,p)/this.midsphere;
        if (ifact == 1) w *= k;
        else w *= (1-ifact)+ifact*k;
    }
    var x = p[0]/w, y = p[1]/w, z = p[2]/w;
    var explode = this.explode;
    if (offset && explode != 0) {
        x += explode*offset[0];
        y += explode*offset[1];
        z += explode*offset[2];
    }
    if (this.transmat) {
        p = [x,y,z];
        x = Vector.dot(this.transmat[0],p);
        y = Vector.dot(this.transmat[1],p);
        z = Vector.dot(this.transmat[2],p);
    }
    return this.drawpoint0(x,y,z);
}

// TBD: encapsulate the various drawing options
// For n > 0, draw a Sierpinski triangle (or some
// variation thereof).
PolyContext.prototype.drawtriangle = function(p,q,r,offset,type,i,n) {
    var Vector = Geometry.Vector;
    if (n == 0) {
        var index0 = this.drawpoint(p,offset);
        var index1 = this.drawpoint(q,offset);
        var index2 = this.drawpoint(r,offset);
        this.drawtriangle0(index0,index1,index2,type,i);
    } else {
        var p1 = Vector.mid(p,q);
        var q1 = Vector.mid(q,r);
        var r1 = Vector.mid(r,p);
        // TBD: Different triangulation styles here
        // k = 0,1,2,3,n%4,random etc.
        //var k = n%4;
        var k = 3;
        if (k != 0) this.drawtriangle(p,p1,r1,offset,type,i,n-1);
        if (k != 1) this.drawtriangle(p1,q,q1,offset,type,i,n-1);
        if (k != 2) this.drawtriangle(r1,q1,r,offset,type,i,n-1);
        if (k != 3) this.drawtriangle(p1,q1,r1,offset,type,i,n-1);
    }
}

PolyContext.prototype.drawface = function(centre,plist,facetype,i,tridepth) {
    var Vector = Geometry.Vector;
    if (this.drawface0) {
        this.drawface0(plist);
    } else {
        if (plist.length == 3) {
            this.drawtriangle(plist[0],plist[1],plist[2],
                              centre,facetype,i,tridepth);
        } else {
            for (var j = 0; j < plist.length; j++) {
                var p1 = plist[j];
                var p2 = plist[(j+1)%plist.length];
                if (Vector.taxi(p1,p2) > 1e-4) {
                    this.drawtriangle(centre,p1,p2,
                                      centre,facetype,i,tridepth);
                }
            }
        }
        if (this.dohemi) {
            for (var j = 0; j < plist.length; j++) {
                var p1 = plist[j];
                var p2 = plist[(j+1)%plist.length];
                if (Vector.taxi(p1,p2) > 1e-4) {
                    this.drawtriangle([0,0,0],p1,p2,
                                      centre,4,i,tridepth);
                }
            }
        }
    }
}

// Draw the Schwarz triangles (fundamental regions) as triangles.
PolyContext.prototype.drawregions = function() {
    var Vector = Geometry.Vector;
    var schwarz = this.schwarz;
    var facedata = this.facedata;
    var points = schwarz.points;
    var regions = schwarz.regions;
    var adjacent = schwarz.adjacent;
    var regionpoints = facedata.regionpoints;
    for (var i = 0; i < regions.length; i++) {
        var region = regions[i];
        var plist = [];
        if (region[3] == 0 || !this.dosnubify) {
            var snubcentre = facedata.snubcentre;
            for (var j = 0; j < 3; j++) {
                plist.push(Vector.mul(points[region[j]], facedata.facedistances[j]));
                if (this.dosnubify && snubcentre) {
                    // This isn't the 'regionpoint' but the 'snubpoint' ie. the
                    // centre of the snub triangle for the face.
                    var adj = adjacent[i][(j+2)%3];
                    plist.push(schwarz.applybary(snubcentre,adj));
                } else if (this.regionstyle == 0) {
                    // Add in edge points for "correct" inversion
                    var edgecentre = schwarz.applybary(facedata.edgecentres[j],i);
                    plist.push(edgecentre);
                }
            }
            for (var j = 0; j < plist.length; j++) {
                plist[j] = Vector.invert(plist[j],this.midsphere);
            }
            var facetype = 4+region[3];
            var centre = Vector.invert(regionpoints[i],this.midsphere);
            this.drawface(centre,plist,facetype,i,this.tridepth)
        }
    }
}

PolyContext.prototype.drawfaces = function() {
    var Vector = Geometry.Vector;
    var schwarz = this.schwarz;
    var facedata = this.facedata;
    var points = schwarz.points;
    var regions = schwarz.regions;
    var regionpoints = facedata.regionpoints;
    // Now go through the facepoint lists for each type of face
    for (var type = 0; type < 3; type++) {
        if (!facedata.faces[type] || this.hideface[type]) continue;
        if (this.dostellate) {
            var facepoints = schwarz.faces[type];
            for (var i = 0; i < facepoints.length; i++) {
                var plist = facepoints[i];
                if (plist) {
                    // Add a centre vertex, this is a corner of a Schwarz triangle.
                    // facedistances[type] tells us how far it is from the origin (the face
                    // centres, ie. the region points are at distance 1.
                    var centre = Vector.mul(points[i],facedata.facedistances[type])                    
                    var facelines = facedata.facelines[type];
                    for (var j = 0; j < facelines.length; j++) {
                        // Triangle from centre to two vertices on plist
                        var p1 = schwarz.applybary(facelines[j][0],plist[0]);
                        var p2 = schwarz.applybary(facelines[j][1],plist[0]);
                        this.drawtriangle(centre,p1,p2,centre,type,i,this.tridepth);
                    }
                }
            }
        } else {
            var facepoints = schwarz.faces[type];
            for (var i = 0; i < facepoints.length; i++) {
                var fpoints = facepoints[i];
                if (fpoints && (!this.dosnubify || fpoints.length > 4)) {
                    // Add a centre vertex, this is a corner of a Schwarz triangle.
                    // facedistances[type] tells us how far it is from the origin (the face
                    // centres, ie. the region points are at distance 1.
                    var centre = Vector.mul(points[i],facedata.facedistances[type]);
                    // For a snub polyhedron, we only include even (+ve) regions.
                    var inc = this.dosnubify ? 2 : 1;
                    var plist = [];
                    for (var j = 0; j < fpoints.length; j+=inc) {
                        plist.push(regionpoints[fpoints[j]]);
                    }
                    this.drawface(centre,plist,type,i,this.tridepth);
                }
            }
        }
    }
    if (this.dosnubify && !this.dostellate && !this.hideface[3] && facedata.snubcentre) {
        // Now draw a snub triangle for -ve regions
        for (var i = 0; i < regions.length; i++) {
            var t = regions[i];
            if (t[3] == 1) {
                var centre = schwarz.applybary(facedata.snubcentre, i);
                var plist = [regionpoints[schwarz.adjacent[i][0]],
                             regionpoints[schwarz.adjacent[i][1]],
                             regionpoints[schwarz.adjacent[i][2]]];
                this.drawface(centre,plist,3,i,this.tridepth);
            }
        }
    }
}

// Set up context for drawing a polyhedron in the geometry object,
// based on the fundamental regions info in schwarz, using
// trilinear coords tri, and with the colorface coloring function.
PolyContext.prototype.setup = function(tri) {
    var Vector = Geometry.Vector;
    var PointSet = Geometry.PointSet;
    this.npoints = 0;
    this.nfaces = 0;
    this.needclone = false;
    this.compound = 0;
    this.pointset = [];
    this.transmat = null;

    var schwarz = this.schwarz;
    var regions = schwarz.regions;
    var points = schwarz.points;
    // Convert trilinear coords to barycentric.
    if (!this.dosnubify) {
        var bary = schwarz.tri2bary(tri);
    } else {
        // For snub polyhedra, we use a bary coord so
        // that 1,1,1 maps through to snub point.
        var bary = [tri[0]*schwarz.snuba[0],
                    tri[1]*schwarz.snuba[1],
                    tri[2]*schwarz.snuba[2]];
        // Scale to unit sphere
        var s = schwarz.applybary(bary,0);
        bary = Vector.div(bary,Vector.length(s));
    }
    var facedata = schwarz.makefacedata(bary);
    if (this.docompound) {
        var regionpoints = facedata.regionpoints;
        // Set up the initial pointset (for compound generation).
        var pointset = [];
        for (var i = 0; i < regionpoints.length; i++) {
            if (!this.dosnubify || regions[i][3] == 0) {
                PointSet.add(regionpoints[i],pointset,1e-6);
            }
        }
        PointSet.sort(pointset);
        this.pointset = pointset;
    }
    this.tri = tri
    this.bary = bary
    this.facedata = facedata;
    if (this.dosnubify) {
        this.midsphere = facedata.snubsphere;
    } else {
        this.midsphere = facedata.midsphere;
    }
    if (this.dostellate) {
        schwarz.stellate(facedata, this.hideface);
    }
}

PolyContext.prototype.drawpolyhedron = function() {
    // Do the actual drawing
    if (this.drawtype <= 1) this.drawfaces();
    if (this.drawtype >= 1) this.drawregions();
}

PolyContext.prototype.drawcompound = function(tripoint) {
    var Vector = Geometry.Vector;
    var PointSet = Geometry.PointSet;
    // We use a group word to describe a particular transformation
    var gens = ["P","Q","R"];
    var rgens = ["PQ","QR","RP"];

    this.setup(tripoint); // Make basic polyhedron
    this.drawpolyhedron();
    if (this.docompound) {
        // Set up for compounding
        var ip = Vector.zrot(this.symmetry[0], this.zrotation);
        var iq = Vector.zrot(this.symmetry[1], this.zrotation);
        var ir = Vector.zrot(this.symmetry[2], this.zrotation);
        // Now we do a transitive closure operation. Keep a list of sets of
        // points, pull sets out and apply all basic operations (from the
        // symmetry group), add back in any sets we haven't previously seen.
        // Initially we just know about the basic pointset just generated.
        var pointsets = [{trans: "", pointset: this.pointset}];
        for (var pindex = 0, count = 0;
             pindex != pointsets.length && count < 200;
             pindex++, count++) {
            var entry = pointsets[pindex];
            var ops = this.dorotonly ? rgens : gens; // Rotation only or all operations?
            for (var i = 0; i < ops.length; i++) {
                var newset = Geometry.mapply(ops[i],entry.pointset,ip,iq,ir);
                PointSet.sort(newset);
                // Have we seen this pointset before?
                var seen = false;
                var parity = (entry.trans.length + ops[i].length)%2;
                for (var j = 0; j < pointsets.length; j++) {
                    // If we care about chirality (ie. are generating a snub figure)
                    // then don't consider reflected pointsets the same).
                    if (this.dosnubify && pointsets[j].trans.length%2 != parity) continue;
                    if (PointSet.equal(newset,pointsets[j].pointset,1e-6)) {
                        seen = true;
                        break;
                    }
                }
                if (!seen) {
                    // A new compound element. Draw it properly and add to queue.
                    // TBD: we could just reuse the geometry with a rotation
                    // and a change of color.
                    var newtrans = entry.trans + ops[i];
                    this.transmat = Geometry.makematrix(newtrans,ip,iq,ir);
                    this.compound++;
                    this.drawpolyhedron();
                    pointsets.push({trans: newtrans, pointset: newset});
                }
            }
        }
    }
}

PolyContext.prototype.runOnCanvas = function(canvas,width,height) {
    var Vector = Geometry.Vector;
    var context = this; // For benefit of subfunctions
    // http://stackoverflow.com/questions/9899807/three-js-detect-webgl-support-and-fallback-to-regular-canvas
    function webglAvailable(canvas) {
        try {
            return Boolean(window.WebGLRenderingContext && 
                           (canvas.getContext("webgl") || 
                            canvas.getContext("experimental-webgl")));
        } catch(e) { 
            return false;
        } 
    }

    // Set up the scene
    var material = new THREE.MeshPhongMaterial({color: 0xffffff}); 
    material.side = THREE.DoubleSide;

    var geometry = new THREE.Geometry();
    geometry.dynamic = true;
    this.geometry = geometry;
    
    var scene = new THREE.Scene(); 
    var mesh = new THREE.Mesh(geometry, material);
    scene.add(mesh);

    var light = new THREE.DirectionalLight(0xffffff, 0.8);
    light.position.set(0,1,1);
    scene.add(light);
    var light = new THREE.DirectionalLight(0xffffff, 0.5);
    light.position.set(0,-1,0);
    scene.add(light);

    var light = new THREE.AmbientLight(0x404040); // soft white light
    scene.add(light);

    var params = { canvas: canvas, antialias: true };
    this.renderer =
        (webglAvailable(canvas)) ? 
        new THREE.WebGLRenderer(params) : 
        new THREE.CanvasRenderer(params);

    this.renderer.setSize(width,height); 

    this.camera = new THREE.PerspectiveCamera(45, width/height, 0.1, 1000); 
    this.camera.position.z = this.initz; 

    var controls = new THREE.OrbitControls( this.camera, canvas );

    // Keep track of time
    this.stopwatch = new PolyContext.Stopwatch(this.initt,
                                               this.animstep,
                                               this.initrunning);

    // Set this if next display loop should recompute geometry
    var needupdate = true;

    var basiccolors = ["red","blue","yellow","lime","yellow","aqua","silver"];

    var compoundColors = [
        new THREE.Color(0xff0000),new THREE.Color(0x00ff00),new THREE.Color(0x0000ff),
        new THREE.Color(0xffff00),new THREE.Color(0x00ffff),new THREE.Color(0xff00ff),
        new THREE.Color(0x808080),new THREE.Color(0xffffff),
    ];
    var faceColors = basiccolors.map(function(c){ return new THREE.Color(c); });
    var white = new THREE.Color(0xffffff);
    var black = new THREE.Color(0x000000);

    this.rotatecolors = function() {
        function rot(a) {
            var t = a[0];
            for (var i = 0; i < a.length-1; i++) {
                a[i] = a[i+1];
            }
            a[a.length-1] = t;
        }
        rot(compoundColors);
        rot(faceColors);
    }

    this.numcolorstyles = 6;
    var vertexcolorstyle = 3;
    this.colorface = function(face,type,index,compound) {
        //face.color = null;
        //face.vertexColors = null;
        var colorstyle = this.colorstyle;
        if (colorstyle == 0) {
            face.color = compoundColors[compound%compoundColors.length];
        } else if (colorstyle == 1) {
            face.color = faceColors[type%faceColors.length];
        } else if (colorstyle == 2) {
            face.color = compoundColors[index%compoundColors.length];
        } else if (colorstyle == 3) {
            face.vertexColors = compoundColors.slice(0,3);
        } else if (colorstyle == 4) {
            var facecolor = faceColors[type%faceColors.length];
            face.vertexColors = [black,facecolor,facecolor];
        } else if (colorstyle == 5) {
            var facecolor = compoundColors[compound%compoundColors.length];
            face.vertexColors = [black,facecolor,facecolor];
        }
    }

    this.drawpoint0 = function (x,y,z) {
        // If we need to add new points to the geometry
        // make sure that we redraw properly.
        if (this.npoints == this.geometry.vertices.length) {
            this.geometry.vertices.push(new THREE.Vector3());
            this.needclone = true;
        }
        var u = this.geometry.vertices[this.npoints];
        u.x = x; u.y = y; u.z = z;
        return this.npoints++;
    }

    this.drawtriangle0 = function(a,b,c,type,index) {
        if (this.nfaces == this.geometry.faces.length) {
            this.geometry.faces.push(new THREE.Face3());
            this.needclone = true;
        }
        var face = this.geometry.faces[this.nfaces];
        face.a = a; face.b = b; face.c = c;
        this.colorface(face,type,index,this.compound);
        this.nfaces++;
    }

    // TBD: Defining a window event handler here seems wrong
    window.addEventListener("keypress",
                            function (event) {
                                var handled = context.handleKey(String.fromCharCode(event.charCode));
                                if (handled) {
                                    needupdate = true;
                                    event.preventDefault();
                                }
                            },
                            false);

    // Set up constant aspects of polyhedra generation
    // schwarz has all the details about the particular set of
    // fundamental regions we are working with.
    var Schwarz = Geometry.Schwarz;
    this.schwarz = new Schwarz(this.angles);
    //this.schwarz.describe(false);

    var sym = Schwarz.makesymmetry(this.initsym);
    this.symmetry = [Vector.normalize(Vector.cross(sym[1],sym[2])),
                     Vector.normalize(Vector.cross(sym[2],sym[0])),
                     Vector.normalize(Vector.cross(sym[0],sym[1]))];

    this.render = function () {
        // Only render at 25fps
        setTimeout(function() {
            requestAnimationFrame( context.render );
        }, 1000 / 25 );
        if (needupdate) {
            // Make sure we have the right color setting in material
            if (context.colorstyle >= vertexcolorstyle) {
                material.vertexColors = THREE.VertexColors;
            } else {
                material.vertexColors = THREE.FaceColors;
            }
            var points = context.tours[context.tournum%context.tours.length];
            var tripoint = context.interpolatepoint(context.stopwatch.getTime(), points);
            context.drawcompound(tripoint);
            // Post process...
            geometry.computeFaceNormals();
            geometry.verticesNeedUpdate = true;
            geometry.elementsNeedUpdate = true;
            geometry.normalsNeedUpdate = true;
            geometry.colorsNeedUpdate = true;
            //geometry.buffersNeedUpdate = true; // If buffer lengths have changed
            if (context.needclone ||
                geometry.vertices.length != context.npoints ||
                geometry.faces.length != context.nfaces) {
                //console.log("Cloning geometry: " + context.npoints + " " + context.nfaces);
                geometry.vertices.length = context.npoints;
                geometry.faces.length = context.nfaces;
                var newgeometry = geometry.clone()
                scene.remove(mesh);
                geometry.dispose();
                geometry = newgeometry;
                context.geometry = geometry;
                // Reusing the mesh object seems OK
                mesh.geometry = geometry; 
                scene.add(mesh);
            }
        }
        needupdate = false;
        // Both z rotation and inversion should be driven
        // from the time rather than however fast we happen
        // to render.
        if (context.dozrot) {
            // Rotate the symmetry frame
            var zrotinc = 0.0025;
            context.zrotation += zrotinc;
            needupdate = true;
        }
        var ifact = context.ifact;
        ifact += context.invertinc;
        if (ifact > 1) ifact = 1;
        if (ifact < 0) ifact = 0;
        if (ifact != context.ifact) {
            context.ifact = ifact;
            needupdate = true;
        }
        // Do we need to update again?
        if (context.stopwatch.running) {
            needupdate = true;
        }
        // This is handled by THREE.js so we don't need to set needupdate.
        if (context.dorotate) {
            mesh.rotation.x += 0.004; 
            mesh.rotation.y += 0.002;
        }
        context.renderer.render(scene, context.camera);
    };
    this.render(0); 
}

// Static function, direct member of PolyContext.
PolyContext.runOnWindow = function() {
    var options = window.location.search;
    if (options.length > 0) {
        // Strip off leading '?'
        options = options.slice(1)
    }

    var canvas = document.createElement("canvas");
    document.body.appendChild(canvas); 
    var width = window.innerWidth;
    var height = window.innerHeight;
    var context = new PolyContext(options);
    window.addEventListener("resize", function() {
        if (context.renderer) {
            var w = window.innerWidth;
            var h = window.innerHeight;
            context.renderer.setSize(w,h);
            context.camera.aspect = w/h;
            context.camera.updateProjectionMatrix();
        }
    });
    context.runOnCanvas(canvas,width,height);
}
