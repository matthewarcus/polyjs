"use strict";

// The MIT License (MIT)

// Copyright (c) 2016 Matthew Arcus

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

// Generate a 3 dimension maze based on a lattice grid.

PolyContext.prototype.maze3 = function (off,options) {
    console.log("3D maze");
    if (off.mazified) return;
    
    var dimx = options.dimx || 5;
    var dimy = options.dimy || 4;
    var dimz = options.dimz || 3;

    var walls = []
    var regions = []
    for (var i = 0; i < dimx*dimy*dimz; i++) regions[i] = i;

    function cell(x,y,z) {
        return x + dimx*(y + dimy*z)
    }

    function coords(n) {
        var x = n%dimx
        var t = Math.floor(n/dimx)
        var y = t%dimy
        var z = Math.floor(t/dimy)
        console.assert(cell(x,y,z) == n)
        return [x,y,z]
    }

    // Union find
    function region(n) {
        if (n == regions[n]) return n;
        var r = region(regions[n])
        regions[n] = r
        return r
    }

    function random(n) {
        return Math.floor(Math.random() * n)
    }

    for (var x = 0; x < dimx; x++) {
        for (var y = 0; y < dimy; y++) {
            for (var z = 0; z < dimz; z++) {
                var n = cell(x,y,z)
                if (x < dimx-1) walls.push([n,cell(x+1,y,z),0])
                if (y < dimy-1) walls.push([n,cell(x,y+1,z),1])
                if (z < dimz-1) walls.push([n,cell(x,y,z+1),2])
            }
        }
    }

    // Construct maze with Kruskal's algorithm
    var maze = []
    for (var i = 0; i < walls.length; i++) {
        var j = i + random(walls.length - i)
        var w = walls[j]
        walls[j] = walls[i]
        walls[i] = w
        var r0 = region(w[0])
        var r1 = region(w[1])
        if (r0 == r1) maze.push(w);
        else regions[r1] = r0;
    }

    var vertices = []
    var faces = []

    var Color = THREE.OFFLoader.Utils.Color;
    var vector = THREE.OFFLoader.Utils.vector
    var color1 = Color.red;
    var color2 = Color.yellow;
    var color3 = Color.green;

    function addvertex(x,y,z) {
        vertices.push(vector(x-dimx/2,y-dimy/2,z-dimz/2));
    }

    function addface(vlist,color) {
        faces.push({ vlist: vlist, color: color })
    }
    for (var i = 0; i < maze.length; i++) {
        var n = maze[i][0]
        var m = maze[i][1]
        var c = coords(maze[i][0])
        var x = c[0]
        var y = c[1]
        var z = c[2]
        var vbase = vertices.length
        var eps = 0.1
        if (maze[i][2] == 0) {
            addvertex(x+1,y+eps,z+eps)
            addvertex(x+1,y+1-eps,z+eps)
            addvertex(x+1,y+1-eps,z+1-eps)
            addvertex(x+1,y+eps,z+1-eps)
            addface([vbase,vbase+1,vbase+2,vbase+3],color1)
        } else if (maze[i][2] == 1){
            addvertex(x+eps,y+1,z+eps)
            addvertex(x+1-eps,y+1,z+eps)
            addvertex(x+1-eps,y+1,z+1-eps)
            addvertex(x+eps,y+1,z+1-eps)
            addface([vbase,vbase+1,vbase+2,vbase+3],color2)
        } else if (maze[i][2] == 2){
            addvertex(x+eps,y+eps,z+1)
            addvertex(x+eps,y+1-eps,z+1)
            addvertex(x+1-eps,y+1-eps,z+1)
            addvertex(x+1-eps,y+eps,z+1)
            addface([vbase,vbase+1,vbase+2,vbase+3],color3)
        } else {
            console.assert(false);
        }
    }
    off.vertices = vertices
    off.faces = faces
    off.mazified = true
    return off
}
