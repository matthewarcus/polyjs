"use strict"

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

// Construct a representation of the Clebsch surface as a Three.js Points object
// The surface is stored as a set of homogeneous coordinates, for display we
// multiply by a varying quaternion to perform an isoclinic rotation in
// projective space before doing the usual projection to Euclidean 3-space.

// Usual disclaimers: this could all be hopelessly incorrect, I'm not an algebraic
// geometrist or even a proper mathematician. Main purpose is to make something
// that looks nice.

// Much use made of:
// http://www.mathcurve.com/surfaces/clebsch/clebsch.shtml
// http://www.mathcurve.com/surfaces/cayley/cayley.shtml

var Clebsch = {};
(function () {
    // Quaternion multiplication (of THREE.Vector4)
    // Multiply p and q as quaternions and put the result in r.
    // r may be p or q (or both).
    function qmul(p,q,r) {
        var p0 = p.x, p1 = p.y, p2 = p.z, p3 = p.w;
        var q0 = q.x, q1 = q.y, q2 = q.z, q3 = q.w;
        var r0 = p0*q0 - p1*q1 - p2*q2 - p3*q3;
        var r1 = p0*q1 + p2*q3 - p3*q2 + p1*q0;
        var r2 = p0*q2 + p3*q1 - p1*q3 + p2*q0;
        var r3 = p0*q3 + p1*q2 - p2*q1 + p3*q0;
        r.x = r0; r.y = r1; r.z = r2; r.w = r3;
    }
    function qconj(p) {
        return new THREE.Vector4(p.x,-p.y,-p.z,-p.w);
    }
    // Solve a quadratic equation az^2 + bz + c = 0
    function quadratic(a,b,c) {
        var disc = b*b - 4*a*c;
        if (disc >= 0) {
            var t = Math.sqrt(disc);
            if (b < 0) t = -t;
            var z1 = (-b - t)/(2*a);
            var z0 = c/(a*z1);
            return [z0,z1];
        }
    }

    function initmatrices() {
        var matrices = []
        matrices.push(null);
        // The "classic" Clebsch view.
        var k = 4 // Vertical elongation
        var l = 3 // Size of passages - smaller is larger
        var m = new THREE.Matrix4;
        var n = new THREE.Matrix4;
        m.set(Math.cos(Math.PI/3),Math.sin(Math.PI/3),-1,0,
              -Math.sin(Math.PI/3),Math.cos(Math.PI/3),0,0,
              k, k, k, 0,
              -l, -l, -l, l);
        n.set(1,0,0,0, 0,0,1,0, 0,1,0,0, 0,0,0,1);
        n.multiply(m);
        matrices.push(n);

        // Vertical rotation
        var m = new THREE.Matrix4;
        var n = new THREE.Matrix4;
        m.set(Math.cos(Math.PI/4),-Math.sin(Math.PI/4), 0, 0,
              Math.sin(Math.PI/4),Math.cos(Math.PI/4), 0, 0,
              0, 0, 1, 0,
              0, 0, 0, 1)
        n.set(1,0,0,0,
              0,Math.cos(Math.PI/6),Math.sin(Math.PI/6), 0,
              0,-Math.sin(Math.PI/6),Math.cos(Math.PI/6), 0,
              0,0,0,1);
        n.multiply(m);
        matrices.push(n);

        var m = new THREE.Matrix4;
        var k = 0.5;
        m.set(1,0,0,0, 0,1,0,0, 0,0,1,0, -k,-k,-k,k);
        m.getInverse(m);
        matrices.push(m);

        var m = new THREE.Matrix4;
        m.set(-1,1,1,1, 1,-1,1,1, 1,1,-1,1, -1,-1,-1, 1);
        m.getInverse(m);
        matrices.push(m);

        var m = new THREE.Matrix4;
        m.set(-1,1,1,1, 1,-1,1,1, 1,1,-1,1, 1,1,1,-1);
        m.getInverse(m);
        matrices.push(m);
        return matrices;
    };

    var Colors = {
        red: new THREE.Color(1,0,0),
        green: new THREE.Color(0,1,0),
        blue: new THREE.Color(0,0,1),
        yellow: new THREE.Color(1,1,0),
        white: new THREE.Color(1,1,1),
        cyan: new THREE.Color(0,1,1),
        magenta: new THREE.Color(1,0,1),
        celadon: new THREE.Color(0xace1af),
    };
    var colorschemes = [
        [Colors.red,Colors.blue,Colors.green,Colors.yellow,Colors.white,Colors.cyan,Colors.magenta],
        [Colors.celadon,Colors.celadon,Colors.celadon,Colors.celadon,Colors.red,Colors.yellow,Colors.blue],
        [Colors.red,Colors.white,Colors.blue,Colors.white,Colors.yellow,Colors.yellow,Colors.yellow]
    ];

    var quinc = [
        new THREE.Vector4(1,0.005,0.01,0.007),
        new THREE.Vector4(1,0.0,0.02,0.0),
    ];
    for (var i = 0; i < quinc.length; i++) {
        quinc[i].normalize();
    }

    // Miscellaneous (potentially) variable parameters
    var radius2;
    //var radius2 = 4;

    var quincindex = 0;
    var colorscheme = 0;
    var matrixindex = 0;
    
    // Accumulate the surface points here.
    var points = [];
    var colors = [];
    var npoints = 0;

    function addpoint(x,y,z,w,color) {
        console.assert(color);
        if (points.length <= npoints) {
            points.push(new THREE.Vector4);
            colors.push(null);
        }
        var point = points[npoints];
        point.x = x;
        point.y = y;
        point.z = z;
        point.w = w;
        // Assume we don't need to copy the color
        colors[npoints] = color;
        npoints++;
    }

    function perm6(x1,x2,x3,x4,color1,color2) {
        color2 = color2 || color1;
        // even permutations -> color1
        // odd permutations -> color2
        addpoint(x1,x2,x3,x4,color1);
        addpoint(x1,x3,x2,x4,color2);
        addpoint(x1,x3,x4,x2,color1);
        addpoint(x3,x1,x2,x4,color1);
        addpoint(x3,x1,x4,x2,color2);
        addpoint(x3,x4,x1,x2,color1);
    }

    function perm3(x1,x2,x3,x4,color1,color2) {
        color2 = color2 || color1;
        // Leave x4 = w fixed
        // even permutations -> color1
        // odd permutations -> color2
        addpoint(x1,x2,x3,x4,color1);
        addpoint(x1,x3,x2,x4,color2);
        addpoint(x3,x1,x2,x4,color1);
    }

    var randomvalues = [];
    var randomindex = 0;

    function getrandom() {
        if (randomindex == randomvalues.length) {
            randomvalues.push(Math.random());
        }
        return randomvalues[randomindex++];
    }

    // Random point in (-Infinity,+Infinity) but
    // clustered around origin.
    function randpoint() {
        var x = 1/(2*getrandom()-1);
        if (x < 0) x += 1;
        if (x > 0) x -= 1;
        return x;
    }

    var time = 0;
    var morphing = true;
    
    // Surface generating functions
    function morph() {
        var a = Math.sin(time);
        var K = 0.25 + 4*a*Math.abs(a);
        time += 0.01;
        
        var N = 100;
        var colors = colorschemes[colorscheme];
        // Almost the same as Clebsch, but with w parameter = K,
        // restricted permutations (since w is now special) and,
        // alas, no lines, as I haven't worked out how to do that
        // yet - Coble has a method based on the hexahedral form
        // that might be usable, otherwise, need to code up more
        // group theory & algebraic geometry than I would like right
        // now, fascinating though that would be.
        function solve(x,y,w) {
            // Find solutions for (x,y,z,w) with w = 1
            // and (x+y+z+w)^3 == x^3+y^3+z^3+K*w^3
            var A = x+y+w;
            var B = x*x*x + y*y*y + K*w;
            var a = 3*A;
            var b = 3*A*A;
            var c = A*A*A - B;
            var t = quadratic(a,b,c);
            if (t) {
                var z0 = t[0];
                var z1 = t[1];
                // All permutations of x,y,z are solutions,
                perm3(x,y,z0,w,colors[0],colors[1]);
                perm3(x,y,z1,w,colors[2],colors[3]);
            }
        }
        // Well, we can have these 3 lines at least.
        for (var i = 0; i < N; i++) {
            var x = randpoint();
            var y = randpoint();
            perm3(x,-x,y,0,colors[4]);
        }
        for (var i = 0; i < N*N; i++) {
            var x1 = randpoint();
            var x2 = randpoint();
            var w = 1;
            solve(x1,x2,w);
        }
        return morphing; // Carry on redrawing
    }

    function clebsch() {
        var colors = colorschemes[colorscheme];
        function check(x,y,z,w) {
            function cube(x) { return x*x*x; }
            var t1 = cube(x+y+z+w);
            var t2 = cube(x)+cube(y)+cube(z)+cube(w);
            var scale = Math.abs(x) + Math.abs(y) + Math.abs(z) + Math.abs(w);
            var diff = Math.abs(t1-t2)/scale;
            if (diff >= 1e-4) {
                console.log("Equation fails:",x,y,z,w,diff);
            }
        }

        function solve(x,y,w) {
            // Find solutions for (x,y,z,w) with w = 1
            // and (x+y+z+w)^3 == x^3+y^3+z^3+w^3
            var A = x+y+w;
            var B = x*x*x + y*y*y + w;
            var a = 3*A;
            var b = 3*A*A;
            var c = A*A*A - B;
            var t = quadratic(a,b,c);
            if (t) {
                var z0 = t[0];
                var z1 = t[1];
                check(x,y,z0,w);
                check(x,y,z1,w);
                // All permutations of x,y,z,w are solutions,
                perm6(x,y,z0,w,colors[0],colors[1]);
                perm6(x,y,z1,w,colors[2],colors[3]);
            }
        }
        function clebsch1(N) {
            for (var i = 0; i < N*N; i++) {
                var x1 = randpoint();
                var x2 = randpoint();
                var w = 1;
                // Solutions for w = 0 form some of the 27 lines
                // which we draw separately.
                solve(x1,x2,w);
            }
        }
        // 15 lines
        function clebsch2(N) {
            for (var i = 0; i < N; i++) {
                var x = randpoint();
                var y = randpoint();
                var color = colors[4];
                perm6(x,-x,y,0,color);
                perm6(x,-x,0,y,color);
                // Need these 3 extra perms for lines 13-15
                addpoint(x,-x,-y,y,color);
                addpoint(x,-y,-x,y,color);
                addpoint(-y,x,-x,y,color);
            }
        }
        // 12 lines
        function clebsch3(N) {
            for (var i = 0; i < N; i++) {
                var p = (1 + Math.sqrt(5))/2;
                var x1 = randpoint();
                var x2 = randpoint();
                var x3 = -(x1 + p*x2);
                var x4 = -(p*x1 + x2);
                // Color to get a "Double Six".
                perm6(x1,x2,x3,x4,colors[5],colors[6]);
                perm6(x1,x2,x4,x3,colors[6],colors[5]);
            }
        }
        var N = 100;
        clebsch3(2*N); // 12 lines
        clebsch2(2*N); // 15 lines
        clebsch1(N); // Main surface
    }
    function cayley() {
        var colors = colorschemes[colorscheme];
        function cayley1(N) {
            for (var i = 0; i < N*N; i++) {
                var x = randpoint();
                var y = randpoint();
                var z = 1;
                var w = -1/(1/x + 1/y + 1/z);
                perm6(x,y,z,w,colors[0],colors[1]);
                perm6(x,y,w,z,colors[2],colors[3]);
            }
        }
        function cayley2(N) {
            for (var i = 0; i < N; i++) {
                var color = colors[4];
                var x1 = randpoint();
                var x2 = randpoint();
                perm6(x1,x2,0,0,color,color);
                addpoint(x1,-x1,x2,-x2,color);
                addpoint(x1,x2,-x1,-x2,color);
                addpoint(x2,x1,-x1,-x2,color);
            }
        }
        var N = 100;
        cayley2(4*N); // 9 lines
        cayley1(N); // Main surface
    }

    var surfaces = [clebsch, cayley, morph];
    var surface = 0;

    var vertexarray;
    var colorarray;
    var geometry;

    function setgeometry() {
        npoints = 0;
        var redraw = surfaces[surface](); // Compute points

        // Reset random so we get the same collection of points
        randomindex = 0;


        // Set up buffers.
        if (!vertexarray || vertexarray.length < npoints*3) {
            geometry.dispose();
            vertexarray = new Float32Array(npoints*3);
            colorarray = new Float32Array(npoints*4);
            geometry.addAttribute('position', new THREE.BufferAttribute(vertexarray,3));
            geometry.addAttribute('color', new THREE.BufferAttribute(colorarray,4));
        }
        geometry.attributes.position.needsUpdate = true;
        geometry.attributes.color.needsUpdate = true;
        geometry.setDrawRange(0,npoints);

        for (var i = 0; i < npoints; i++) {
            var color = colors[i];
            colorarray[4*i+0] = color.r;
            colorarray[4*i+1] = color.g;
            colorarray[4*i+2] = color.b;
            colorarray[4*i+3] = color.a;
        }
        return redraw;
    }

    Clebsch.runOnWindow = function(canvas) {
        var renderer;
        var camera;
        var quat;
        var matrices = initmatrices();

        var running = true;

        function processoptions(options) {
            var args = options.split('&');
            args.forEach(function (arg) {
                if (arg == '') return;
                //console.log("Doing parameter '" + arg + "'");
                var matches;
                if (matches = arg.match(/^morphing$/)) {
                    morphing = true;
                } else if (matches = arg.match(/^run$/)) {
                    running = true;
                } else if (matches = arg.match(/^stop$/)) {
                    running = false;
                } else if (matches = arg.match(/^clebsch$/)) {
                    surface = 0;
                } else if (matches = arg.match(/^cayley$/)) {
                    surface = 1;
                } else if (matches = arg.match(/^morph$/)) {
                    surface = 2;
                } else if (matches = arg.match(/^surface=([\d]+)$/)) {
                    surface = Number(matches[1])%surfaces.length;
                } else if (matches = arg.match(/^projection=([\d]+)$/)) {
                    matrixindex = Number(matches[1])%matrices.length;
                } else if (matches = arg.match(/^colors=([\d]+)$/)) {
                    colorscheme = Number(matches[1])%colorschemes.length;
                } else if (matches = arg.match(/^q=([\d]+)$/)) {
                    quincindex = Number(matches[1])%quinc.length;
                } else {
                    console.log("Unknown option " + arg);
                }
            })
        }

        var options = window.location.search;
        if (options.length > 0) {
            // Strip off leading '?'
            options = options.slice(1)
        }
        processoptions(options);
        window.addEventListener("resize", function() {
            if (renderer) {
                var w = window.innerWidth;
                var h = window.innerHeight;
                renderer.setSize(w,h);
                camera.aspect = w/h;
                camera.updateProjectionMatrix();
            }
        });
        window.addEventListener("keypress", function (event) {
            if (!event.ctrlKey) {
                // Ignore event if control key pressed.
                var c = String.fromCharCode(event.charCode);
                switch(c) {
                case 'a':
                    // Decrease object scale
                    scale /= 1.1;
                    break;
                case 'c':
                    // Change color scheme
                    colorscheme = (colorscheme+1)%colorschemes.length;
                    needupdate = true;
                    break;
                case 'm':
                    // Toggle morphing
                    morphing = !morphing;
                    if (morphing) needupdate = true;
                    break;
                case 'p':
                    // Change projection
                    matrixindex = (matrixindex+1)%matrices.length;
                    break;
                case 'q':
                    // Change quaternion
                    quincindex = (quincindex+1)%quinc.length;
                    quat = new THREE.Vector4(1,0,0,0);
                    break;
                case 'r':
                    // Reset
                    quat = new THREE.Vector4(1,0,0,0);
                    break;
                case 's':
                    // Increase object scale
                    scale *= 1.1;
                    break;
                case 't':
                    // Change surface
                    surface = (surface+1)%surfaces.length;
                    needupdate = true;
                    break;
                case ' ':
                    running = !running;
                    event.preventDefault();
                    break;
                }
            }
        }, false);
        
        var scene = new THREE.Scene();
        var material = new THREE.PointsMaterial({
            size: 0.1,
            vertexColors: true,
            map: THREE.ImageUtils.loadTexture("images/ball.png"),
            alphaTest: 0.5, // Clip to texture
        });
        geometry = new THREE.BufferGeometry();
        geometry.dynamic = true;
        scene.add(new THREE.Points(geometry,material));
        
        var params = { canvas: canvas, antialias: true };
        var width = window.innerWidth;
        var height = window.innerHeight;
        renderer = new THREE.WebGLRenderer(params);
        renderer.setPixelRatio(window.devicePixelRatio);
        renderer.setSize(width,height); 
        renderer.setClearColor(new THREE.Color(0,0,0.1));
        camera = new THREE.PerspectiveCamera(45, width/height, 0.1, 1000);
        camera.position.x = 0; 
        camera.position.y = 0; 
        camera.position.z = 6; 

        var controls = new THREE.OrbitControls( camera, canvas );

        // A quaternion represents an (isoclinic) rotation in
        // 4-space - effectively changing the 'plane at infinity'.
        // Start with identity,
        quat = new THREE.Vector4(1,0,0,0);
        quat.normalize();
        // and rotate the quaternion each time around the loop
        // by multiplication with quinc.

        var needupdate = true;
        var scale = 1.0;
        
        // Positions are dynamic - we could do all this in the shader, but this
        // will do for now.
        var qtmp = new THREE.Vector4;
        function setvertices() {
            for (var i = 0; i < npoints; i++) {
                qtmp.copy(points[i]);
                if (matrices[matrixindex]) {
                    qtmp.applyMatrix4(matrices[matrixindex]);
                }
                qmul(quat,qtmp,qtmp);
                //qmul(qtmp,quat,qtmp);
                var x = qtmp.x, y = qtmp.y, z = qtmp.z, w = qtmp.w;
                var eps = 1e-5;
                // w might be 0, but apart from making bounding box
                // calculations difficult, nothing too awful seems to happen.
                // All the same, clamp to a small non-zero value.
                if (0 <= w && w < eps) w = eps;
                if (-eps < w && w < 0) w = -eps;
                w /= scale;
                x /= w; y /= w; z /= w;
                if (radius2 && x*x+z*z > radius2) x = 10000;
                vertexarray[3*i+0] = x;
                vertexarray[3*i+1] = y;
                vertexarray[3*i+2] = z;
            }
            geometry.attributes.position.needsUpdate = true;
        }
        
        var dorender = function () {
            //console.log(npoints, points.length);
            // Do a full render whether we need to or not.
            setTimeout(function() { requestAnimationFrame(dorender); }, 1000 / 25 );
            if (needupdate) {
                needupdate = setgeometry();
            }
            setvertices();
            if (running) {
                qmul(quat,quinc[quincindex],quat);
            }
            renderer.render(scene, camera);
        }
        requestAnimationFrame(dorender);
    }
})();
