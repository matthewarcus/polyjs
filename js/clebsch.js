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
// that looks nice & I hope it succeeds in that at least.

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

    var Colors = {
        red: new THREE.Color(1,0,0),
        green: new THREE.Color(0,1,0),
        blue: new THREE.Color(0,0,1),
        yellow: new THREE.Color(1,1,0),
        white: new THREE.Color(1,1,1),
        cyan: new THREE.Color(0,1,1),
        magenta: new THREE.Color(1,0,1),
    };

    // Accumulate the surface points here.
    var points = [];
    var colors = [];

    function cube(x) { return x*x*x; }
    function check(x,y,z,w,color) {
        var t1 = cube(x+y+z+w);
        var t2 = cube(x)+cube(y)+cube(z)+cube(w);
        var scale = Math.abs(x) + Math.abs(y) + Math.abs(z) + Math.abs(w);
        var diff = Math.abs(t1-t2)/scale;
        if (diff >= 1e-5) {
            console.log("Equation fails:",x,y,z,w,diff);
        }
        points.push(new THREE.Vector4(x,y,z,w));
        colors.push(color);
    }

    function perm6(x1,x2,x3,x4,color1,color2) {
        color2 = color2 || color1;
        // even permutations -> color1
        // odd permutations -> color2
        check(x1,x2,x3,x4,color1);
        check(x1,x3,x2,x4,color2);
        check(x1,x3,x4,x2,color1);
        check(x3,x1,x2,x4,color1);
        check(x3,x1,x4,x2,color2);
        check(x3,x4,x1,x2,color1);
    }

    var color1 = Colors.red;
    var color2 = Colors.blue;
    var color3 = Colors.green;
    var color4 = Colors.yellow;
    var color5 = Colors.cyan;
    var color6 = Colors.magenta;
    var color7 = Colors.white;
    
    function solve(x,y) {
        // Find solutions for (x,y,z,w) with w = 1
        // and (x+y+z+w)^3 == x^3+y^3+z^3+w^3
        // Solutions for w = 0 form some of the 27 lines
        // which we draw separately.
        var w = 1;
        var A = x+y+w;
        var B = x*x*x + y*y*y + w;
        // if A is very small, z0 and z1 are very large;
        if (Math.abs(A) < 1e-5) return;
        var a = 3*A;
        var b = 3*A*A;
        var c = A*A*A - B;
        var t = quadratic(a,b,c);
        if (t) {
            var z0 = t[0];
            var z1 = t[1];
            // All permutations of x,y,z,w are solutions,
            // But these 12 give good coverage of the surface
            perm6(x,y,z0,w,color1,color2);
            perm6(x,y,z1,w,color3,color4);
        }
    }
    // Random point in (-Infinity,+Infinity) but
    // clustered around origin.
    function randpoint() {
        var x = 1/(2*Math.random()-1);
        if (x < 0) x += 1;
        if (x > 0) x -= 1;
        return x;
    }
    // Main surface.
    function generate1(N) {
        for (var i = 0; i < N*N; i++) {
            var x1 = randpoint();
            var x2 = randpoint();
            solve(x1,x2);
        }
    }
    // 15 lines
    function generate2(N) {
        for (var i = 0; i < N; i++) {
            var x = randpoint();
            var y = randpoint();
            perm6(x,-x,y,0,color7);
            perm6(x,-x,0,y,color7);
            // Need these 3 extra perms for lines 13-15
            check(x,-x,-y,y,color7);
            check(x,-y,-x,y,color7);
            check(-y,x,-x,y,color7);
        }
    }
    // 12 lines
    function generate3(N) {
        for (var i = 0; i < N; i++) {
            var p = (1 + Math.sqrt(5))/2;
            var x1 = randpoint();
            var x2 = randpoint();
            var x3 = -(x1 + p*x2);
            var x4 = -(p*x1 + x2);
            // Color to get a "Double-Six".
            perm6(x1,x2,x3,x4,color5,color6);
            perm6(x1,x2,x4,x3,color6,color5);
        }
    }

    var N = 100;
    var ranges = [];
    generate3(2*N); // 12 lines
    ranges.push(points.length);
    generate2(2*N); // 15 lines
    ranges.push(points.length);
    generate1(N); // Main surface
    ranges.push(points.length);
    var npoints = points.length;
    var nranges = ranges.length;
    var range = nranges-1;
    
    var renderer;
    var camera;
    var geometry;
    Clebsch.runOnWindow = function(canvas) {
        window.addEventListener("resize", function() {
            if (renderer) {
                var w = window.innerWidth;
                var h = window.innerHeight;
                renderer.setSize(w,h);
                camera.aspect = w/h;
                camera.updateProjectionMatrix();
            }
        });
        var running = true;
        window.addEventListener("keypress", function (event) {
            if (!event.ctrlKey) {
                // Ignore event if control key pressed.
                var c = String.fromCharCode(event.charCode);
                switch(c) {
                case ' ':
                    running = !running;
                    event.preventDefault();
                    break;
                case '[':
                    range = (range+nranges-1)%nranges;
                    geometry.setDrawRange(0,ranges[range]);
                    break;
                case ']':
                    range = (range+1)%nranges;
                    geometry.setDrawRange(0,ranges[range]);
                    break;
                    de
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
        camera.position.z = 6; 

        var controls = new THREE.OrbitControls( camera, canvas );

        // A quaternion represents an (isoclinic) rotation in
        // 4-space - effectively changing the 'plane at infinity'.
        // Start with identity,
        var quat = new THREE.Vector4(1,0,0,0);
        quat.normalize();
        // and rotate the quaternion each time around the loop
        // by multiplication with quinc.
        var quinc = new THREE.Vector4(1,0.005,0.01,0.007);
        quinc.normalize();

        // Set up buffers.
        var vertexarray = new Float32Array(npoints*3);
        var colorarray = new Float32Array(npoints*4);
        // Colors are static
        for (var i = 0; i < npoints; i++) {
            var color = colors[i];
            colorarray[4*i+0] = color.r;
            colorarray[4*i+1] = color.g;
            colorarray[4*i+2] = color.b;
            colorarray[4*i+3] = color.a;
        }
        geometry.addAttribute('position', new THREE.BufferAttribute(vertexarray,3));
        geometry.addAttribute('color', new THREE.BufferAttribute(colorarray,4));

        // Positions are dynamic - we could do all this in the shader, but this
        // will do for now.
        var qtmp = new THREE.Vector4;
        function setvertices() {
            for (var i = 0; i < npoints; i++) {
                qmul(quat,points[i],qtmp);
                var x = qtmp.x, y = qtmp.y, z = qtmp.z, w = qtmp.w;
                var eps = 1e-5;
                // w might be 0, but apart from making bounding box
                // calculations difficult, nothing too awful seems to happen.
                // All the same, clamp to a small non-zero value.
                if (0 <= w && w < eps) w = eps;
                if (-eps < w && w < 0) w = -eps;
                x /= w; y /= w; z /= w;
                vertexarray[3*i+0] = x;
                vertexarray[3*i+1] = y;
                vertexarray[3*i+2] = z;
            }
            geometry.attributes.position.needsUpdate = true;
        }
        
        var dorender = function () {
            // Do a full render whether we need to or not.
            setTimeout(function() { requestAnimationFrame(dorender); }, 1000 / 25 );
            setvertices();
            if (running) qmul(quat,quinc,quat);
            renderer.render(scene, camera);
        }
        requestAnimationFrame(dorender);
    }
})();
