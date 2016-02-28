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

// Construct a representation of the Clebsch surface as a THREE.js particle system.
// The surface is stored as a set of homogeneous coordinates, for display we
// multiply by a varying quaternion to perform an isoclinic rotation in
// projective space (before doing the usual projection to 3-space).

// The Clebsch surface has the equation (in projective 4-space):

// x0+x1+x2+x3+x4 = 0, x0^3+x1^3+x2^3+x3^3+x4^3 = 0

// Eliminating x0 and renaming a little, we get:

// (x+y+z+w)^3 = x^3 + y^3 + z^3 + w^3  [***]

// for projective 3-space (ie. with 4-element homogeneous coordinates).

// Since equations are homogeneous, can just consider case of w = 1 and w = 0 (plane
// at infinity). Plane at infinity solutions aren't that interesting. For w = 1, we have:

// (x+y+z+1)^3 = x^3 + y^3 + z^3 + 1

// Given values for x and y, we can solve for z easily - the cubes drop out and we just
// have a quadratic equation that can be solved in the usual way:

// 3Az^2 + 3A^2z + A^3 - B = 0

// where A = x+y+1, B = x^3 + y^3 + 1

// This give a set of homogeneous points [x,y,z,w] satisfying [***] and we can
// cyclically permute the coordinates to get further solutions (any permutation
// is a solution in fact).

// All that remains is to project into 3-space - as usual we divide by the w-coordinate,
// but to get different projections, before doing this we rotate in projective space
// by multiplying by a quaternion & varying the quaternion varies the projection.
// (Quaternion [d,-a,-b,-c] puts plane [a,b,c,d] at infinity - alternatively, rotates
// [a,b,c,d] to [0,0,0,1] - it's well known that quaternions can be used to represent
// rotations in 3-space, but they also work for 4-space (with 3-space as a special
// case) - a 4-space rotation is uniquely represented (up to sign) by x -> pxq where p
// and q are unit quaternions (representing 'isoclinic' or 'Clifford' rotations -
// rotations by the same amount about different planes, giving a twisty-turny effect).

// The Clebsch surface of course is famous for its 27 lines, these are visible
// in our animation, though some may be hiding in the plane at infinity or not have
// enough points displayed to be easily discerned.

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

    // Accumulate the surface points here.
    var points = [];

    function cube(x) { return x*x*x; }
    function check(x,y,z,w) {
        var t1 = cube(x+y+z+w);
        var t2 = cube(x)+cube(y)+cube(z)+cube(w);
        var diff = Math.abs(t1-t2);
        if (diff >= 1e-6) {
            console.log("Equation fails:",x,y,z,w,diff);
        }
        points.push(new THREE.Vector4(x,y,z,w));
    }

    function solve(x,y) {
        // Find solutions for (x,y,z,w) with w = 1
        // and (x+y+z+w)^3 == x^3+y^3+z^3+w^3
        // There are solutions for w = 0 but they
        // aren't very interesting & we ignore them.
        var w = 1;
        var A = x+y+1;
        var B = x*x*x + y*y*y + 1;
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
            // we just use the cyclic shifts. Of course,
            // the solutions overlap, but when the points are
            // colored differently, a nice patchwork arises.
            check(x,y,z0,w);
            check(x,y,z1,w);
            check(w,x,y,z0);
            check(w,x,y,z1);
            check(z0,w,x,y);
            check(z1,w,x,y);
            check(y,z0,w,x);
            check(y,z1,w,x);
        }
    }
    // Construct a square grid of solutions.
    function generate(A,inc) {
        var eps = 1e-4;
        for (var x = -A; x < A+eps; x += inc) {
            for (var y = -A; y < A+eps; y += inc) {
                solve(x,y);
            }
        }
    }
    generate(5,0.1);
    var npoints = points.length;
    console.log("Generated:",npoints);

    var renderer;
    var camera;
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
                switch(String.fromCharCode(event.charCode)) {
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
        var geometry = new THREE.BufferGeometry();
        geometry.dynamic = true;
        scene.add(new THREE.Points(geometry,material));
        
        var params = { canvas: canvas, antialias: true };
        var width = window.innerWidth;
        var height = window.innerHeight;
        renderer = new THREE.WebGLRenderer(params);
        renderer.setPixelRatio(window.devicePixelRatio);
        renderer.setSize(width,height); 

        camera = new THREE.PerspectiveCamera(45, width/height, 0.1, 1000);
        camera.position.z = 6; 

        var controls = new THREE.OrbitControls( camera, canvas );

        // Don't need lights for a Points object...
        
        var stats; // = new Stats();
        if (stats) {
	    stats.domElement.style.position = 'absolute';
	    stats.domElement.style.top = '0px';
	    canvas.parentNode.appendChild(stats.domElement);
        }

        // A quaternion represents an (isoclinic) rotation in
        // 4-space - effectively changing the 'plane at infinity'.
        // Start with identity,
        var quat = new THREE.Vector4(1,0,0,0);
        quat.normalize();
        // and rotate the quaternion each time around the loop
        // by multiplication with quinc.
        var quinc = new THREE.Vector4(1,0,0.01,0);
        quinc.normalize();
        var colors = [ new THREE.Color(1,0,0),
                       new THREE.Color(0,1,0),
                       new THREE.Color(0,0,1),
                       new THREE.Color(1,1,0),
                     ];

        // Set up buffers.
        var vertexarray = new Float32Array(npoints*3);
        var colorarray = new Float32Array(npoints*4);
        // Colors are static
        for (var i = 0; i < npoints; i++) {
            var color = colors[i%colors.length];
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
                // w might be 0, but apart from making bounding box
                // calculations difficult, nothing too awful seems to happen.
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
            if (stats) stats.update();
        }
        requestAnimationFrame(dorender);
    }
})();
