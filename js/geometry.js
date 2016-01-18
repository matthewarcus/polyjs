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

// Miscellaneous functions used by the main polyhedra generation code

// The usual namespace object
var Geometry = {};

(function() {
    // Vector utility functions. We could use THREE.js vectors here.
    var Vector = {
        dot: function(u,v) {
            if (!v) v = u;
            var w = 0;
            for (var i = 0; i < u.length; i++) {
                w += u[i]*v[i];
            }
            return w;
        },
        cross: function(u,v) {
            var w = [];
            w[0] = u[1]*v[2]-u[2]*v[1];
            w[1] = u[2]*v[0]-u[0]*v[2];
            w[2] = u[0]*v[1]-u[1]*v[0];
            return w;
        },
        triple: function(u,v,w) {
            return this.dot(u,this.cross(v,w));
        },
        length: function(u) {
            return Math.sqrt(Vector.dot(u,u));
        },
        dist: function(u,v) {
            var s = 0;
            for (var i = 0; i < u.length; i++) {
                var t = u[i]-v[i];
                s += t*t;
            }
            return Math.sqrt(s);
        },
        taxi: function(u,v) {
            var s = 0;
            for (var i = 0; i < u.length; i++) {
                s += Math.abs(u[i]-v[i]);
            }
            return s;
        },
        add: function(u,v) {
            var w = [];
            for (var i = 0; i < u.length; i++) {
                w[i] = u[i]+v[i];
            }
            return w;
        },
        add3: function(u,v,w) { 
            return Vector.add(u,Vector.add(v,w));
        },
        sub: function(u,v) {
            var w = [];
            for (var i = 0; i < u.length; i++) {
                w[i] = u[i]-v[i];
            }
            return w;
        },
        mul: function(u,x) {
            var w = [];
            for (var i = 0; i < u.length; i++) {
                w[i] = u[i]*x;
            }
            return w;
        },
        div: function(u,x) {
            return Vector.mul(u,1/x);
        },
        negate: function(u) {
            return Vector.mul(u,-1);
        },
        normalize: function(u) {
            return Vector.div(u,Vector.length(u));
        },
        mid: function(u,v) {
            var w = new Array(u.length);
            for (var i = 0; i < u.length; i++) {
                w[i] = (u[i]+v[i])/2;
            }
            return w;
        },
        interp: function(u,v,k) {
            var w = new Array(u.length);
            for (var i = 0; i < u.length; i++) {
                w[i] = k*u[i]+(1-k)*v[i];
            }
            return w;
        },
        copy: function(u) {
            var w = [];
            for (var i = 0; i < u.length; i++) {
                w[i] = u[i];
            }
            return w;
        },
        set: function(u,v) {
            console.assert(u.length == 0 || u.length == v.length,
                           "Incompatible array lengths");
            for (var i = 0; i < v.length; i++) {
                u[i] = v[i];
            }
        },
        approxeq: function(u,v,eps) {
            //console.assert(eps);
            //console.assert(u.length == v.length, "Incompatible array lengths");
            for (var i = 0; i < u.length; i++) {
                if (Math.abs(u[i] - v[i]) > eps) return false;
            }
            return true;
        },
        eq: function(u,v) {
            console.assert(u.length == v.length,
                           "Incompatible array lengths");
            for (var i = 0; i < u.length; i++) {
                if (u[i] != v[i]) return false;
            }
            return true;
        },
        // Rotate about z-axis (for compounds)
        zrot: function(p,theta) {
            return [Math.cos(theta)*p[0] + Math.sin(theta)*p[1],
                    -Math.sin(theta)*p[0] + Math.cos(theta)*p[1],
                    p[2]];
        },
        // Rotate about y-axis (for compounds)
        yrot: function(p,theta) {
            return [Math.cos(theta)*p[0] + Math.sin(theta)*p[2],
                    p[1],
                    -Math.sin(theta)*p[0] + Math.cos(theta)*p[2]]
        },
        intersection: function(p,q,r) {
            var P = Vector.cross(q,r);
            var Q = Vector.cross(r,p);
            var R = Vector.cross(p,q);
            var k = Vector.dot(p,P);
            if (Math.abs(k) < 1e-4) return null;
            var A = Vector.dot(p,p);
            var B = Vector.dot(q,q);
            var C = Vector.dot(r,r);
            // Maybe use a 4-vector here?
            var s = Vector.div(Vector.applybary([A,B,C],P,Q,R),k);
            return s;
        },
        reflect: function(p,q,r) {
            // Reflect r in the plane of p and q (and the origin)
            var n = Vector.normalize(Vector.cross(p,q));
            return Vector.sub(r, Vector.mul(n,2*Vector.dot(n,r)));
        },
        planereflect: function(P,r) {
            // Reflect r in the plane (through origin) with normal P
            return Vector.sub(r, Vector.mul(P,2*Vector.dot(P,r)));
        },
        reflectionmatrix: function(p,q) {
            var v0 = reflect(p,q,[1,0,0]);
            var v1 = reflect(p,q,[0,1,0]);
            var v2 = reflect(p,q,[0,0,1]);
            return [v0,v1,v2];
        },
        // "barycentric" here just means interpreting
        // coordinates [a,b,c] as the weighted sum of aA+bB+cC
        // where A,B,C are the vertices of some triangle.
        // "trilinear" means defining a point by its distance
        // from 3 different planes, or equivalently, by the dot product
        // with 3 different vectors.
        // if p,q,r is a (spherical) triangle, then P = q x r, Q = r x p, R = p x q
        // are normals to its sides
        tri2bary: function(tri,p,q,r) {
            var P = Vector.cross(q,r);
            var Q = Vector.cross(r,p);
            var R = Vector.cross(p,q);
            var bary = [Vector.length(P)*tri[0],
                        Vector.length(Q)*tri[1],
                        Vector.length(R)*tri[2]];
            var s = Vector.applybary(bary,p,q,r);
            return Vector.div(bary,Vector.length(s));
        },
        getbary: function(p,q,r,s) {
            var P = Vector.cross(q,r);
            var Q = Vector.cross(r,p);
            var R = Vector.cross(p,q);
            var k = Vector.dot(p,P); // Triple product (p,q,r)
            var A = Vector.dot(s,P);
            var B = Vector.dot(s,Q);
            var C = Vector.dot(s,R);
            return [A/k,B/k,C/k];
        },
        applybary: function(bary,p,q,r) {
            return [bary[0]*p[0]+bary[1]*q[0]+bary[2]*r[0],
                    bary[0]*p[1]+bary[1]*q[1]+bary[2]*r[1],
                    bary[0]*p[2]+bary[1]*q[2]+bary[2]*r[2]];
        },
        // Invert relative to an origin-centred sphere.
        invert: function(p,r2) {
            // return kp where kpp = r2
            var k = Vector.dot(p,p)/r2;
            //return [p[0],p[1],p[2],k];
            return Vector.div(p,k);
        }
}

    // Manage sets of points, ordered on first coordinate with
    // reasonably efficient equality (up to epsilon).
    var PointSet = {};

    // Simple set equality, mainly to check correctness of the more
    // efficient setequal1.
    // There are various pathological cases since approxeq isn't an
    // equivalence relation, but I don't think they arise much in
    // our application.
    function setequal0(s,t,eps) {
        console.assert(s.length == t.length, "Incompatible set lengths");
        for (var i = 0; i < s.length; i++) {
            var found = false;
            for (var j = 0; j < t.length; j++) {
                if (Vector.approxeq(s[i],t[j],eps)) {
                    found = true;
                    break;
                }
            }
            if (!found) return false;
        }
        return true;
    }

    function setequal1(s,t,eps) {
        for (var i = 0, j = 0; i < s.length; i++) {
            while(j < s.length && t[j][0] + eps < s[i][0]) j++;
            if (j == t.length) return false;
            for (var k = j; ; k++) {
                if (k == t.length) return false;
                if (s[i][0] + eps < t[k][0]) return false;
                if (Vector.approxeq(s[i],t[k],eps)) break;
            }
        }
        return true;
    }

    PointSet.equal = function(s,t,eps) {
        var res = setequal1(s,t,eps);
        // if (res !== setequal0(s,t,eps)) {
        //     console.log(JSON.stringify(s));
        //     console.log(JSON.stringify(t));
        //     console.assert(false);
        // }
        return res;
    }

    // Add an element to a set if it's not already there.
    // Don't assume set is ordered at this point.
    PointSet.add = function(a,s,eps) {
        for (var i = 0; i < s.length; i++) {
            if (Vector.approxeq(a,s[i],eps)) return false;
        }
        s.push(a);
        return true;
    }

    PointSet.sort = function(s) {
        s.sort(function(u,v) { return u[0] - v[0]; });
    }

    // Now for some serious 3D vector geometry...

    function makesymmetry(angles) {
        // Spherical triangle cosine rule
        function cosine (A,B,C)
        {
            var acos = Math.acos, cos = Math.cos, sin = Math.sin;
            var a = acos((cos(A) + cos(B)*cos(C))/(sin(B)*sin(C)));
            var b = acos((cos(B) + cos(C)*cos(A))/(sin(C)*sin(A)));
            var c = acos((cos(C) + cos(A)*cos(B))/(sin(A)*sin(B)));
            if (isNaN(a) || isNaN(b) || isNaN(c)) return null;
            else return [a,b,c];
        }

        // Set p,q,r to the vertices of a spherical triangle with
        // sides subtending angles a,b,c.
        function triangle1(a,b,c)
        {
            var p = [0,0,1];
            var q = [0, Math.sin(c), Math.cos(c)];
            var rz = Math.cos(b);
            var ry = (Math.cos(a)-q[2]*rz)/q[1];
            var t = 1-rz*rz-ry*ry;
            if (t < 0) t = 0; // Clamp to avoid rounding problems
            var r = [Math.sqrt(t),ry,rz];
            // This used to reverse q and r to get a +ve OpenGL triangle to start
            // with, but that mucks up the tri coordinates, and anyway the complicated
            // stuff needs drawing both sides.
            return [p,q,r];
        }

        var a = angles[0], b = angles[1], c = angles[2];
        var A = Math.PI/a, B = Math.PI/b, C = Math.PI/c;
        var t = cosine(A,B,C);
        if (!t) return null;
        return triangle1(t[0],t[1],t[2]);
    }

    // Find bary coords of point whose 3 reflections form an equilateral triangle.
    // Fairly standard application of 2-dimensional Newton-Raphson.
    function getsnub(p,q,r) {
        // Solve f(a,b,c) = g(a,b,c) = h(a,b,c)
        // Here f,g,h are distances to 3 sides of triangle. a,b,c are bary coords
        // In fact, we can set a+b+c = 1, so only 2 variables really.
        // Have a vector quantity: [f-g,h-g], which we want to set to [0,0].
        // f(x+dx) = f(x) + F(dx)
        // ie. f(x) + F(dx) = 0 => dx = -inv(F)(f(x))
        function jacobian(f,a,b,eps) {
            // f(a+eps) = f(a-eps) + 2*eps*f'(a) => f'(a) =  (f(a+eps)-f(a-eps))/(2*eps)
            var s0 = f(a+eps,b);
            var s1 = f(a-eps,b);
            var s2 = f(a,b+eps);
            var s3 = f(a,b-eps);
            // df[0]/da df[0]/db
            // df[1]/da df[1]/db
            return [ (s0[0]-s1[0])/(2*eps), (s2[0]-s3[0])/(2*eps),
                     (s0[1]-s1[1])/(2*eps), (s2[1]-s3[1])/(2*eps) ];
        }
        function inv(m) {
            var a = m[0], b = m[1], c = m[2], d = m[3];
            var det = a*d - b*c;
            return [d/det,-b/det,-c/det,a/det];
        }
        function mapply(m,s) {
            var a = m[0], b = m[1], c = m[2], d = m[3];
            var x = s[0], y = s[1];
            return [a*x+b*y,c*x+d*y];
        }
        function refine(f,s) {
            var a = s[0], b = s[1];
            // 0 = f(a+dx) = f(a)+A(dx)
            // f(a) = -A(dx)
            // dx = -inv(A)(f(a))
            var m = inv(jacobian(f,a,b,1e-6));
            var s = f(a,b);
            var dx = mapply(m,s);
            return [a-dx[0],b-dx[1]];
        }
        function f(a,b) {
            var s = Vector.applybary([a,b,1-a-b],p,q,r);
            var p0 = Vector.reflect(q,r,s);
            var q0 = Vector.reflect(r,p,s);
            var r0 = Vector.reflect(p,q,s);
            var d0 = Vector.dist(p0,q0);
            var d1 = Vector.dist(q0,r0);
            var d2 = Vector.dist(r0,p0);
            return [d1-d0,d2-d1];
        }
        // Middle of the triangle
        var s = [1/3,1/3];
        // 6 iterations is more than enough
        for (var i = 0; i < 6; i++) {
            s = refine(f,s);
            //console.log(" -- " + s);
        }
        var t = f.apply(null,s);
        var check = Math.abs(t[0]) + Math.abs(t[1]);
        // Did we find a solution?
        if (check > 1e-6) console.log("Snubification failure: " + t);
        return [s[0],s[1],1-s[0]-s[1]];
    }

    function Schwarz(angles) {
        var t = makesymmetry(angles);
        if (!t) return null;
        // p,q,r are the points of the initialSchwarz triangle
        var p = t[0], q = t[1], r = t[2];
        // P,Q,R are the barycentric coordinates of reflections of p,q,r.
        var P = Vector.getbary(p,q,r,Vector.reflect(q,r,p));
        var Q = Vector.getbary(p,q,r,Vector.reflect(r,p,q));
        var R = Vector.getbary(p,q,r,Vector.reflect(p,q,r));
        this.p = p; this.q = q; this.r = r; 
        this.P = P; this.Q = Q; this.R = R; 

        var points = [p,q,r];
        var regions = [[0,1,2,0]];
        var faces = [[],[],[]];
        var adjacent = [];
        function getpoint(p) {
            var eps = 1e-6;
            for (var i = 0; i < points.length; i++) {
                if (Vector.approxeq(points[i],p,eps)) {
                    return i;
                }
            }
            points.push(p)
            return points.length-1;
        }
        function getregion(t) {
            for (var i = 0; i < regions.length; i++) {
                if (Vector.eq(regions[i],t)) return i;
            }
            regions.push(t)
            return regions.length-1;
        }
        // We generate points in order so add to list only if greater
        // than the current last element (or if the list is empty).
        function addfacepoint(faces,point,triangle)
        {
            if (!faces[point]) faces[point] = []
            faces[point].push(triangle);
        }
        function sortfaces(pfaces,triangles,q,r)
        {
            function swap(a,i,j) {
                var t = a[i];
                a[i] = a[j];
                a[j] = t;
            }
            pfaces.map(function(regions) {
                if (regions) {
                    // Find first +ve element and move to front
                    for (var i = 0; i < regions.length; i++) {
                        if (triangles[regions[i]][3] == 0) {
                            swap(regions,0,i);
                            break;
                        }
                    }
                    var nextp = q; // Swap q and r to go round other way
                    for (var i = 0; i < regions.length-1; i++) {
                        // Try to swap points[i+1] with a later point that shares two
                        // vertexes with points[i]
                        for (var j = i+2; j < regions.length; j++) {
                            if (triangles[regions[i]][nextp] == triangles[regions[j]][nextp]) {
                                swap(regions,i+1,j);
                                break;
                            }
                        }
                        nextp = q+r-nextp;
                    }
                }
            });
        }

        for (var i = 0; i != regions.length; i++) {
            if (regions.length > 200) {
                alert("Schwarz triangle generation failure");
                return {};
            }
            var t = regions[i];
            var ttype = t[3];
            var px = points[t[0]];
            var qx = points[t[1]];
            var rx = points[t[2]];
            // Face i is part of the p-face associated with p, etc
            addfacepoint(faces[0],t[0],i);
            addfacepoint(faces[1],t[1],i);
            addfacepoint(faces[2],t[2],i);
            var pi = getpoint(Vector.applybary(P,px,qx,rx));
            var qi = getpoint(Vector.applybary(Q,px,qx,rx));
            var ri = getpoint(Vector.applybary(R,px,qx,rx));
            var t1 = getregion([pi,t[1],t[2],1-ttype]);
            var t2 = getregion([t[0],qi,t[2],1-ttype]);
            var t3 = getregion([t[0],t[1],ri,1-ttype]);
            adjacent[i] = [t1,t2,t3];
        }
        sortfaces(faces[0],regions,1,2);
        sortfaces(faces[1],regions,2,0);
        sortfaces(faces[2],regions,0,1);
        this.points = points;
        this.regions = regions;
        this.adjacent = adjacent;
        this.faces = faces;
        this.snuba = getsnub(p,q,r);
    }
    Schwarz.makesymmetry = makesymmetry;
    Schwarz.prototype.describe = function (verbose) {
        var points = this.points;
        var regions = this.regions;
        var adjacent = this.adjacent;
        console.log(points.length + " points; " + regions.length + " regions");
        if (!verbose) return;
        var pfaces = this.faces[0];
        var qfaces = this.faces[1];
        var rfaces = this.faces[2];
        console.log(pfaces.length + " p faces; " + 
                    qfaces.length + " q faces; " +
                    rfaces.length + " r faces");
        for (var i = 0; i < regions.length; i++) {
            console.log(i + ": " + regions[i]);
        }
        for (var i = 0; i < adjacent.length; i++) {
            console.log(adjacent[i]);
        }
        for (var i = 0; i < points.length; i++) {
            console.log(points[i]);
        }
    }
    Schwarz.prototype.tri2bary = function(tri) {
        return Vector.tri2bary(tri,this.p,this.q,this.r);
    }
    Schwarz.prototype.applybary = function(bary,region) {
        var points = this.points;
        var regions = this.regions;
        var p = points[regions[region][0]];
        var q = points[regions[region][1]];
        var r = points[regions[region][2]];
        return Vector.applybary(bary,p,q,r);
    }

    // Compute the distances of the Schwarz triangle vertices from the origin.
    // It's just the dot product with the region point (which is at distance 1).
    // Since we use barycentric coords, we just calculate for region 0.
    Schwarz.prototype.makefacedata = function(bary)
    {
        var points = this.points;
        var regions = this.regions;
        var p = points[0], q = points[1], r = points[2];
        var s = Vector.applybary(bary,p,q,r);
        var s0 = Vector.reflect(q,r,s);
        var s1 = Vector.reflect(r,p,s);
        var s2 = Vector.reflect(p,q,s);
        // Take triple product of both kinds of face triangle.
        // If both are zero, then face is degenerate.
        // FIXME: please explain this
        // TBD: If one is zero, we can omit that sector in drawing face
        var t0 = Math.abs(Vector.triple(p,s,s1)) + Math.abs(Vector.triple(p,s,s2));
        var t1 = Math.abs(Vector.triple(q,s,s0)) + Math.abs(Vector.triple(q,s,s2));
        var t2 = Math.abs(Vector.triple(r,s,s0)) + Math.abs(Vector.triple(r,s,s1));
        var eps = 1e-4;
        // Map of regions to the corresponding points
        var regionpoints = [];
        for (var i = 0; i < regions.length; i++) {
            var region = regions[i];
            regionpoints[i] = Vector.applybary(bary,
                                               points[region[0]],
                                               points[region[1]],
                                               points[region[2]]);
        }
        // Compute barycentric coords of centre of snub triangle.
        // ie. the centre of s0,s1,s2. The centre here means the face
        // normal (so polar reciprocation will work).
        var norm = Vector.cross(Vector.sub(s1,s0),Vector.sub(s2,s0));
        var normlen = Vector.length(norm);
        if (normlen < eps) {
            // If the triangle is degenerate, the calculation fails.
            if (0) {
                // A suitable approximation
                var newbary = [bary[0],bary[1],bary[2]];
                if (Math.abs(newbary[0] < eps)) newbary[0] = eps;
                else if (Math.abs(newbary[1] < eps)) newbary[1] = eps;
                else if (Math.abs(newbary[2] < eps)) newbary[2] = eps;
                var ns = Vector.applybary(newbary,p,q,r);
                var ns0 = Vector.reflect(q,r,ns);
                var ns1 = Vector.reflect(r,p,ns);
                var ns2 = Vector.reflect(p,q,ns);
                norm = Vector.cross(Vector.sub(ns1,ns0),Vector.sub(ns2,ns0));
            } else {
                var tmp;
                // But this is better
                // So, if the snub triangle is degenerate (because the face
                // point is in the corner of the Schwarz triangle), then we want
                // the limit, which turns out to be the mid point of
                // line from the face point in the opposite side.
                if (Math.abs(bary[0] > eps)) tmp = s0;
                else if (Math.abs(bary[1] > eps)) tmp = s1;
                else if (Math.abs(bary[2] > eps)) tmp = s2;
                else tmp = s0;
                norm = Vector.mid(s,tmp);
            }
            normlen = Vector.length(norm);
        }
        norm = Vector.div(norm,normlen);
        // Set the height correctly (to be on the s0,s1,s2 plane)
        norm = Vector.mul(norm,Vector.dot(s0,norm));
        var snubcentre = Vector.getbary(p,q,r,norm);

        // Taking the middle edge seems to give good results
        var snubsphere = Math.max(Vector.dot(Vector.mid(s0,s1)),
                                  Math.min(Vector.dot(Vector.mid(s1,s2)),
                                           Vector.dot(Vector.mid(s2,s0))));
        // Coords of the edge centres of the polyhedron
        var e0 = Vector.mid(s,s0);
        var e1 = Vector.mid(s,s1);
        var e2 = Vector.mid(s,s2);
        var midsphere = Math.min(Vector.dot(e0,e0),
                                 Vector.dot(e1,e1),
                                 Vector.dot(e2,e2));
        var edgecentres = [Vector.getbary(p,q,r,e2),
                           Vector.getbary(p,q,r,e0),
                           Vector.getbary(p,q,r,e1)];
        return {
            facedistances: [Vector.dot(p,s), Vector.dot(q,s), Vector.dot(r,s)],
            faces: [t0 > eps,t1 > eps,t2 > eps], // TBD: think of a better name
            regionpoints: regionpoints,
            snubcentre: snubcentre,
            edgecentres: edgecentres,
            snubsphere: snubsphere,
            midsphere: midsphere,
        };
    }

    // Maybe use a "map over all faces" function
    // TBD: include snub faces as well if required
    Schwarz.prototype.makefacelines = function(facedata,hideface,centre,type,i) {
        //TBD: Also add snub faces
        var points = this.points;
        var regions = this.regions;
        var regionpoints = facedata.regionpoints;
        var facelines = [];
        var stelllimit = 1e4;
        for (var type1 = 0; type1 < 3; type1++) {
            if (!facedata.faces[type1] || hideface[type1]) continue;
            var faces1 = this.faces[type1];
            var rcentre1 = facedata.facedistances[type1];
            for (var j = 0; j < faces1.length; j++) {
                if ((type != type1 || i != j) && faces1[j]) {
                    var p1 = Vector.mul(points[j],rcentre1);
                    var end1 = null, end2 = null, llength = null;
                    for (var type2 = 0; type2 < 3; type2++) {
                        if (!facedata.faces[type2] || hideface[type2]) continue;
                        var faces2 = this.faces[type2];
                        var rcentre2 = facedata.facedistances[type2];
                        for (var k = 0; k < faces2.length; k++) {
                            if ((type2 != type || k != i) &&
                                (type2 != type1 || k != j) &&
                                faces2[k]) {
                                var p2 = Vector.mul(points[k],rcentre2);
                                var s = Vector.intersection(centre,p1,p2);
                                if (s != null && Vector.dist(centre,s) < stelllimit) {
                                    if (!end1) end1 = s;
                                    else if (!end2) {
                                        end2 = s;
                                        llength = Vector.dist(end1,end2);
                                    } else {
                                        // eg: end1 ... end2 ... s
                                        var llen1 = Vector.dist(s,end1);
                                        var llen2 = Vector.dist(s,end2);
                                        if (llen1 > llen2) {
                                            if (llen1 > llength) { 
                                                end2 = s; llength = llen1; 
                                            }
                                        } else if (llen2 > llength) {
                                            end1 = s; llength = llen2;
                                        }
                                    }
                                }
                            }
                        }
                        if (llength) {
                            facelines.push([end1,end2,[i,j,k]]);
                        }
                    }
                }
            }
        }
        return facelines;
    }

    Schwarz.prototype.stellate = function(facedata,hideface) {
        var points = this.points;
        var regions = this.regions;
        facedata.facelines = [];
        for (var type = 0; type < 3; type++) {
            if (facedata.faces[type] && !hideface[type]) {
                var faces = this.faces[type];
                for (var i = 0; i < faces.length; i++) {
                    if (faces[i]) break;  // Find a face
                }
                var rcentre = facedata.facedistances[type];
                var centre = Vector.mul(points[i],rcentre);
                var facelines = this.makefacelines(facedata,hideface,centre,type,i);
                // Now change to barycentric coords
                // faces[i] is the list of regions in the face
                // use coords based on first region
                var p = points[regions[faces[i][0]][0]];
                var q = points[regions[faces[i][0]][1]];
                var r = points[regions[faces[i][0]][2]];
                for (i = 0; i < facelines.length; i++) {
                    var end1 = facelines[i][0];
                    var end2 = facelines[i][1];
                    var bary1 = Vector.getbary(p,q,r,end1);
                    var bary2 = Vector.getbary(p,q,r,end2);
                    facelines[i][0] = bary1;
                    facelines[i][1] = bary2;
                }
                facedata.facelines[type] = facelines;
            }
        }
    }

    // Maybe move these to Vector?
    function applyop(trans,s,P,Q,R) {
        for (var i = 0; i < trans.length; i++) {
            var t = trans[i];
            if (t == 'P') s = Vector.planereflect(P,s);
            else if (t == 'Q') s = Vector.planereflect(Q,s);
            else if (t == 'R') s = Vector.planereflect(R,s);
            else console.log ("Bad operation: " + trans);
        }
        return s;
    }

    function makematrix(trans,P,Q,R) {
        var m0 = applyop(trans,[1,0,0],P,Q,R);
        var m1 = applyop(trans,[0,1,0],P,Q,R);
        var m2 = applyop(trans,[0,0,1],P,Q,R);
        return [ [ m0[0],m1[0],m2[0] ],
                 [ m0[1],m1[1],m2[1] ],
                 [ m0[2],m1[2],m2[2] ] ];
    }

    // Apply single operation trans over 'set' s (an array)
    // Returns new array
    function mapply(trans,s,P,Q,R) {
        var res = new Array(s.length);
        for (var i = 0; i < s.length; i++) {
            res[i] = applyop(trans,s[i],P,Q,R);
        }
        return res;
    }

    function applymatrix(m,u) {
        return [Vector.dot(m[0],u),
                Vector.dot(m[1],u),
                Vector.dot(m[2],u)];
    }

    function test() {
        var triangles = [
            [5,3,2],
            [3,3,5/2],
            [5,5,3/2],
            [5,5/2,2],
            [5,3,5/3],
            [5/2,5/2,5/2],
            [5,3,3/2],
            [5,5,5/4],
            [3,5/2,2],
            [5,5/2,3/2],
            [5,2,5/3],
            [3,5/2,5/3],
            [5,3,5/4],
            [5,2,3/2],
            [3,2,5/3],
            [5/2,5/2,3/2],
            [3,3,5/4],
            [3,5/2,5/4],
            [5/2,2,3/2],
            [5/2,5/3,5/3],
            [3,5/3,3/2],
            [3,2,5/4],
            [5/2,2,5/4],
            [5/2,3/2,3/2],
            [2,5/3,3/2],
            [5/3,5/3,3/2],
            [2,5/3,5/4],
            [2,3/2,5/4],
            [5/3,3/2,5/4],
            [3/2,3/2,5/4],
            [3/2,5/4,5/4],
            [5/4,5/4,5/4],
        ];
        for (var i = 0; i < triangles.length; i++) {
            console.log((i+1) + ":")
            console.log(triangles[i])
            var a = triangles[i][0];
            var b = triangles[i][1];
            var c = triangles[i][2];
            var schwarz = new Schwarz([a,b,c]);
            console.log(schwarz.snuba)
            var schwarz = new Schwarz([a,c,b]);
            console.log(schwarz.snuba)
            var schwarz = new Schwarz([b,a,c]);
            console.log(schwarz.snuba)
            var schwarz = new Schwarz([b,c,a]);
            console.log(schwarz.snuba)
            var schwarz = new Schwarz([c,a,b]);
            console.log(schwarz.snuba)
            var schwarz = new Schwarz([c,b,a]);
            console.log(schwarz.snuba)
        }
    }
    
    // Set up the external symbols.
    Geometry.Vector = Vector;
    Geometry.PointSet = PointSet;
    Geometry.Schwarz = Schwarz;
    Geometry.makematrix = makematrix;
    Geometry.mapply = mapply;
    Geometry.test = test;
}())

// Geometry.test();
