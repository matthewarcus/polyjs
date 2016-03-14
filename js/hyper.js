//var geometry = require("./geometry")

//console.log(geometry)
//console.log(geometry.Geometry)
//var Geometry = geometry.Geometry
(function() {
    var Vector = Geometry.Vector;
    var vprint = Vector.vprint;

    var verbose = false;

    function applybary4(bary,p,q,r,s) {
        var Vector = Geometry.Vector;    
        return Vector.add4(Vector.mul(p,bary[0]),
                           Vector.mul(q,bary[1]),
                           Vector.mul(r,bary[2]),
                           Vector.mul(s,bary[3]));
    }

    function VertexSet(eps) {
        this.vertices = []
        this.eps = eps || 1e-4;
        this.add = function(v) {
            for (var i = 0; i < this.vertices.length; i++) {
                if (Vector.approxeq(v,this.vertices[i],this.eps)) {
                    return i;
                }
            }
            this.vertices.push(v)
            return this.vertices.length-1;
        }
        this.get = function(i) {
            console.assert(i >= 0);
            console.assert(i < this.vertices.length);
            return this.vertices[i];
        }
        this.length = function() {
            return this.vertices.length;
        }
    }

    function qconj(q) {
        return [q[0],-q[1],-q[2],-q[3]];
    }

    function qmul(q0,q1) {
        return [q0[0]*q1[0] - q0[1]*q1[1] - q0[2]*q1[2] - q0[3]*q1[3],
                q0[0]*q1[1] + q0[2]*q1[3] - q0[3]*q1[2] + q0[1]*q1[0],
                q0[0]*q1[2] + q0[3]*q1[1] - q0[1]*q1[3] + q0[2]*q1[0],
                q0[0]*q1[3] + q0[1]*q1[2] - q0[2]*q1[1] + q0[3]*q1[0]]
    }

    var quat1 = Vector.normalize([1,0,0,0]);
    var quat2 = Vector.normalize([1,0,0,0]);
    var dquat1 = Vector.normalize([1,0.005,0.005,0.005]);
    //var dquat2 = Vector.normalize([1,-0.005,-0.005,-0.005]);
    var dquat2 = Vector.normalize([1,0,0,0]);

    //console.log("#",quat1,quat2,Vector.dot(quat1,quat2));

    // Dump set of vertices and edges as an OFF file
    // We do the projection to 3-space here rather than having
    // 4-points in the off file, but it might better to do it that way.
    function makeoff(vertices,edges,wcam) {
        console.assert(wcam >= 0);
        //console.log(wcam);
        var doperspective = true;
        var camera = [0,0,0,-wcam]
        var project = function(v) {
            // Project on to w = 1 from camera at [0,0,0,-wcam], ie.
            // so from the origin, get the gnomonic projection).
            // (camera + k(v-camera)).w = 1
            var k = (1+wcam)/(v[3]+wcam)
            return Vector.add(camera,Vector.mul(Vector.sub(v,camera), k))
        }
        // Handle duplicate vertices, null edges.
        // Perhaps we should remove unused vertices too.
        var newvertices = new VertexSet;
        var vertexmap = [];
        for (var i = 0; i < vertices.length; i++) {
            vertices[i] = qmul(quat1,qmul(vertices[i],quat2));
            //console.log(vertices[i])
            if (!doperspective || vertices[i][3] > -wcam) {
                vertexmap[i] = newvertices.add(vertices[i]);
            }
        }
        var newedges = []
        var checkedge = function(p1,p2,type) {
            if (p1 == p2) return;
            for (var i = 0; i < newedges.length; i++) {
                var edge = newedges[i];
                if (p1 == edge[0] && p2 == edge[1] & type == edge[2]) {
                    return;
                }
            }
            newedges.push([p1,p2,type]);
        }
        var intersect = function(p1,p2) {
            var v1 = vertices[p1];
            var v2 = vertices[p2];
            // Find k such that (v1 + k(v2-v1)).w = -0.9*wcam
            var k = (v1[3]+0.99*wcam)/(v1[3]-v2[3]);
            var newv = Vector.add(v1, Vector.mul(Vector.sub(v2,v1),k));
            //console.log(newv,0.9*D);
            return newv;
        }
        for (var i = 0; i < edges.length; i++) {
            var p1 = vertexmap[edges[i][0]]
            var p2 = vertexmap[edges[i][1]]
            if (p1 != undefined && p2 == undefined) {
                p2 = newvertices.add(intersect(edges[i][0], edges[i][1]))
            } else if (p2 != undefined && p1 == undefined) {
                p1 = newvertices.add(intersect(edges[i][0], edges[i][1]))
            }
            if (p1 != undefined && p2 != undefined) {
                var type = edges[i][2]
                checkedge(p1,p2,type)
            }
        }
        if (0) {
            console.log("OFF");
            console.log("" + newvertices.length() + " " + newedges.length + " " + 0);
            for (var i = 0; i < newvertices.length(); i++) {
                var v = project(newvertices.get(i));
                v = v.map(function(x) {
                    return x.toPrecision(6);
                })
                console.log(v[0] + " " + v[1] + " " + v[2])
            }
            for (var i = 0; i < newedges.length; i++) {
                if (verbose) {
                    console.log("#", Vector.length(Vector.sub(newvertices.get(newedges[i][0]),
                                                              newvertices.get(newedges[i][1]))));
                }
                console.log(2 + " " + newedges[i][0] + " " + newedges[i][1] + " " + newedges[i][2]);
            }
        } else {
            var vertices = newvertices.vertices.map(function(v) {
                if (doperspective) {
                    var v1 = project(v);
                } else {
                    var v1 = v;
                }
                return new THREE.Vector3(v1[0],v1[1],v1[2]);
            });
            var faces = newedges.map(function(edge) {
                return { vlist: [edge[0],edge[1]], type: edge[2] }
            })
            return { vertices: vertices, faces: faces }
        }
    }

    // x,y,z,w are normals
    // Now form closure of reflecting among themselves.
    function close(p,q,r,s,P,Q,R,S,C,camera) {
        var Vector = Geometry.Vector;
        var dot = Vector.dot;
        var vprint = Vector.vprint;
        var vertices = [p,q,r,s]
        var vmap = [[],[],[],[]];
        function getvertex(p) {
            for (var i = 0; i < vertices.length; i++) {
                if (Vector.approxeq(vertices[i],p,1e-4)) {
                    return i;
                }
            }
            vertices.push(p);
            return vertices.length-1;
        }
        for (var i = 0; i < vertices.length; i++) {
            var v = vertices[i];
            var pi = getvertex(Vector.planereflect(P,v)); // NOTE ORDER!!
            var qi = getvertex(Vector.planereflect(Q,v));
            var ri = getvertex(Vector.planereflect(R,v));
            var si = getvertex(Vector.planereflect(S,v));
            vmap[0][i] = pi;
            vmap[1][i] = qi;
            vmap[2][i] = ri;
            vmap[3][i] = si;
        }
        var region0 = [0,1,2,3]
        var rmap = [[region0]]
        var regions = [region0]
        function check(region,k) {
            // Reflect region in plane[k] 
            var region2 = [vmap[k][region[0]],
                           vmap[k][region[1]],
                           vmap[k][region[2]],
                           vmap[k][region[3]]];
            var j = region2[0]
            if(!rmap[j]) rmap[j] = [];
            for (var i = 0; i < rmap[j].length; i++) {
                if (Vector.eq(region2,rmap[j][i])) return;
            }
            rmap[j].push(region2);
            regions.push(region2);
        }
        for (var i = 0; i < regions.length; i++) {
            var region = regions[i];
            check(region,0)
            check(region,1)
            check(region,2)
            check(region,3)
        }
        //console.log("# regions:", regions.length)
        //console.log("# vertices:", vertices.length)
        if (verbose) {
            for (var i = 0; i < vertices.length; i++) {
                console.log(i, vprint(vertices[i]));
            }
            for (var i = 0; i < regions.length; i++) {
                console.log(regions[i]);
            }
        }
        // Now, we can calculate edges. We have an edge
        // between all region pairs that differ in just one
        // digit eg. ([a,b,c,d],[a,b,c',d]), eg. the first
        // such pair, [ 0, 1, 2, 3 ] and [ 0, 1, 2, 7 ].

        // And sets of regions that differ in just two digits
        // form a face, eg:
        // [ 0, 1, 2, 3 ]
        // [ 0, 1, 2, 7 ]
        // [ 0, 1, 12, 7 ]
        // [ 0, 1, 12, 13 ]
        // [ 0, 1, 6, 13 ]
        // [ 0, 1, 6, 3 ]
        // Here we can order the regions so that they correspond
        // to a series of alternate reflections.
        // And sets that differ in 3 digits form a single cell, etc.
        // For now, let's just generate the edges:
        var edges = []
        for (var i = 0; i < regions.length; i++) {
            for (var j = i+1; j < regions.length; j++) {
                var count = 0;
                var type = 0;
                for (var k = 0; k < regions[i].length; k++) {
                    if (regions[i][k] != regions[j][k]) {
                        count++;
                        type = k;
                    }
                }
                if (count == 1) edges.push([i,j,type])
            }
        }
        if (verbose) {
            console.log("edges",edges.length);
            for (var i = 0; i < edges.length; i++) {
                console.log("e",regions[edges[i][0]],regions[edges[i][1]]);
            }
        }
        // This should be the end of the function.
        // return { vertices: vertices, regions: regions, edges: edges }

        // It's a relief that the number of edges is twice the number of
        // regions - just what we would expect.
        // Now let's try and get just the edges of hypercube (we all know
        // what one of them looks like, so it's a good test. We need to pick
        // a vertex of the fundamental region, but which one? According to
        // the Dynkin diagram (Coxeter-Dynkin is a bit long & Coxeter has plenty
        // of other things named after him (deservedly so of course) whereas
        // Dynkin doesn't seem to have much else [this is totally
        // untrue, Dynkin is a figure of some interest - see
        // Wikipedia], it should be the first vertex of the region,
        // anyway (the ringed node in the diagram).
        //makeoff(vertices,hypercube);
        
        // Also, want "barycentric" coords:
        // If p,q,r,s are vertices of tetrahedron, then want:
        // a,b,c,d such that:
        // ap + bq + cr + ds = v
        // Once again, this is just a matrix inverse
        // A[a,b,c,d] = v

        // Now I want barycentric coordinates for C:
        // Express C as a weighted sum of the vertices p,q,r,s
        // C = ap + bq + cr + ds
        // ie. C = M [a,b,c,d] where
        var M = [ p[0], q[0], r[0], s[0],
                  p[1], q[1], r[1], s[1],
                  p[2], q[2], r[2], s[2],
                  p[3], q[3], r[3], s[3] ];
        // Once again, use the magic of linear algebra:
        var N = Geometry.invert4(M);
        var bary = Geometry.apply4(N,C)
        if (verbose) {
            console.log("C", vprint(C));
            console.log("Check", vprint([dot(C,P),dot(C,Q),dot(C,R),dot(C,S)]));
            console.log("bary",vprint(bary));
            console.log(vprint(Geometry.apply4(M,bary)));
            console.log(applybary4(bary,p,q,r,s));
        }
        var regionpoints = regions.map(function(region) {
            return applybary4(bary,
                              vertices[region[0]],
                              vertices[region[1]],
                              vertices[region[2]],
                              vertices[region[3]]);
        });
        // Check each region has a unique point
        if (0) {
            for (var i = 0; i < regionpoints.length; i++) {
                for (var j = i+1; j < regionpoints.length; j++) {
                    if (Vector.approxeq(regionpoints[i],
                                        regionpoints[j],
                                        1e-4)) {
                        console.log("# Error: shared regionpoints",i,j);
                    }
                }
            }
        }
        return makeoff(regionpoints,edges,camera);
    }

    // Angles are P.Q, Q.R, R.P, R.S, P.S, Q.S making Schwarz triangles.
    // [A,B,C,a,b,c] form Schwarz triangles [A,B,C],[A,B,c],[A,b,C],[a,B,C]
    // representing the 4 elements of the symmetry group.

    // Given 6 dihedral angles (expressed as reciprocal multiples of pi)
    // construct 4 planes with those angles. Return null if no solution
    // found.
    function solveangles(angles,verbose) {
        var Vector = Geometry.Vector
        var dot = Vector.dot
        var cross = Vector.cross
        function dihedral(x) { return Math.cos(Math.PI*(1-1/x)); }
        var A = dihedral(angles[0])
        var B = dihedral(angles[1])
        var C = dihedral(angles[2])
        var D = dihedral(angles[3])
        var E = dihedral(angles[4])
        var F = dihedral(angles[5])
        var P = [1,0,0,0]
        // Q = [a,b,0,0], P.Q = A
        var a = A;
        if (a*a >= 1) {
            if (verbose) console.log("No solution for Q")
            return null;
        }
        var b = Math.sqrt(1 - a*a)
        //console.log(a*a,b)
        var Q = [a,b,0,0]
        // R = [c,d,e,0]
        // R.P = C = c
        // Q.R = B = ac + bd
        var c = C, d = (B-a*c)/b
        if (c*c + d*d >= 1) {
            if (verbose) console.log("# No solution for R")
            return null;
        }
        var e = Math.sqrt(1 - c*c - d*d)
        var R = [c,d,e,0]
        if (verbose) {
            console.log("# P", vprint(P))
            console.log("# Q", vprint(Q))
            console.log("# R", vprint(R))
            console.log("# P.Q",dot(P,Q),A)
            console.log("# Q.R",dot(Q,R),B)
            console.log("# R.P",dot(R,P),C)
        }

        var X = cross(P,Q)
        var Y = cross(Q,R)
        var Z = cross(R,P)

        // u = iX + jY + kZ
        // R.S = R.u = D = iR.X
        // P.S = P.u = E = jP.Y
        // Q.S = Q.u = F = kQ.Z

        var i = D/dot(R,X)
        var j = E/dot(P,Y)
        var k = F/dot(Q,Z)

        var u = Vector.add3(Vector.mul(X,i),Vector.mul(Y,j),Vector.mul(Z,k))
        if (verbose) {
            console.log("#", i,j,k)
            console.log("#", u)
        }
        // We might have dot(u,u) = 1, in which case the fundamental region
        // is 3 dimensional and we have an infinite 3d honeycomb. We can't
        // do honeycombs (yet) so exclude this case.
        if (dot(u,u) >= 1-1e-4) {
            if (verbose) console.log("No solution for S", dot(u,u))
            return null;
        }
        var k = Math.sqrt(1 - dot(u,u))
        var S = Vector.add(u,Vector.mul([0,0,0,1],k))
        if (verbose) {
            console.log("# S", vprint(S))
            console.log("# R.S",dot(R,S),D)
            console.log("# P.S",dot(P,S),E)
            console.log("# Q.S",dot(Q,S),F)
        }
        return [P,Q,R,S];
    }

    function polychoron(angles,quad,camera) {
        var Vector = Geometry.Vector
        var dot = Vector.dot
        var planes = solveangles(angles,verbose);
        var P = planes[0], Q = planes[1], R = planes[2], S = planes[3];
        // Now solve the trilinear (quadriplanar) equations, eg:
        // P.p = a
        // Q.p = b
        // R.p = c
        // S.p = d
        // Construct 4*4 matrix and invert, then p = (x,y,z,w) = Pinv(a,b,c,d)
        // The 4 face centres are then:

        // Norm(Pinv([1,0,0,0]))
        // Norm(Pinv([0,1,0,0]))
        // Norm(Pinv([0,0,1,0]))
        // Norm(Pinv([0,0,0,1]))

        // And we can find the face points as eg:
        // Norm(Pinv([1,1,0,0]))
        // Norm(Pinv([1,1,1,1]))
        // etc.

        var m = [
            P[0],P[1],P[2],P[3],
            Q[0],Q[1],Q[2],Q[3],
            R[0],R[1],R[2],R[3],
            S[0],S[1],S[2],S[3]];
        var n = Geometry.invert4(m)

        // Now n magically solves our trilinear (or rather, following Coxeter, quadriplanar) coordinates:
        // For honeycombs (ie. with a 3d region), m isn't invertible, alas.
        var p = Vector.normalize(Geometry.apply4(n,[1,0,0,0]))
        var q = Vector.normalize(Geometry.apply4(n,[0,1,0,0]))
        var r = Vector.normalize(Geometry.apply4(n,[0,0,1,0]))
        var s = Vector.normalize(Geometry.apply4(n,[0,0,0,1]))
        if (verbose) {
            console.log(vprint(m));
            console.log(vprint(n));
            console.log("p", vprint(p), vprint([dot(p,P),dot(p,Q),dot(p,R),dot(p,S)]));
            console.log("q", vprint(q), vprint([dot(q,P),dot(q,Q),dot(q,R),dot(q,S)]));
            console.log("r", vprint(r), vprint([dot(r,P),dot(r,Q),dot(r,R),dot(r,S)]));
            console.log("s", vprint(s), vprint([dot(s,P),dot(s,Q),dot(s,R),dot(s,S)]));
        }
        var regionpoint = /*Vector.normalize*/(Geometry.apply4(n,quad));
        return close(p,q,r,s,P,Q,R,S,regionpoint,camera);
    }

    // Find all integral solutions for angles, up to some limit..
    // First three numbers are a Schwarz triangle and all permutations
    // are valid, so generate in lexical order.
    // Last three numbers aren't symmetric (though we can permute if
    // corresponding numbers in first triple are duplicated).
    function solveall(N) {
        for (var i = 2;  i < N; i++) 
            for (var j = i;  j < N; j++) 
                for (var k = j;  k < N; k++) 
                    for (var l = 2;  l < N; l++) 
                        for (var m = 2;  m < N; m++) 
                            for (var n = 2;  n < N; n++) {
                                if (i == j && l > m) continue;
                                if (j == k && m > n) continue;
                                var p = solveangles([i,j,k,l,m,n]);
                                if (p) {
                                    console.log(i,j,k,l,m,n);
                                }
                            }
    }

    PolyContext.prototype.hyper = function(off,options,running) {
        var angles = [];
        var quad = [];
        options.w = Number(options.w) || 4;
        options.angles = options.angles || "4:3:2:3:2:2"
        options.quad = options.quad || "1:0:0:0"
        var matches;
        if (matches = options.angles.match(/^([\d]+)(?:\/([\d]+))?:([\d]+)(?:\/([\d]+))?:([\d]+)(?:\/([\d]+))?:([\d]+)(?:\/([\d]+))?:([\d]+)(?:\/([\d]+))?:([\d]+)(?:\/([\d]+))?$/)) {
            //console.log("Got angles",matches);
            for (var i = 0; i < 6; i++) {
                var a = Number(matches[2*i+1]);
                if (matches[2*i+2]) a /= Number(matches[2*i+2]);
                angles.push(a);
            }
        } else {
            alert("Invalid off.angles: " + options.angles);
            return null;
        }
        if (matches = options.quad.match(/^([\d]+):([\d]+):([\d]+):([\d]+)$/)) {
            //console.log("Got quad",matches);
            for (var i = 0; i < 4; i++) {
                var a = Number(matches[i+1]);
                quad.push(a);
            }
        } else {
            alert("Invalid quad: " + options.quad);
            return null;
        }
        //console.log(options.w);
        //var angles = ;
        //var angles = [5,3,2,3,2,2];
        //var angles = [5/2,3,2,3,2,2];
        //var angles = [3,4,2,3,2,2]
        //var quad = [1,0,0,0];
        //var quad = [0,1,1,0];
        //var quad = [1,0,0,0];
        var camera = options.w;
        if (running) {
            quat1 = qmul(quat1,dquat1);
            quat2 = qmul(quat2,dquat2);
        }
        return polychoron(angles,quad,camera);
    }
    PolyContext.prototype.quat = function(off,options,running) {
        let quat = [2,0,0,0];
        var dquat1 = Vector.normalize([1,0.2,0.0,0.0]);
        var dquat2 = Vector.normalize([1,0.0,0.2,0.0]);
        var dquat3 = Vector.normalize([1,0.0,0.0,0.2]);
        var dquat4 = Vector.normalize([1,0.2,0.2,0.2]);
        var vertices = [];
        var faces = [];
        var N = 1000;
        for (var i = 0; i < N; i++) {
            vertices.push(new THREE.Vector3(quat[1],quat[2],quat[3]));
            quat = qmul(quat,dquat1);
            quat = qmul(dquat2,quat);
            quat = qmul(quat,dquat4);
            //quat = qmul(dquat4,quat);
        }
        for (var i = 0; i < N-1; i++) {
            faces.push({ vlist: [i,i+1] });
        }
        return {vertices: vertices, faces: faces }
    }
})()

//module.exports = { solveall: solveall, polychoron: polychoron }

