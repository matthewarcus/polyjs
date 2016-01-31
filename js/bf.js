PolyContext.prototype.loadbf = function(off,options) {
    var vector = THREE.OFFLoader.Utils.vector
    var vnegate = THREE.OFFLoader.Utils.vnegate
    var vadd = THREE.OFFLoader.Utils.vadd
    var vadd3 = THREE.OFFLoader.Utils.vadd3
    var vsub = THREE.OFFLoader.Utils.vsub
    var vdiv = THREE.OFFLoader.Utils.vdiv
    var vmul = THREE.OFFLoader.Utils.vmul
    var vdot = THREE.OFFLoader.Utils.vdot
    var vcross = THREE.OFFLoader.Utils.vcross
    var Color = THREE.OFFLoader.Utils.Color;
    if (options.bfheight === undefined) {
        options.bfheight = 0.0;
        options.bfheightinc = 0.02;
    }
    function yrot(p,theta) {
        return vector(Math.cos(theta)*p.x + Math.sin(theta)*p.z,
                      p.y,
                      -Math.sin(theta)*p.x + Math.cos(theta)*p.z)
    }
    function solve(a,b,j,k,r) {
        //console.log("Solve",a,b,j,k,r)
        // Solve v.a = j, v.b = k, |v| = r
        var c = vcross(a,b)
        // Find A,B,C such that v = Aa + Bb + Cc
        // v.a = A*aa + B*ab = j
        // v.b = A*ab + B*bb = k
        // A = (j - B*ab)/aa
        // (j - B*ab)*ab/aa + B*bb = k
        // (j - B*ab)*ab + B*bb*aa = k*aa
        // j*ab - B*ab*ab + B*bb*aa = k*aa
        // B(bb*aa - ab*ab) = k*aa - j*ab 
        // A = (j - B*ab)/aa
        //console.log(c)
        var aa = vdot(a,a)
        var bb = vdot(b,b)
        var cc = vdot(c,c)
        var ab = vdot(a,b)
        //console.log("# Params1:", aa,bb,cc,ab)
        var B = (k*aa-j*ab)/(bb*aa - ab*ab)
        var A = (j - B*ab)/aa
        // v.v = aa + 2ab + bb + cc (since c orthogonal to a and b)
        // C2 = 0 just when |A*a + B*b| = r
        var C2 = (r*r - A*A*aa - 2*A*B*ab - B*B*bb)/cc;
        var v = vadd(vmul(a,A),vmul(b,B));
        //console.log("# Params2:", A,B,C2);
        if (C2 < 0) {
            // Clamp?
            console.log("# impossible length");
            return;
        }
        var C = Math.sqrt(C2);
        var v1 = vadd(vadd(vmul(a,A),vmul(b,B)),vmul(c,C));
        var v2 = vsub(vadd(vmul(a,A),vmul(b,B)),vmul(c,C));
        //console.log("# Result:", v,vdot(v,a),j,vdot(v,b),k,vlength(v),r);
        //console.log("# Result: ",C2)
        return [v1,v2];
    }

    for (var attempt = 0; attempt < 3; attempt++) {
        var h = options.bfheight
        options.bfheight += options.bfheightinc;
        var vertices = []
        var faces = []
        function addvertex(v,color) {
            var i = vertices.length
            vertices.push(v);
            if (color) faces.push({ vlist: [i], color: color })
            return i
        }
        function addface(vlist,color) {
            faces.push({ vlist: vlist, color: color });
        }
        if (Math.abs(h) < 1e-5) h = 1e-5; // A bit of a cop out
        var A = vector(0,h,0)
        var B = vector(2,h,0)
        var j = -2*vdot(A,B)/3
        var k = -4
        var r = 2
        // Need v.A = j, v.B = k, |v| = r
        var s = solve(A,B,j,k,r)
        //console.log(s[0],s[1])
        if (s == undefined) {
            // Reverse direction and carry on
            options.bfheightinc *= -1
            options.bfheight += 2*options.bfheightinc
        } else {
            var C0 = s[0]
            var C1 = s[1]
            var C = C1
            //console.log("# Check 1: " + vdot(B,C))
            //console.log("# Check 2: " + vdot(vadd(B,C), C))
            addvertex(vector(0,0,0),Color.red);
            var v1 = vadd(B,C);
            var v2  = B;
            var v3 = vadd(B,vmul(C,1.5));
            var v4 = vcross(v1,C)
            var k = Math.sqrt(vdot(C,C)*0.75/vdot(v4,v4))
            var v5 = vadd3(B,vmul(C,1.5),vmul(v4,k));
            var v6 = vadd3(B,vmul(C,1.5),vmul(v4,-k));
            var PI = 3.14159
            var upper = []
            var lower = []
            var N = 3
            for (var i = 0; i < N; i++) {
                var rotation = 2*PI*i/N;
                var n1 = addvertex(yrot(v1,rotation),Color.blue)
                var n2 = addvertex(yrot(v2,rotation),Color.red);
                var n3 = addvertex(yrot(v3,rotation),Color.green);
                var n7 = addvertex(yrot(v5,rotation),Color.cyan);
                var n8 = addvertex(yrot(v6,rotation),Color.cyan);
                addface([0,n1],Color.cyan);
                addface([n2,n3],Color.white);
                addface([n2,n7],Color.white);
                addface([n7,n8],Color.white);
                addface([n8,n2],Color.white);
                addface([n2,n7,n8],Color.orange)
                lower.push(n2);
                var n4 = addvertex(yrot(vnegate(v1),rotation),Color.blue)
                var n5 = addvertex(yrot(vnegate(v2),rotation),Color.red);
                var n6 = addvertex(yrot(vnegate(v3),rotation),Color.green);
                var n9 = addvertex(yrot(vnegate(v5),rotation), Color.cyan);
                var n10 = addvertex(yrot(vnegate(v6),rotation),Color.cyan);
                addface([0,n4],Color.cyan);
                addface([n5,n6],Color.white);
                addface([n5,n9],Color.white);
                addface([n9,n10],Color.white);
                addface([n10,n5],Color.white);
                addface([n5,n9,n10],Color.orange)
                upper.push(n5);
            }
            addface(upper,Color.yellow);
            addface(lower,Color.yellow);
            for (var i = 0; i < N; i++) {
                addface([lower[i],lower[(i+1)%N]],Color.white);
                addface([upper[i],upper[(i+1)%N]],Color.white);
            }
            return {vertices: vertices, faces: faces}
        }
    }
    console.log("loadbf failed", options.bfheight, options.bfheightinc)
    return {vertices: [], faces: [] }
}

PolyContext.prototype.desargues = function(off,options) {
    var vector = THREE.OFFLoader.Utils.vector
    var vadd = THREE.OFFLoader.Utils.vadd
    var vsub = THREE.OFFLoader.Utils.vsub
    var vdiv = THREE.OFFLoader.Utils.vdiv
    var vmul = THREE.OFFLoader.Utils.vmul
    var vdot = THREE.OFFLoader.Utils.vdot
    var vcross = THREE.OFFLoader.Utils.vcross
    var vnormalize = THREE.OFFLoader.Utils.vnormalize
    var Color = THREE.OFFLoader.Utils.Color;
    var vertices = []
    var faces = []
    function addvertex(v,color) {
        var i = vertices.length
        vertices.push(v);
        if (color) faces.push({ vlist: [i], color: color })
        return i
    }
    function addface(vlist,color) {
        faces.push({ vlist: vlist, color: color });
    }
    function addline(p,q,r) {
        addface([p,q],Color.white)
        addface([q,r],Color.white)
        addface([r,p],Color.white)
    }
    function intersect(p0,p1,q0,q1) {
        var n = vsub(p1,p0)
        var m = vsub(q1,q0)
        // Want p0 + kn = q0 + km
        // or
        // kn = q0 - p0 + km
        // k(n x m) = (q0 - p0) x m
        var t1 = vcross(n,m)
        var t2 = vcross(vsub(q0,p0),m)
        console.log(vcross(t1,t2).length())
        var k = vdot(t1,t2)/vdot(t1,t1)
        return vadd(p0,vmul(n,k))
    }
    var a0 = vector(-2,0,-2)
    var b0 = vector(0,-1,-2)
    var c0 = vector(-2,-2,-2)
    var k0 = 1.333;
    var k1 = 1.618;
    var k2 = 2;
    var a1 = vmul(a0,k0)
    var b1 = vmul(b0,k1)
    var c1 = vmul(c0,k2)
    var d = intersect(a0,b0,a1,b1)
    var e = intersect(b0,c0,b1,c1)
    var f = intersect(c0,a0,c1,a1)
    var O = addvertex(vector(0,0,0),Color.blue)
    var A0 = addvertex(a0,Color.red);
    var B0 = addvertex(b0,Color.red);
    var C0 = addvertex(c0,Color.red);
    var A1 = addvertex(a1,Color.yellow);
    var B1 = addvertex(b1,Color.yellow);
    var C1 = addvertex(c1,Color.yellow);
    var D = addvertex(d,Color.green);
    var E = addvertex(e,Color.green);
    var F = addvertex(f,Color.green);
    addline(O,A0,A1)
    addline(O,B0,B1)
    addline(O,C0,C1)
    addline(A0,B0,D)
    addline(A1,B1,D)
    addline(B0,C0,E)
    addline(B1,C1,E)
    addline(C0,A0,F)
    addline(C1,A1,F)
    addline(D,E,F)
    return { vertices: vertices, faces: faces }
}

PolyContext.prototype.theorem = function(off,options) {
    var vector = THREE.OFFLoader.Utils.vector
    var vadd = THREE.OFFLoader.Utils.vadd
    var vdiv = THREE.OFFLoader.Utils.vdiv
    var Color = THREE.OFFLoader.Utils.Color;
    var cos = Math.cos
    var sin = Math.sin
   
    var vertices = []
    var faces = []
    function addvertex(v,color) {
        var i = vertices.length
        vertices.push(v);
        if (color) faces.push({ vlist: [i], color: color })
        return i
    }
    function addface(vlist,color) {
        faces.push({ vlist: vlist, color: color });
    }
    function bisect(v0,v1,color) {
        return addvertex(vdiv(vadd(vertices[v0],vertices[v1]),2),color)
    }
    if (options.theta === undefined) options.theta = 0
    var theta = options.theta
    options.theta += 0.01
    var pi = 3.14159;
    var A = addvertex(vector(1,-0.5,1), Color.green)
    var B = addvertex(vector(-0.5, 0, -0.5), Color.green)
    var C = addvertex(vector(-1, -0.5, 1), Color.green)
    var D = addvertex(vector(sin(theta),1,-cos(theta)), Color.green)

    var E = bisect(A,C,Color.yellow)
    var F = bisect(C,B,Color.yellow)
    var G = bisect(A,D,Color.yellow)
    var H = bisect(D,B,Color.yellow)
    var I = bisect(G,F,Color.red)

    var edges = [
        [A,B],[B,C],[C,A],[A,D],[B,D],
        [G,H],[H,F],[F,E],[E,G],
        [G,F],[H,E]
    ]
    edges.forEach(function(vlist) { addface(vlist,Color.white); });
    return { vertices: vertices, faces: faces }
}

PolyContext.prototype.zonohedron = function(off,options) {
    var vcross = THREE.OFFLoader.Utils.vcross;
    var vtaxi = THREE.OFFLoader.Utils.vtaxi;
    var vntaxi = THREE.OFFLoader.Utils.vntaxi;
    var star;
    if (off) {
        //console.log("Making star from poly")
        star = [];
        // Make a star from the vertices in off
        off.vertices.forEach(function (v) {
            // This linear scan will do for now.
            // Check for equality first
            for (var i = 0; i < star.length; i++) {
                if (vtaxi(star[i],v) < 1e-5 || vntaxi(star[i],v) < 1e-5) {
                    return;
                }
            }
            // Then for collinearity
            for (var i = 0; i < star.length; i++) {
                if (vcross(star[i],v).length() < 1e-5) {
                    return;
                }
            }
            star.push(v);
        })
    } else {
        var N = options.n || 3;
        var M = options.m || 1;
        var scale = 1;
        if (options.theta != undefined) {
            console.log(options.theta)
            scale = 2*Math.sin(options.theta);
            options.theta += 0.01;
        }
        star = THREE.OFFLoader.makeStar(N,M,scale)
    }
    var k = options.k || 0
    var newstar = []
    for (var i = 0; i < k; i++) {
        THREE.OFFLoader.starZonohedron(star,null,newstar)
        star = newstar; newstar = []
    }
    var off = { vertices: [], faces: [] }
    off = THREE.OFFLoader.starZonohedron(star,off)
    if (0) {
        for (var i = 0; i < off.vertices.length; i++) {
            console.log(off.vertices[i])
        }
        for (var i = 0; i < off.faces.length; i++) {
            console.log(off.faces[i].vlist)
        }
    }
    return off
}

// Draw a semi-regular polygon.
// Assumes first 2 edges are typical of all and that the edge centre
// is closest point to origin (both are the case with Wythoff polyhedra). 
// Polygon sides are 2pi*n/m-theta, 2pi*n/m+theta, theta may be
// negative to give a retroflex side which should be drawn as a sort
// of ear shape.
PolyContext.prototype.polygon = function(off,options) {
    var vector = THREE.OFFLoader.Utils.vector
    var vadd = THREE.OFFLoader.Utils.vadd
    var vsub = THREE.OFFLoader.Utils.vsub
    var vdiv = THREE.OFFLoader.Utils.vdiv
    var vmul = THREE.OFFLoader.Utils.vmul
    var vdot = THREE.OFFLoader.Utils.vdot
    var Color = THREE.OFFLoader.Utils.Color;
    var cos = Math.cos
    var sin = Math.sin
   
    var vertices = []
    var faces = []
    function vertex(i) {
        while (i < 0) i += 2*n
        while (i >= 2*n) i -= 2*n
        return vertices[i]
    }
    function addvertex(v,color) {
        var i = vertices.length
        vertices.push(v);
        if (color) faces.push({ vlist: [i], color: color })
        return i
    }
    function addface(vlist,color) {
        faces.push({ vlist: vlist, color: color });
    }
    var theta = options.theta || 0
    options.theta = theta + 0.01
    var n = options.n || 5
    var m = options.m || 1
    //console.log("polygon",n,m,theta)
    for (var i = 0; i < n; i++) {
        var a = i*Math.PI*2*m/n;
        addvertex(vector(sin(a-theta),cos(a-theta),0),Color.red)
        addvertex(vector(sin(a+theta),cos(a+theta),0),Color.green)
    }
    var origin = addvertex(vector(0,0,0))
    // Compute intersection of adjoining edges
    // If edge side of centre, use that (might be
    // further out from edge or other side of centre,
    // either of which is OK.
    var efacts = [], eindex;
    for (var i = 0; i < 2; i++) {
        var p0 = vertex(i-1) // previous vertex
        var p1 = vertex(i)
        var p2 = vertex(i+1) // Edge vertices
        var c = vdiv(vadd(p1,p2),2) // Edge centre
        var e = vsub(p2,p1) // Edge vector
        var v = vsub(p1,p0) // Previous edge
        var k = -vdot(p0,e)/vdot(v,e) 
        var q = vadd(p0,vmul(v,k)) // Intersection with centre line
        // q = kc => q.c = k*c.c => k = q.c/c.c
        var efact = vdot(q,c)/vdot(c,c)
        efacts.push(efact)
        if (0 < efact && efact < 1) eindex = i; // Need an ear
    }
    // New list of vertices
    var vlist = []
    for (var i = 0; i < 2*n; i++) {
        var p1 = vertex(i)
        var p2 = vertex(i+1)
        var c = vdiv(vadd(p1,p2),2)
        var q = vmul(c,efacts[i%2])
        addvertex(c,Color.blue)
        var qi = addvertex(q,Color.cyan)
        addface([i,(i+1)%(2*n)],Color.yellow)
        if (i%2 == eindex) vlist.push(qi)
    }
    console.assert(vlist.length == 0 || vlist.length == n)
    if (vlist.length > 0) {
        // Draw the ear as a separate triangle
        for (var i = 0; i < n; i++) {
            var p0 = (2*i+eindex)%(2*n)
            var p1 = (2*i+eindex+1)%(2*n)
            var q0 = vlist[i]
            var q1 = vlist[(i+1)%n]
            addface([q0,p0,p1],Color.white)
            addface([origin,q0,q1],Color.white)
        }
    } else {
        // Easy case, draw straight out from centre
        for (var i = 0; i < 2*n; i++) {
            var p0 = i
            var p1 = (i+1)%(2*n)
            addface([origin,p0,p1],Color.white)
        }
    }
    return { vertices: vertices, faces: faces }
}

PolyContext.prototype.dipolygonid = function(off,options) {
    if (off) {
        off = THREE.OFFLoader.dipolygonid(off)
        THREE.OFFLoader.twistertransform(off,options)
    }
    return off
}

PolyContext.prototype.hoberman = function(off,options) {
    if (!off) return
    var vadd = THREE.OFFLoader.Utils.vadd
    var vsub = THREE.OFFLoader.Utils.vsub
    var vdiv = THREE.OFFLoader.Utils.vdiv
    var vmul = THREE.OFFLoader.Utils.vmul
    var vdot = THREE.OFFLoader.Utils.vdot
    var vcross = THREE.OFFLoader.Utils.vcross
    var vdist = THREE.OFFLoader.Utils.vdist
    var Color = THREE.OFFLoader.Utils.Color
    var vertices = off.vertices
    var nvertices = vertices.length
    var faces = off.faces
    var newvertices = []
    var newfaces = []
    var lambda = options.lambda || 1
    options.lambdainc = options.lambdainc || 0.01
    options.lambda = lambda + options.lambdainc
    var eps = 1e-4
    if (1-eps < lambda && lambda < 1+eps) { lambda = 1-eps; }
    if (-1-eps < lambda && lambda < -1+eps) { lambda = -1+eps }
    function addedge(p1,p2,color) {
        newfaces.push({ vlist: [p1,p2], color: color })
    }
    function addvertex(p,color) {
        //console.log("Vertex",p,color)
        var i = newvertices.length
        newvertices.push(p)
        if (color) {
            newfaces.push({ vlist: [i], color: color })
        }
        return i;
    }
    // Edge linkage transformation
    var lambdavertices;
    var thetavertices;
    var kappa;
    for (var faceindex = 0; faceindex < faces.length; faceindex++) {
        var face = faces[faceindex]
        var vlist = face.vlist
        if (vlist.length == 2) {
            var color = face.color
            var p = vertices[vlist[0]];
            var q = vertices[vlist[1]];
            var e = vdiv(vadd(p,q),2)    //edge vector
            if (lambdavertices === undefined) {
                var len = vdist(p,q)
                var A = vdot(q,q)
                var B = -2*lambda*vdot(p,q)
                var C = lambda*lambda*vdot(p,p) - len*len
                // TBD: check on stable quadratic solutions
                var disc = B*B - 4*A*C
                //console.log("solution to edge link",lambda,disc,A,B,C)
                if (disc < 0) {
                    options.lambdainc = -options.lambdainc
                    options.lambda += options.lambdainc
                    options.lambda += options.lambdainc
                    return;
                }
                // Either + or - is realizable of course.
                var theta = (-B + Math.sqrt(disc))/(2*A)
                var p1 = vmul(p,lambda)
                var p2 = vmul(p,theta)
                var q1 = vmul(q,lambda)
                var q2 = vmul(q,theta)
                var en = vdiv(vadd(p1,q2),2) // new edge centre
                var n = vdiv(vsub(p1,q2),len) // normal to bisecting plane
                kappa = vdot(en,n)/vdot(e,n) // distance to plane
                //console.log("Plane ", n, n.length(), kappa, vdot(en,n), vdot(e,n))
                //console.log(p,q,len)
                //console.log(p1,q2,vdist(p1,q2))
                //console.log(q1,p2,vdist(q1,p2))
                //showedge(new THREE.Vector3(0,0,0),p);
                //showedge(new THREE.Vector3(0,0,0),q);
                //showedge(new THREE.Vector3(0,0,0),e);
                //console.log(vdot(vsub(pivot,p1),vsub(q2,pivot)))
                //console.log("Adding vertices", vertices.length)
                lambdavertices = []
                thetavertices = []
                for (var i = 0; i < nvertices; i++) {
                    var n = addvertex(vmul(vertices[i],lambda),Color.red);
                    lambdavertices.push(n)
                    n = addvertex(vmul(vertices[i],theta), Color.green);
                    thetavertices.push(n)
                }
            }
            var p1 = lambdavertices[vlist[0]]
            var p2 = thetavertices[vlist[0]]
            var q1 = lambdavertices[vlist[1]]
            var q2 = thetavertices[vlist[1]]
            var pivot = addvertex(vmul(e,kappa),Color.blue)
            addedge(p1,pivot,Color.cyan);
            addedge(pivot,q2,Color.cyan);
            addedge(q1,pivot,Color.yellow);
            addedge(pivot,p2,Color.yellow);
        }
    }
    // I guess 'this' here is THREE.OFFLoader
    if (newfaces.length == 0 && !this.alerted) {
        alert("No edges found for edge linkage");
        this.alerted = true;
    }
    return { vertices: newvertices, faces: newfaces }
}

PolyContext.prototype.invert = function (off, options) {
    if (!off) return;
    // Invert in unit sphere
    var vdot = THREE.OFFLoader.Utils.vdot
    for (var i = 0; i < off.vertices.length; i++) {
        var k = vdot(off.vertices[i],off.vertices[i])
        off.vertices[i].divideScalar(k)
    }
    return off
}

//eqtest()
// function doit() {
//     vertices.push([0,0,0]);
//     for (var h = -3; h < 3; h += 0.05) {
//         bftest(h)
//     }
//     //bftest(0.0001)
//     console.log("OFF")
//     console.log(vertices.length + " " + edges.length + " " + 0)
//     for (var i = 0; i < vertices.length; i++) {
//         console.log(vertices[i][0] + " " + vertices[i][1] + " " + vertices[i][2])
//     }
//     for (var i = 0; i < edges.length; i++) {
//         console.log(2 + " " + edges[i][0] + " " + edges[i][1]);
//     }
// }

//doit()
