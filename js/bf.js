"use strict"

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
    options.bfheight = Number(options.bfheight);
    options.bfheightinc = Number(options.bfheightinc);
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
    var qmul = THREE.OFFLoader.Utils.qmul
    var vcross = THREE.OFFLoader.Utils.vcross
    var vnormalize = THREE.OFFLoader.Utils.vnormalize
    var Color = THREE.OFFLoader.Utils.Color;
    var vertices = []
    var faces = []
    // Apply a quaternion for a Clifford rotation
    if (!options.quat) {
        options.quat = new THREE.Vector4(1,0,0,0);
        options.quat.normalize();
        options.quinc = new THREE.Vector4(1,0.002,0.01,0.005);
        options.quinc.normalize();
    }
    var quat = options.quat;
    var quinc = options.quinc;
    qmul(quat,quinc,quat);
    function addvertex(v,color) {
        var i = vertices.length
        var qv = new THREE.Vector4(v.x,v.y,v.z,0)
        qmul(quat,qv,qv);
        vertices.push(qv);
        if (color) faces.push({ vlist: [i], color: color })
        return i
    }
    function addface(vlist,color) {
        faces.push({ vlist: vlist, color: color });
    }
    function addline(p,q,r) {
        addface([p,q],Color.white)
        if (r) {
            addface([q,r],Color.white)
            addface([r,p],Color.white)
        }
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
    var a = vector(0,-1,-1);
    var b = vector(-1,-1,-1);
    var c = vector(0,-1,0);
    var d = vector(0,1,-1);
    var e = vector(0,0,-1);
    var f = vector(1,-1,-1);
    var g = vector(0,-1,1);
    var h = intersect(e,f,b,d);
    var i = intersect(b,c,f,g);
    var j = intersect(g,e,c,d);
    var A = addvertex(a,Color.red);
    var B = addvertex(b,Color.blue);
    var C = addvertex(c,Color.blue);
    var D = addvertex(d,Color.blue);
    var E = addvertex(e,Color.green);
    var F = addvertex(f,Color.green);
    var G = addvertex(g,Color.green);
    var H = addvertex(h,Color.yellow);
    var I = addvertex(i,Color.yellow);
    var J = addvertex(j,Color.yellow);
    addline(A,D)
    addline(A,G)
    addline(B,D)

    addline(B,F)
    addline(B,I)
    addline(C,D)

    addline(E,G)
    addline(F,H)
    addline(F,G)

    addline(H,I)
    /*
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
*/
    return { vertices: vertices, faces: faces }
}

PolyContext.prototype.theorem = function(off,options) {
    var vector = THREE.OFFLoader.Utils.vector;
    var vadd = THREE.OFFLoader.Utils.vadd;
    var vsub = THREE.OFFLoader.Utils.vsub;
    var vdiv = THREE.OFFLoader.Utils.vdiv;
    var qmul = THREE.OFFLoader.Utils.qmul;
    var Color = THREE.OFFLoader.Utils.Color;
    var cos = Math.cos
    var sin = Math.sin
    var tan = Math.tan
   
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
    function bisect(v0,v1) {
        return vdiv(vadd(v0,v1),2);
    }
    // Construct a tetrahedron in 4-space
    var AQ = new THREE.Vector4(1,0,0,0);
    var BQ = new THREE.Vector4(0,1,0,0);
    var CQ = new THREE.Vector4(0,0,1,0);
    var DQ = new THREE.Vector4(0,0,0,1);
    // Apply a quaternion for a Clifford rotation
    if (!options.quat) {
        options.quat = new THREE.Vector4(1,0,-0.618,1);
        options.quat.normalize();
        options.quinc = new THREE.Vector4(1,0,0.01,0);
        options.quinc.normalize();
    }
    var quat = options.quat;
    var quinc = options.quinc;
    qmul(AQ,quat,AQ);
    qmul(BQ,quat,BQ);
    qmul(CQ,quat,CQ);
    qmul(DQ,quat,DQ);
    qmul(quat,quinc,quat);
    // And project into 3-space in the usual way.
    var A = vector(AQ.x,AQ.y,AQ.z);
    var B = vector(BQ.x,BQ.y,BQ.z);
    var C = vector(CQ.x,CQ.y,CQ.z);
    var D = vector(DQ.x,DQ.y,DQ.z);

    var E = bisect(A,C);
    var F = bisect(C,B);
    var G = bisect(A,D);
    var H = bisect(D,B);
    var I = bisect(G,F);
    // Now convert to indices (the joys of dynamic typing!)
    // Subtract I from everything to keep the centre at the centre.
    A = addvertex(vsub(A,I),Color.green);
    B = addvertex(vsub(B,I),Color.green);
    C = addvertex(vsub(C,I),Color.green);
    D = addvertex(vsub(D,I),Color.green);
    E = addvertex(vsub(E,I),Color.yellow);
    F = addvertex(vsub(F,I),Color.yellow);
    G = addvertex(vsub(G,I),Color.yellow);
    H = addvertex(vsub(H,I),Color.yellow);
    I = addvertex(vsub(I,I),Color.red);

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
    if (off.vertices) {
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
        var N = Number(options.n) || 3;
        var M = Number(options.m) || 1;
        var scale = 1;
        if (options.theta != undefined) {
            options.theta = Number(options.theta);
            console.log(options.theta)
            scale = 2*Math.sin(options.theta);
            options.theta += 0.01;
        }
        star = THREE.OFFLoader.makeStar(N,M,scale)
    }
    var k = Number(options.k) || 0
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
    var theta = Number(options.theta) || 0
    options.theta = theta + 0.01
    var n = Number(options.n) || 5
    var m = Number(options.m) || 1
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
    // There seems to be a problem in r74 where the indexed
    // elements don't get updated properly. Should move to
    // BufferGeometry, but in the meantime, force a clone
    // each time.
    return { vertices: vertices, faces: faces, needclone: true }
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
    var lambda = Number(options.lambda) || 1
    options.lambdainc = Number(options.lambdainc) || 0.01
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
    var vector = THREE.OFFLoader.Utils.vector
    var vdot = THREE.OFFLoader.Utils.vdot
    var vadd = THREE.OFFLoader.Utils.vadd
    var vsub = THREE.OFFLoader.Utils.vsub
    var vdiv = THREE.OFFLoader.Utils.vdiv
    if (!off) return;
    options.t = Number(options.t) || 0;
    var t = options.t;
    var icentre = vector(0.6*Math.sin(t),0.1*Math.cos(t/3),0);
    options.t += 0.01;
    // Invert in unit sphere
    var newvertices = [];
    for (var i = 0; i < off.vertices.length; i++) {
        var v = vsub(off.vertices[i],icentre);
        var k = vdot(v,v)*10;
        newvertices.push(vadd(icentre,vdiv(v,k)));
    }
    return { vertices: newvertices, faces: off.faces }
}

PolyContext.prototype.zoomer = (function () {
    var zoom = 2;
    var theta = 0;
    var types;
    return function(off,options) {
        var Color = THREE.OFFLoader.Utils.Color;
        var vmul = THREE.OFFLoader.Utils.vmul
        var vector = THREE.OFFLoader.Utils.vector
        var vertices = [];
        var faces = [];
        var z = 100;
        var depth = 15;
        var width = 0.1;
        var p = vector(50,50,-1000);
        vertices.push(p);
        var pindex = 0;
        var scale = zoom;
        if (!types) {
            types = [];
            for (var n = 0; n < 8*depth; n++) {
                types[n] = Math.floor(8*Math.random());
            }
        }
        var a = Math.cos(theta);
        var b = -Math.sin(theta);
        var c = Math.sin(theta);
        var d = Math.cos(theta);
        theta += 0.02;
        function rotate(v) {
            var x = v.x;
            var y = v.y;
            v.x = a*x + b*y;
            v.y = c*x + d*y;
            return v;
        }
        var index = 0;
        var color = Color.straw;
        for (var n = 0; n < depth; n ++) {
            for (var i = -1; i <= +1; i += 1) {
                for (var j = -1; j <= +1; j += 1) {
                    if (i == 0 && j == 0) continue;
                    var v0 = rotate(vmul(vector(i-width,j-width,0),scale));
                    var v1 = rotate(vmul(vector(i+width,j-width,0),scale));
                    var v2 = rotate(vmul(vector(i+width,j+width,0),scale));
                    var v3 = rotate(vmul(vector(i-width,j+width,0),scale));
                    var vindex = vertices.length;
                    vertices.push(v0);
                    vertices.push(v1);
                    vertices.push(v2);
                    vertices.push(v3);
                    faces.push({ vlist: [vindex+0,vindex+1,vindex+2,vindex+3], type: types[index++] });
                    faces.push({ vlist: [vindex+0,vindex+1,pindex], color: color });
                    faces.push({ vlist: [vindex+1,vindex+2,pindex], color: color });
                    faces.push({ vlist: [vindex+2,vindex+3,pindex], color: color });
                    faces.push({ vlist: [vindex+3,vindex+0,pindex], color: color });
                }
            }
            scale /= 2;
        }
        if (1) {
            zoom *= 1.05
            if (zoom > 4) {
                zoom /= 2;
                for (var n = 0; n < 8*depth; n++) {
                    if (n < 8*depth-8) {
                        types[n] = types[n+8];
                    } else {
                        types[n] = Math.floor(6*Math.random());
                    }
                }
            }
        } else {
            zoom *= 0.95
            if (zoom < 2) {
                zoom *= 2;
                for (n = 8*depth; n > 0; n--) {
                    if (n-1 >= 8) {
                        types[n-1] = types[n-9];
                    } else {
                        types[n-1] = Math.floor(6*Math.random());
                    }
                }
            }
        }
        return { vertices: vertices, faces: faces }
    }
})()

PolyContext.prototype.slerp = function(off,options) {
    var Vector3 = THREE.Vector3
    var Vector4 = THREE.Vector4
    var vnegate = THREE.OFFLoader.Utils.vnegate
    var vadd = THREE.OFFLoader.Utils.vadd
    var vadd3 = THREE.OFFLoader.Utils.vadd3
    var vsub = THREE.OFFLoader.Utils.vsub
    var vdiv = THREE.OFFLoader.Utils.vdiv
    var vmul = THREE.OFFLoader.Utils.vmul
    var vdot = THREE.OFFLoader.Utils.vdot
    var qmul = THREE.OFFLoader.Utils.qmul
    var vcross = THREE.OFFLoader.Utils.vcross
    var Color = THREE.OFFLoader.Utils.Color;

    var vertices = []
    var faces = []
    // We should rotate q1 with SLERP as well...
    var q0 = new Vector4(1,0.5,0.666,1);
    var qrot = new Vector4(1,0.01,0.011,0.008);
    q0.normalize(); qrot.normalize();
    var q1 = options.q1;
    if (!q1) {
        var q1 = new Vector4(0,1,0,0);
        q1.normalize();
        options.q1 = q1;
    }
    qmul(q1,qrot,q1);
    //console.log(q1);
    var k = vdot(q0,q1);
    var q2 = vsub(q1,vmul(q0,vdot(q0,q1)));
    q2.normalize();
    //console.log(q1);
    //console.log(q2);
    //console.log(dot(q0,q2));
    //console.log(dot(q1,q2));
    var zz = 1;
    var s = [ new Vector4(0,1,1,zz), new Vector4(0,1,-1,zz), new Vector4(0,-1,-1,zz), new Vector4(0,-1,1,zz) ];
    var colors = [Color.red,Color.green,Color.blue,Color.yellow];
    var N = 32;
    // A +iB is rotation by pi/N.
    var A = Math.cos(Math.PI/N), B = Math.sin(Math.PI/N);
    var a = 1, b = 0;

    var vbase = vertices.length;
    for (var j = 0; j < s.length; j++) {
        var t = s[j];
        console.assert(Math.abs(t.x) < 1e-5);
        vertices.push(new Vector3(t.y,t.z,t.w));
    }
    if (0) {
        for (var j = 0; j < s.length; j++) {
            faces.push({ vlist: [vbase+j], color: colors[j] });
            faces.push({ vlist: [vbase+j,vbase+(j+1)%s.length] });
        }
    }
    for(var i = 0; i < N; i++) {
        var quat = vadd(vmul(q0,a),vmul(q2,b));
        var cquat = new Vector4(quat.x,-quat.y,-quat.z,-quat.w);
        if (0) {
            var theta = i*Math.PI/N;
            console.assert(Math.abs(a-Math.cos(theta)) < 1e-5);
            console.assert(Math.abs(b-Math.sin(theta)) < 1e-5);
        }
        var a1 = a*A - b*B;
        var b1 = a*B + A*b;
        a = a1; b = b1;
        // Show rotation axis
        if (0) {
            var vbase = vertices.length;
            var vtmp = new Vector3(quat.y,quat.z,quat.w);
            //vtmp.normalize();
            //console.log(vtmp);
            vtmp.multiplyScalar(2);
            vertices.push(vtmp);
            vertices.push(vnegate(vtmp));
            faces.push({vlist:[vbase,vbase+1], color: Color.orange});
        }
        var vbase = vertices.length;
        for (var j = 0; j < s.length; j++) {
            var t = new Vector3();
            qmul(s[j],quat,t);
            qmul(cquat,t,t); // Omit this for a nice single-sided isoclinic rotation.
            console.assert(Math.abs(t.x) < 1e-5);
            vertices.push(new Vector3(t.y,t.z,t.w));
        }
        for (var j = 0; j < s.length; j++) {
            faces.push({ vlist: [vbase+j], color: colors[j] });
            faces.push({ vlist: [vbase+j,vbase+(j+1)%s.length] });
        }
    }
    return { vertices: vertices, faces: faces }
}

PolyContext.prototype.desargues2 = function(off,options) {
    var Vector3 = THREE.Vector3
    var Vector4 = THREE.Vector4
    var vnegate = THREE.OFFLoader.Utils.vnegate
    var vadd = THREE.OFFLoader.Utils.vadd
    var vadd3 = THREE.OFFLoader.Utils.vadd3
    var vsub = THREE.OFFLoader.Utils.vsub
    var vdiv = THREE.OFFLoader.Utils.vdiv
    var vmul = THREE.OFFLoader.Utils.vmul
    var vdot = THREE.OFFLoader.Utils.vdot
    var qmul = THREE.OFFLoader.Utils.qmul
    var vcross = THREE.OFFLoader.Utils.vcross
    var Color = THREE.OFFLoader.Utils.Color;
    var x = Math.sqrt(5);
    var pentatope =
        [Vector4(1,1,1,-1/x),
         Vector4(1,-1,-1,-1/x),
         Vector4(-1,1,-1,-1/x),
         Vector4(-1,-1,1,-1/x),
         Vector4(0,0,0,x-1/x)];
    var n = Vector4(0,0,0,1);
    var k = 1;
    var a = pentatope.map(function(p) { n.dot(p) });
    var p = [];
}

PolyContext.prototype.desargues2 = function(off,options) {
    var Vector3 = THREE.Vector3
    var Vector4 = THREE.Vector4
    var vnegate = THREE.OFFLoader.Utils.vnegate
    var vadd = THREE.OFFLoader.Utils.vadd
    var vadd3 = THREE.OFFLoader.Utils.vadd3
    var vsub = THREE.OFFLoader.Utils.vsub
    var vdiv = THREE.OFFLoader.Utils.vdiv
    var vmul = THREE.OFFLoader.Utils.vmul
    var vdot = THREE.OFFLoader.Utils.vdot
    var qmul = THREE.OFFLoader.Utils.qmul
    var vcross = THREE.OFFLoader.Utils.vcross
    var Color = THREE.OFFLoader.Utils.Color;
    var x = Math.sqrt(5);

    function random4() {
        var args = []
        for (var i = 0; i < 4; i++) {
            args.push(-1 + 2*Math.random());
        }
        return new Vector4(args[0],args[1],args[2],args[3]);
    }
    var vertices = [];
    var faces = [];
    var pentatope = options.pentatope;
    if (!pentatope) {
        if (0) {
            pentatope = [];
            for (var i = 0; i < 5; i++) {
                pentatope.push(random4());
            }
        } else {
            pentatope =
                [ new Vector4(1,1,1,-1/x),
                  new Vector4(1,-1,-1,-1/x),
                  new Vector4(-1,1,-1,-1/x),
                  new Vector4(-1,-1,1,-1/x),
                  new Vector4(0,0,0,x-1/x)
                ];
        }
        options.pentatope = pentatope;
    }
    var n = new Vector4(4,3,2,1);
    n.normalize();
    var lambda = 0;
    for (var i = 0; i < 5; i++) {
        vertices.push(pentatope[i]);
    }
    var pairs = [];
    for (var i = 0; i < 5; i++) {
        faces.push({ vlist: [i], color: Color.red });
        for (var j = i+1; j < 5; j++) {
            faces.push({ vlist: [i,j], color: Color.yellow });
            // Find a such that (a*v[i] + (1-a)v[j]).n = k
            // ie. a*(v[i].n - v[j].n) = k - v[j].n
            var ti = n.dot(pentatope[i])
            var tj = n.dot(pentatope[j])
            if (Math.abs(ti-tj) < 1e-5) {
                // Avoid division by zero or very small.
                var a = 10000; // Approximately infinity
            } else {
                var a = (lambda - tj)/(ti - tj);
            }
            var p = vadd(vmul(pentatope[i],a),vmul(pentatope[j],1-a));
            var k = vertices.length;
            pairs.push([i,j,k]);
            vertices.push(p);
            //console.log(p,a,ti,tj);
            faces.push({ vlist: [k], color: Color.green });
        }
    }
    function findpair(i,j) {
        for (var n = 0; n < pairs.length; n++) {
            if (pairs[n][0] == i && pairs[n][1] == j) {
                return pairs[n][2];
            }
        }
        console.assert(false);
    }
    for (var i = 0; i < 5; i++) {
        for (var j = i+1; j < 5; j++) {
            for (var k = j+1; k < 5; k++) {
                // I know the line segments overlap, I'm being lazy
                var p = findpair(i,j);
                var q = findpair(i,k);
                var r = findpair(j,k);
                var color = Color.white;
                faces.push({ vlist: [p,q], color: color });
                faces.push({ vlist: [q,r], color: color });
                faces.push({ vlist: [r,p], color: color });
            }
        }
    }
    for (var i = 0; i < 5; i++) {
        var q = new Vector4(1,0,0.005,0.01);
        q.normalize();
        qmul(q,pentatope[i],pentatope[i]);
        //var r = new Vector4(1,0.01,0,0);
        //r.normalize();
        //qmul(pentatope[i],r,pentatope[i]);
    }
    return { vertices: vertices, faces: faces }
}

PolyContext.prototype.tughra = function(off,options) {
    let A = 1, B = 1, C = 1;
    let D = 1, E = 1, F = 1;
    let G = 1, H = 1, I = 1;
    let J = 3;
    let N = 1000;
    if (options.A != undefined) A = options.A;
    if (options.B != undefined) B = options.B;
    if (options.C != undefined) C = options.C;
    if (options.D != undefined) D = options.D;
    if (options.E != undefined) E = options.E;
    if (options.F != undefined) F = options.F;
    if (options.G != undefined) G = options.G;
    if (options.H != undefined) H = options.H;
    if (options.I != undefined) I = options.I;
    if (options.J != undefined) J = options.J;
    if (options.N != undefined) N = options.N;

    var vector = THREE.OFFLoader.Utils.vector
    const sin = Math.sin;
    const cos = Math.cos;
    const vertices = [];
    const vlist = [];
    let x0 = 0, y0 = 0; z0 = 0;
    for (let i = 0; i < N; i++) {
        let theta0 = 2*Math.PI*i/N, theta = theta0;
        let x = 0, y = 0, z = 0;
        // For LR-symmetry, want x(-theta) = -x(theta), yz(-theta) = yz(theta)
        // So x should have odd powers of sin(n*theta)
        // y and z should have even powers of sin(n*theta)
        for (let j = 1; j <= J; j++) {
            // Could use eg. B*y*sin(theta) here, etc.
            let x1 = A*x + B*sin(theta) + C*sin(theta)*cos(theta);
            let y1 = D*y + E*cos(theta) + F*cos(theta)*cos(theta);
            let z1 = G*z + H*cos(theta) + I*sin(theta)*sin(theta);
            x = x1; y = y1; z = z1;
            theta += theta;
        }
        vertices.push(vector(x,y,z));
        x0 += x; y0 += y; z0 += z;
        vlist.push(i);
    }
    x0 /= N; y0 /= N; z0 /= N; 
    for (let i = 0; i < N; i++) {
        vertices[i].x -= x0;
        vertices[i].y -= y0;
        vertices[i].z -= z0;
    }
    return { vertices: vertices, faces: [{ vlist: vlist }] }
};

(function () {
        var vector = THREE.OFFLoader.Utils.vector
        var add = THREE.OFFLoader.Utils.vadd
        var sub = THREE.OFFLoader.Utils.vsub
        var mul = THREE.OFFLoader.Utils.vmul
        var dot = THREE.OFFLoader.Utils.vdot
        var cross = THREE.OFFLoader.Utils.vcross
        var dist = THREE.OFFLoader.Utils.vdist
        function length(v) { return Math.sqrt(dot(v,v)); }
        function solve1(a,b,j,k,r,parity) {
            // Solve v.a = j, v.b = k, |v| = r
            parity = parity || 1
            var c = cross(a,b)
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
            if (length(c) < 1e-4) {
                console.log("cross(a,b) vanishes " + a + b);
                return;
            }
            var aa = dot(a,a)
            var bb = dot(b,b)
            var cc = dot(c,c)
            var ab = dot(a,b)
            //console.log("Params1:", aa,bb,cc,ab)
            var B = (k*aa-j*ab)/(bb*aa - ab*ab)
            var A = (j - B*ab)/aa
            // v.v = aa + 2ab + bb + cc (since c orthogonal to a and b)
            var C2 = (r*r - A*A*aa - 2*A*B*ab - B*B*bb)/cc;
            //var v = add(mul(a,A),mul(b,B));
            //console.log("Params2:", A,B,C2);
            if (C2 < -1e-3) {
                console.log("# impossible length: ", C2);
                return;
            } else if (C2 < 0) {
                C2 = 0;
            }
            var C = parity*Math.sqrt(C2);
            var v = add(mul(a,A),add(mul(b,B),mul(c,C)));
            //console.log("Result:", v,dot(v,a),j,dot(v,b),k,length(v),r);
            return v;
        }

        function fold0(A,B,a,b,r,parity) {
            // Find P such that |P-A| = a, |P-B| = b, |P| = r
            // P.P = rr
            // (P-A).(P-A) = P.P - 2A.P + A.A = aa
            // A.P = 0.5*(rr + AA - aa)
            // B.P = 0.5*(rr + BB - bb)
            const P = solve1(A,B,
                             0.5*(r*r + dot(A,A) - a*a),
                             0.5*(r*r + dot(B,B) - b*b),
                             r, parity);
            //console.log(A,B,P);
            if (P) {
                const eps = 1e-3;
                console.assert(Math.abs(dist(P,A)-a) < eps);
                console.assert(Math.abs(dist(P,B)-b) < eps);
                console.assert(length(P)-r < eps);
            }
            return P;
        }
        function fold1(A,B,C,a,b,c,parity) {
            // return P such that |P-A| = a, etc.
            const P = fold0(sub(A,C),sub(B,C),a,b,c,parity);
            if (P) return add(C,P);
        }

    PolyContext.prototype.origami = function(off,options) {
        const vertices = []
        const faces = []
        function fold(A,B,C,a,b,c,parity) {
            A = vertices[A];
            B = vertices[B];
            C = vertices[C];
            return fold1(A,B,C,a,b,c,parity);
        }
        function face(p,q,r) {
            faces.push({ vlist:[p,q,r] });
        }
        function point(n,r,g,b) {
            faces.push({ vlist:[n], color: new THREE.Color(r,g,b) });
            //console.log(1,n,r,g,b);
        }
        function vertex(v) {
            //console.log(v[0],v[1],v[2]);
            console.assert(v);
            if (v) {
                const index = vertices.length;
                vertices.push(v);
                return index;
            } else {
                return 0;
            }
        }
        const t = options.t || 0;
        options.t = t + 0.015;
        const parity = 1;
        let a = 0.5*(1+Math.cos(t));

        const k = Math.sqrt(2);
        const j = 1/Math.cos(Math.PI/8);
        const m = Math.tan(Math.PI/8);
        //console.log("# ",j*j,m*m+1);
        if (Math.abs(a) < 1e-2) a = 1e-2;
        const x = -0.25;
        const A = vertex(vector(a,x,a)) // 0
        const B = vertex(vector(-a,x,a)) // 1
        const C = vertex(vector(a,x,-a)) // 2
        const O = vertex(fold(B,A,C,k,k,k,parity)); // 3
        const D = vertex(fold(C,O,A,j,1-m,j)); // 4
        const E = vertex(fold(A,O,B,j,1-m,j)); // 5
        const F = vertex(fold(A,D,C,1,m,1)); // 6
        const G = vertex(fold(B,E,A,1,m,1)); // 7
        //const H = fold(C,O,B,k,k,k,parity); // 8
        const H = vertex(vector(-a,x,-a)); //fold(C,O,B,k,k,k,parity); // 8
        const I = vertex(fold(B,O,H,j,1-m,j)); // 9
        const J = vertex(fold(H,O,C,j,1-m,j)); // 10
        const K = vertex(fold(H,I,B,1,m,1)); // 11
        const L = vertex(fold(C,J,H,1,m,1)); // 12
        //console.log("#",fold([t,0,0],[0,t,0],[0,0,0],k,k,k));
        //console.log("OFF");
        //console.log("8 16 0");
        //console.log("# ", dist(A,O), dist(B,O), dist(A,B));
        //console.log("# ", dist(O,D), dist(D,F));
        face(0,6,4);
        face(0,4,3);
        face(0,5,3);
        face(0,7,5);
        face(2,4,6);
        face(2,3,4);
        face(1,3,5);
        face(1,5,7);
        face(O,C,J);
        face(O,J,H);
        face(O,H,I);
        face(O,I,B);
        face(J,C,L);
        face(J,L,H);
        face(I,H,K);
        face(I,K,B);
        point(0,1,0,0);
        point(1,0,1,0);
        point(2,0,0,1);
        point(3,1,1,0);
        point(4,0,1,1);
        point(5,1,0,1);
        point(6,1,1,1);
        point(7,0.5,0.5,0.5);
        point(8,1,0.5,0);
        point(9,1,1,0.5);
        point(10,0.5,0,1);
        point(11,0,0.5,1);
        point(12,1,0,0.5);
        return { vertices: vertices, faces: faces }
    }
    PolyContext.prototype.miura = function(off,options) {
        const vertices = []
        const faces = []
        function fold(A,B,C,a,b,c,parity) {
            A = vertices[A];
            B = vertices[B];
            C = vertices[C];
            console.assert(A);
            console.assert(B);
            console.assert(C);
            const p = fold1(A,B,C,a,b,c,parity);
            if (!p) throw("Invalid");
            return p;
        }
        function face(p,q,r,s) {
            faces.push({ vlist:[p,q,r,s] });
        }
        function point(n,r,g,b) {
            faces.push({ vlist:[n], color: new THREE.Color(r,g,b) });
            //console.log(1,n,r,g,b);
        }
        function vertex(v) {
            //console.log(v[0],v[1],v[2]);
            console.assert(v);
            if (v) {
                const index = vertices.length;
                vertices.push(v);
                return index;
            } else {
                return 0;
            }
        }
        if (!options.tinc) {
            options.tinc = 0.01;
            options.t = 0;
        }
        //const t = Math.PI/4;
        const d = 1;
        const e = Math.sqrt(2)+0.3;
        const xbase = 0;
        let trying = true;
        while (trying) {
            try {
                let t = options.t;
                if (t < 0.01) {
                    t = 0.01;
                    options.tinc *= -1;
                }
                options.t = t + options.tinc;
                console.log("Running",options.t,options.tinc);
                let B0 = vertex(vector(0-3,xbase,0));
                let A0 = vertex(vector(-Math.cos(t)-3,xbase,-Math.sin(t)));
                let C0 = vertex(vector(-Math.cos(t)-3,xbase,Math.sin(t)));

                let M = 20;
                for (let i = 0; i < M; i++) {
                    const parity = i%2 ? 1 : -1;
                    const B1 = vertex(fold(A0,B0,C0,e,d,e,parity));
                    const A1 = vertex(add(vertices[B1],sub(vertices[A0],vertices[B0])));
                    const C1 = vertex(add(vertices[B1],sub(vertices[C0],vertices[B0])));
                    face(A0,B0,B1,A1);
                    face(C1,B1,B0,C0);
                    A0 = A1; B0 = B1; C0 = C1;
                }
                let row = [];
                const f = dist(vertices[0],vertices[4]); // B0,A1
                //console.log("f: ", f);
                for (let i = 0; i < M+1; i++) {
                    row[i] = 3*i+2; // All the C's
                    //point(row[i],1,0,1);
                }
                let N = 20;
                for (let j = 0; j < N; j++) {
                    let newrow = [];
                    // Now add a new row, connecting with the C's
                    for (let i = 1; i < M; i++) {
                        const parity = (i%2) ? -1 : 1;
                        newrow[i] = vertex(fold(row[i-1],row[i],row[i+1],(j%2)?f:e,1,(j%2)?e:f,parity));
                        //point(newrow[i],1,0,0);
                    }
                    newrow[0] = vertex(add(vertices[newrow[1]],sub(vertices[row[0]],vertices[row[1]])));
                    //point(newrow[0],1,1,0);
                    newrow[M] = vertex(add(vertices[newrow[M-1]],sub(vertices[row[M]],vertices[row[M-1]])));
                    //point(newrow[M],0,1,0);
                    for (let i = 0; i < M; i++) {
                        face(row[i],row[i+1],newrow[i+1],newrow[i]);
                    }
                    row = newrow;
                }
                trying = false;
            } catch (err) {
                if (err !== "Invalid") throw(err);
                console.log("Caught ", err);
                options.tinc *= -1;
                options.t += 2*options.tinc;
                faces = [];
                vertices = [];
            }
        }
                
        return { vertices: vertices, faces: faces }
    }
} ())
