"use strict";
(function () {
    // Construct a "desmic configuration" of 3 tetrahedra,
    // each edge of which intersects an edge of the other two,
    // and with each pair in perspective from each vertex of
    // the third.
    PolyContext.prototype.desmic = function(off,options) {
        var vector = THREE.OFFLoader.Utils.vector
        var qmul = THREE.OFFLoader.Utils.qmul
        var Color = THREE.OFFLoader.Utils.Color;
        var vertices = []
        var faces = []
        var time = this.stopwatch.getTime();
        // Show lines extending outwards in "projective" style.
        options.projective = true;
        
        // Apply a quaternion for a Clifford translation
        if (!options.quat) {
            options.quat = new THREE.Vector4(1,0,0,0);
            options.quat.normalize();
            options.quinc = new THREE.Vector4(1,0.002,0.002,0.002);
            options.quinc.normalize();
        }
        var quat = options.quat;
        var quinc = options.quinc;
        qmul(quat,quinc,quat);
        function addvertex(v,color) {
            var i = vertices.length
            var qv = new THREE.Vector4(v[0],v[1],v[2],v[3])
            qmul(quat,qv,qv);
            var w = qv.w;
            var eps = 1e-3;
            if (0 < w && w < eps) w = eps;
            if (0 > w && w > -eps) w = -eps;
            vertices.push(new THREE.Vector4(qv.x, qv.y, qv.z, qv.w));
            if (color) faces.push({ vlist: [i], color: color })
            return i
        }
        function addface(vlist,color) {
            faces.push({ vlist: vlist, color: color });
        }
        function addline(p,q,color) {
            color = color || Color.white;
            addface([p,q],color)
        }
        function minus(X) {
            return [-X[0],-X[1],-X[2],-X[3]];
        }
        function add(X,Y,Z,T) {
            return [X[0]+Y[0]+Z[0]+T[0],
                    X[1]+Y[1]+Z[1]+T[1],
                    X[2]+Y[2]+Z[2]+T[2],
                    X[3]+Y[3]+Z[3]+T[3]];
        }
        function tetra(x,y,z,t,color,facecolor) {
            facecolor = facecolor || Color.straw;
            addline(x,y,color); addline(y,z,color); addline(z,x,color);
            addline(x,t,color); addline(y,t,color); addline(z,t,color);
            // return; // Uncomment to just show lines.
            addface([x,y,z],facecolor);
            addface([x,y,t],facecolor);
            addface([x,z,t],facecolor);
            addface([y,z,t],facecolor);
        }
        var X = [1,0,0,0];
        var Y = [0,1,0,0];
        var Z = [0,0,1,0];
        var T = [0,0,0,1];
        var x = addvertex(X,Color.yellow);
        var y = addvertex(Y,Color.yellow);
        var z = addvertex(Z,Color.yellow);
        var t = addvertex(T,Color.yellow);
        tetra(x,y,z,t,Color.red);
        if (1) {
            addvertex([1,1,0,0],Color.orange);
            addvertex([1,0,1,0],Color.orange);
            addvertex([1,0,0,1],Color.orange);
            addvertex([0,1,1,0],Color.orange);
            addvertex([0,1,0,1],Color.orange);
            addvertex([0,0,1,1],Color.orange);

            addvertex([1,-1,0,0],Color.orange);
            addvertex([1,0,-1,0],Color.orange);
            addvertex([-1,0,0,1],Color.orange);
            addvertex([0,1,-1,0],Color.orange);
            addvertex([0,-1,0,1],Color.orange);
            addvertex([0,0,-1,1],Color.orange);
        }
        if (0) return { vertices: vertices, faces: faces }
        // a+b = [0,0,-2,2]
        // a-b = [2,-2,0,0]
        // a+d = [2,0,0,2]
        // a-d = [0,-2,-2,0]
        var a = addvertex([1,-1,-1,1],Color.yellow);
        var b = addvertex([-1,1,-1,1],Color.yellow);
        var c = addvertex([-1,-1,1,1],Color.yellow);
        var d = addvertex([1,1,1,1],Color.yellow);
        tetra(a,b,c,d,Color.green);
        var s = addvertex(add(X,Y,Z,minus(T)),Color.yellow);
        var p = addvertex(add(minus(X),Y,Z,T),Color.yellow);
        var q = addvertex(add(X,minus(Y),Z,T),Color.yellow);
        var r = addvertex(add(X,Y,minus(Z),T),Color.yellow);
        tetra(s,p,q,r,Color.blue);
        if (1) {
            addline(a,s,Color.white);
            addline(a,p,Color.white);
            addline(a,q,Color.white);
            addline(a,r,Color.white);

            addline(b,s,Color.white);
            addline(b,p,Color.white);
            addline(b,q,Color.white);
            addline(b,r,Color.white);

            addline(c,s,Color.white);
            addline(c,p,Color.white);
            addline(c,q,Color.white);
            addline(c,r,Color.white);

            addline(d,s,Color.white);
            addline(d,p,Color.white);
            addline(d,q,Color.white);
            addline(d,r,Color.white);
        }
        return { vertices: vertices, faces: faces }
    }
} ())
