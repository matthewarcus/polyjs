"use strict";

// The MIT License (MIT)

// New material: Copyright (c) 2016 Matthew Arcus

// Based on THREE.js OBJLoader.js
// Original attribution:
// @author mrdoob / http://mrdoob.com/

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

THREE.OFFLoader = function ( manager ) {
    this.manager = ( manager !== undefined ) ? manager : THREE.DefaultLoadingManager;
};

THREE.OFFLoader.prototype = {
    constructor: THREE.OFFLoader,
    load: function ( url, geometry, options, onLoad, onProgress, onError ) {
	var scope = this;
	var loader = new THREE.XHRLoader( scope.manager );
	loader.setCrossOrigin( this.crossOrigin );
	loader.load( url, function ( text, other ) {
            var off = scope.parse(text)
            if (options.theta != undefined) {
                off = THREE.OFFLoader.dipolygonid(off)
                THREE.OFFLoader.twistertransform(off,options.theta)
            } else if (options.lambda != undefined) {
                off = THREE.OFFLoader.edgelinkage(off,options.lambda)
            }
	    onLoad(THREE.OFFLoader.display(off, geometry, options));
	}, onProgress, onError );
    },

    setCrossOrigin: function ( value ) {
	this.crossOrigin = value;
    },

    verbose: false,
}

THREE.OFFLoader.defaultoptions = function () {
    return {
        vertexstyle: "sphere",
        vertexwidth: 0.06,
        edgewidth: 0.02,
        alledges: false,
        allvertices: false,
        novertices: false,
        noedges: false,
        nofaces: false,
        dospline: false,
        reversefaces: false,
    }
};

// Need to put this somewhere better
THREE.OFFLoader.Utils = {
    vector: function (x,y,z) {
        return new THREE.Vector3(x,y,z)
    },
    vnegate: function (v) {
        var t = v.clone()
        return t.negate()
    },
    vadd: function (v0, v1) {
        var t = v0.clone();
        t.add(v1);
        return t;
    },
    vadd3: function (v0, v1, v2) {
        var t = v0.clone();
        t.add(v1);
        t.add(v2);
        return t;
    },
    vsub: function (v0, v1) {
        var t = v0.clone();
        t.sub(v1);
        return t;
    },
    vmul: function (v0,theta) {
        var v = v0.clone();
        v.multiplyScalar(theta);
        return v;
    },
    vdiv: function (v0,theta) {
        var v = v0.clone();
        v.divideScalar(theta);
        return v;
    },
    vcross: function (v0, v1) {
        var t = v0.clone();
        t.cross(v1);
        return t;
    },
    vdot: function (v0, v1) {
        return v0.dot(v1);
    },
    vdist: function (v0, v1) {
        return v0.distanceTo(v1);
    },
    interp: function (v0,v1,theta) {
        var v = vmul(v0,theta);
        v.add(vmul(v1,1-theta));
        return v;
    },
    makenormal: function (v0,v1,v2) {
        var vsub = THREE.OFFLoader.Utils.vsub;
        var n = vsub(v2,v0);
        n.cross(vsub(v1,v0));
        n.normalize();
        return n;
    },
    perpdist: function (v,v0,v1) {
        var n = v1.clone(); // n = v1
        n.sub(v0); // n = v1 - v0
        n.divideScalar(n.length()); // n = n0 = (v1 - v0)/|v1 - v0|
        // n now along line from v0 to v1
        var t = v.clone();
        t.sub(v0); // t = v - v0
        n.multiplyScalar(n.dot(t)); // n = n0(n0.(v-v0))
        n.add(v0); // n = v0 + n0(n0.v)
        // n is point on v0,v1 nearest v
        //console.log(vsub(n,v0).dot(vsub(v,n)))
        return n.distanceTo(v);
    },
    Color: {
        white: new THREE.Color(1,1,1),
        red: new THREE.Color(1,0,0),
        yellow:new THREE.Color(1,1,0),
        green: new THREE.Color(0,1,0),
        cyan: new THREE.Color(0,1,1),
        blue: new THREE.Color(0,0,1),
        magenta:new THREE.Color(1,0,1),
        orange: new THREE.Color(1,0.5,0),
        raspberry: new THREE.Color(1,0,0.5),
        purple: new THREE.Color(1,0.5,1),
    }
}


THREE.OFFLoader.prototype.parse = function ( text, geometry, options ) {
    // Remove all occurences of a from s
    function remove(a,s) {
        while(true) {
            var i = s.indexOf(a);
            if (i < 0) break;
            s.splice(i,1);
        }
        return s;
    }
    if (this.verbose) console.time( 'OFFLoader' );
    var vertices = [];
    var faces = [];
    var state = 0, count = 0, nvertices = 0, nfaces = 0;

    var lines = text.split( '\n' );
    for (var lineindex = 0; lineindex < lines.length; lineindex++) {
	var line = lines[ lineindex ];
	line = line.trim();
	if ( line.length === 0 || line.charAt( 0 ) === '#' ) continue;
        var items = remove("",line.split(" "));
        if (state == 0) {
            if (items.length != 1 || items[0] != "OFF") {
                alert("File load error: invalid OFF line: " + line);
                return;
            }
            state = 1;
        } else if (state == 1) {
            if (items.length == 3) {
                nvertices = parseInt(items[0]);
                nfaces = parseInt(items[1]);
                state = 2;
            }
        } else if (state == 2) {
            var x = parseFloat(items[0]);
            var y = parseFloat(items[1]);
            var z = parseFloat(items[2]);
            if (items[3]) {
                // Crude attempt to deal with homogenous coordinates.
                // w actually zero is problematic, so just make it very small.
                // eps can't be too big or texture interpolation goes awry
                var eps = 1e-4;
                var w = parseFloat(items[3]);
                if (w < 0) {
                    w = -w; x = -x; y = -y; z = -z;
                }
                if (w < eps) w = eps;
                x /= w; y /= w; z /= w;
            }
            vertices.push(new THREE.Vector3(x,y,z));
            nvertices--;
            if (nvertices == 0) {
                state = 3;
            }
        } else if (state == 3) {
            if (nfaces > 0) {
                nfaces--;
                var color;
                var n = parseInt(items[0]);
                // Color can also just be an index
                if (items.length == n+4) {
                    // Reuse color objects?
                    color = new THREE.Color(parseFloat(items[n+1]),
                                            parseFloat(items[n+2]),
                                            parseFloat(items[n+3]));
                }
                var vlist = [];
                for (var i = 0; i < n; i++) {
                    var v = parseInt(items[i+1]);
                    vlist.push(v);
                }
                faces.push({ vlist: vlist, color: color });
            }
        }
    }
    if (this.verbose) console.timeEnd( 'OFFLoader' );
    return { vertices: vertices, faces: faces }
}

// Copy of Geometry.merge, but updates the geometry in place as far as possible
THREE.Geometry.prototype.mlaMerge = function ( geometry, matrix, materialIndexOffset ) {
    if ( geometry instanceof THREE.Geometry === false ) {
	console.error( 'THREE.Geometry.mlaMerge(): geometry not an instance of THREE.Geometry.', geometry );
	return;
    }
    var normalMatrix,
	//vertexOffset = this.vertices.length,
	vertexOffset = this.mlaVertexOffset,
	vertices1 = this.vertices,
	vertices2 = geometry.vertices,
	//faceOffset = this.faces.length,
	faceOffset = this.mlaFaceOffset,
	faces1 = this.faces,
	faces2 = geometry.faces,
	//uvOffset = this.faceVertexUvs[ 0 ].length,
	uvOffset = this.mlaUvOffset,
	uvs1 = this.faceVertexUvs[ 0 ],
	uvs2 = geometry.faceVertexUvs[ 0 ];
    //console.log(this.mlaVertexOffset, vertices1.length, this.mlaFaceOffset, faces1.length)
    this.mlaVertexOffset += vertices2.length
    this.mlaFaceOffset += faces2.length
    this.mlaUvOffset += uvs2.length
    if ( materialIndexOffset === undefined ) materialIndexOffset = 0;
    if ( matrix !== undefined ) {
	normalMatrix = new THREE.Matrix3().getNormalMatrix( matrix );
    }
    // vertices
    while (vertices1.length < vertexOffset+vertices2.length ) {
        vertices1.push(new THREE.Vector3);
    }
    for ( var i = 0, il = vertices2.length; i < il; i ++ ) {
	var vertex2 = vertices2[ i ];
        var vertex1 = vertices1[vertexOffset + i]
        vertex1.x = vertex2.x
        vertex1.y = vertex2.y
        vertex1.z = vertex2.z
	if ( matrix !== undefined ) vertex1.applyMatrix4( matrix );
    }
    // faces
    while (faces1.length < faceOffset + faces2.length) {
        faces1.push(new THREE.Face3);
    }
    for ( i = 0, il = faces2.length; i < il; i ++ ) {
	var face2 = faces2[ i ], face1, color,
	    face2VertexNormals = face2.vertexNormals,
	    face2VertexColors = face2.vertexColors;
        face1 = faces1[faceOffset+i]
	var face1VertexNormals = face1.vertexNormals
	var face1VertexColors = face1.vertexColors

	//faceCopy = new THREE.Face3( face.a + vertexOffset, face.b + vertexOffset, face.c + vertexOffset );
        face1.a = face2.a + vertexOffset;
        face1.b = face2.b + vertexOffset;
        face1.c = face2.c + vertexOffset;
	face1.normal.copy( face2.normal );
	if ( normalMatrix !== undefined ) {
	    face1.normal.applyMatrix3( normalMatrix ).normalize();
	}
        while (face1VertexNormals.length < face2VertexNormals.length) {
            face1VertexNormals.push(new THREE.Vector3);
        }
        console.assert(face2VertexNormals.length == 0 || face2VertexNormals.length == 3)
	for ( var j = 0, jl = face2VertexNormals.length; j < jl; j ++ ) {
	    var normal1 = face1VertexNormals[ j ];
	    var normal2 = face2VertexNormals[ j ];
            normal1.copy(normal2)
	    if ( normalMatrix !== undefined ) {
		normal1.applyMatrix3( normalMatrix ).normalize();
	    }
	    //face1.vertexNormals.push( normal );
	}
	face1.color.copy( face2.color );
        while (face1VertexColors.length < face2VertexColors.length) {
            face1VertexColors.push(new THREE.Color)
        }
	for ( var j = 0, jl = face2VertexColors.length; j < jl; j ++ ) {
	    var color1 = face1VertexColors[ j ];
	    var color2 = face2VertexColors[ j ];
            color1.copy(color2)
	    //face1.vertexColors.push( color.clone() );
	}
	face1.materialIndex = face2.materialIndex + materialIndexOffset;
    }
    // uvs
    while (uvs1.length < uvOffset+uvs2.length) {
        uvs1.push( [new THREE.Vector2,new THREE.Vector2,new THREE.Vector2] );
    }
    for ( i = 0, il = uvs2.length; i < il; i++ ) {
        var uv1 = uvs1[uvOffset + i];
        //console.assert(uv1)
        //console.assert(uv1.length == 3)
	var uv2 = uvs2[ i ]
        //console.assert(uv2)
        //console.assert(uv2.length == 3) // Check is triangle
	for ( var j = 0, jl = uv2.length; j < jl; j ++ ) {
            uv1[j].x = uv2[j].x
            uv1[j].y = uv2[j].y
	    //uvCopy.push( new THREE.Vector2( uv[ j ].x, uv[ j ].y ) );
	}
	//uvs1.push( uvCopy );
    }
}

THREE.OFFLoader.display = function ( off, geometry, options ) {
    var vadd = THREE.OFFLoader.Utils.vadd
    var vsub = THREE.OFFLoader.Utils.vsub
    var vdiv = THREE.OFFLoader.Utils.vdiv
    var vmul = THREE.OFFLoader.Utils.vmul
    var vdot = THREE.OFFLoader.Utils.vdot
    var vcross = THREE.OFFLoader.Utils.vcross
    var makenormal = THREE.OFFLoader.Utils.makenormal

    // Options
    var vertexstyle = options.vertexstyle;
    var vertexwidth = options.vertexwidth;
    var edgewidth = options.edgewidth;
    var novertices = options.novertices;
    var noedges = options.noedges;
    var nofaces = options.nofaces;
    var alledges = options.alledges;
    var allvertices = options.allvertices;
    var dospline = options.dospline;
    var reversefaces = options.reversefaces;
    var docut = false;
    
    //var geometry = new THREE.Geometry();
    geometry.mlaVertexOffset = 0
    geometry.mlaFaceOffset = 0
    geometry.mlaUvOffset = 0
    var vertices = geometry.vertices;
    var faces = geometry.faces;
    var faceVertexUvs = geometry.faceVertexUvs[0];
    //faceVertexUvs[0] = [];
    function addface(a, b, c, uva, uvb, uvc, normal, color) {
        console.assert(geometry.mlaFaceOffset == geometry.mlaUvOffset)
        if (faces.length <= geometry.mlaFaceOffset) {
            faces.push(new THREE.Face3)
            faceVertexUvs.push([])
        }
        var i = geometry.mlaFaceOffset;
        geometry.mlaFaceOffset++;
        geometry.mlaUvOffset++;
	faces[i].a = a
	faces[i].b = b
	faces[i].c = c
	faces[i].normal = normal
	if (color) faces[i].color = color
        faceVertexUvs[i] = [uva,uvb,uvc];
    }
    function addvertex(vertex) {
        if (vertices.length <= geometry.mlaVertexOffset) {
            vertices.push(new THREE.Vector3);
        }
        var i = geometry.mlaVertexOffset;
        vertices[i].copy(vertex);
        geometry.mlaVertexOffset++;
        return i;
    }
    var vertexmodel;
    var edgemodel;
    var defaultvertexcolor = new THREE.Color(0.9,0.9,0.9);
    var defaultedgecolor = defaultvertexcolor;
    var m = new THREE.Matrix4;
    var tt = new THREE.Matrix4;
    var yaxis = new THREE.Vector3(0,1,0);
    var xaxis = new THREE.Vector3(1,0,0);
    function showvertex(v0, color) {
        //console.log(v0,color)
        if (!vertexmodel) {
            if (vertexstyle === "cylinder") {
                vertexmodel = new THREE.CylinderGeometry( vertexwidth, vertexwidth, vertexwidth );
            } else if (vertexstyle === "icosahedron") {
                vertexmodel = new THREE.IcosahedronGeometry( vertexwidth, 2 );
            } else {
                if (vertexstyle !== "sphere") alert("Unknown vertexstyle: " + vertexstyle);
                vertexmodel = new THREE.SphereGeometry( vertexwidth );
            }
        }
        m.identity();
        var vlen = v0.length();
        if (vlen > 1e-4) {
            var raxis = vcross(yaxis,v0);
            //console.log("raxis",v0,raxis,raxis.length())
            var theta = Math.acos(vdot(yaxis,v0)/vlen);
            if (raxis.length() < 1e-4) {
                // v0 near parallel to yaxis
                // so use xaxis as rotation
                // theta we hope is 0 or pi
                raxis = xaxis;
            } else {
                raxis.normalize();
            }
            tt.makeRotationAxis(raxis,theta);
            m.multiply(tt);
            tt.makeTranslation(0,v0.length(),0);
            m.multiply(tt);
        }
        for (var i = 0; i < vertexmodel.faces.length; i++) {
            vertexmodel.faces[i].color = color || defaultvertexcolor;
        }
        geometry.mlaMerge(vertexmodel,m);
    }
    var edge = new THREE.Vector3;  // Probably an absurd microptimization
    var midpoint = new THREE.Vector3; // Ditto.
    function showedge(v0, v1, color) {
        // Make cylinder, stretch to right length, rotate
        // to the correct orientation, then translate into position.
        edge.subVectors(v1,v0);
        var edgelen = edge.length();   // Total height of cylinder
        midpoint.addVectors(v0,v1); // Map cylinder origin to here.
        midpoint.divideScalar(2);

        m.identity();
        tt.makeTranslation(midpoint.x,midpoint.y,midpoint.z);
        m.multiply(tt);

        var raxis = vcross(yaxis,edge);
        if (raxis.length() > 1e-4) {
            var theta = Math.acos(vdot(yaxis,edge)/edgelen);
            raxis.normalize();
            tt.makeRotationAxis(raxis,theta);
            m.multiply(tt);
        }
        tt.makeScale(1,edgelen,1);
        m.multiply(tt);
        if (!edgemodel) {
            edgemodel = new THREE.CylinderGeometry( edgewidth, edgewidth, 1);
        }
        for (var i = 0; i < edgemodel.faces.length; i++) {
            edgemodel.faces[i].color = color || defaultedgecolor;
        }
        geometry.mlaMerge(edgemodel,m);
    }
    var cut = [];
    var cutlim = -2;
    var cutnormal = new THREE.Vector3(1,2,3);
    cutnormal.normalize();

    var nvertices = off.vertices.length
    var nfaces = off.faces.length
    for (var i = 0; i < nvertices; i++) {
        addvertex(off.vertices[i]);
    }
    // Mustn't call showvertex until after all vertices
    // have been set up or the indexing goes awry.
    if (allvertices) {
        for (var i = 0; i < nvertices; i++) {
            showvertex(vertices[i]); // This extends vertices
        }
    }
    for (var faceindex = 0; faceindex < nfaces; faceindex++) {
        var n = off.faces[faceindex].vlist.length
        var vlist = off.faces[faceindex].vlist
        var color = off.faces[faceindex].color
        for (var i = 0; i < vlist.length; i++) {
            console.assert(vlist[i] < vertices.length,
                           "Vertex index out of range");
        }
        if (n == 1) {
            if (!novertices) {
                var v0 = vertices[vlist[0]];
                showvertex(v0,color);
            }
        } else if (n == 2) {
            if (!noedges) {
                var p = vertices[vlist[0]];
                var q = vertices[vlist[1]];
                showedge(p,q,color);
            }
        } else {
            if (docut) {
                var t = [];
                var cutcount = 0;
                for (var i = 0; i < n; i++) {
                    var v0 = vertices[vlist[i]];
                    var v1 = vertices[vlist[(i+1)%n]];
                    var test0 = v0.dot(cutnormal);
                    var test1 = v1.dot(cutnormal);
                    // Should do something about division by 0 here
                    if (test0 >= cutlim) t.push(vlist[i]);
                    if ((test0 < cutlim && test1 >= cutlim) ||
                        (test0 >= cutlim && test1 < cutlim)) {
                        var k = (test1 - cutlim)/(test1 - test0);
                        //console.log(v0);
                        //console.log(v1);
                        //console.log(k);
                        newvertex = addvertex(interp(v0,v1,k));
                        t.push(newvertex);
                        cut.push(newvertex);
                        cutcount++;
                        //console.log("Added 2 " + t[t.length-1] + " " + v0.y + " " + v1.y + " " + k);
                    }
                }
                vlist = t;
                n = vlist.length;
                if (n <= 2) continue;
            }
            var centroid = new THREE.Vector3();
            for (var i = 0; i < n; i++) {
                centroid.add(vertices[vlist[i]]);
            }
            centroid.divideScalar(n);
            if (alledges) {
                for (var i = 0; i < n; i++) {
                    var v0 = vertices[vlist[i]];
                    var v1 = vertices[vlist[(i+1)%n]];
                    showedge(v0,v1);
                }
            }
            if (dospline) {
                var s = [];
                for (var i = 0; i < n; i++) {
                    s.push(vertices[vlist[i]]);
                }
                var curve = new THREE.ClosedSplineCurve3(s);
                var geom = new THREE.TubeGeometry(curve, n, edgewidth, 8, true);
                if (color) {
                    for (var i = 0; i < geom.faces.length; i++) {
                        geom.faces[i].color = color;
                    }
                }
                geometry.mlaMerge(geom);
            }
            if (!nofaces) {
                var icentroid = addvertex(centroid);
                var normal = new THREE.Vector3;
                // Average the normals over the whole face.
                // We could chance it and just check one pair of vertices.
                for (var i = 0; i < n; i++) {
                    var v1 = vertices[vlist[i]];
                    var v2 = vertices[vlist[(i+1)%n]];
                    normal.add(makenormal(centroid,v2,v1)); // Note order
                }
                normal.normalize();
                if (reversefaces) normal.negate();
                // Now we have a normal, let's work out texture coordinates
                var uaxis = vertices[vlist[0]].clone();
                uaxis.sub(centroid);
                uaxis.normalize();
                var vaxis = vcross(uaxis,normal);
                var uvs = [];
                for (var i = 0; i < n; i++) {
                    var tmp = vsub(vertices[vlist[i]], centroid);
                    uvs.push(new THREE.Vector2(0.5+vdot(tmp,uaxis), 0.5+vdot(tmp,vaxis)));
                }
                var uvcentroid = new THREE.Vector2(0.5,0.5);
                for (var i = 0; i < n; i++) {
                    var p0 = vlist[i];
                    var p1 = vlist[(i+1)%n];
                    var uv0 = uvs[i];
                    var uv1 = uvs[(i+1)%n];
                    //normal = makenormal(centroid,vertices[p1],vertices[p0]) // Per edge normal
                    if (reversefaces) addface(icentroid, p1, p0, uvcentroid, uv1, uv0, normal, color);
                    else addface(icentroid, p0, p1, uvcentroid, uv0, uv1, normal, color);
                }
            }
        }
    }
    if (docut) {
        var cutcentre = new THREE.Vector3;
        for (var i = 0; i < cut.length; i++) {
            //console.log(JSON.stringify(vlist[cut[i]]));
            cutcentre.add(vertices[cut[i]]);
        }
        cutcentre.divideScalar(cut.length);
        icutcentre = addvertex(cutcentre);
        color = new THREE.Color(0.9,0.5,0.5);
        console.log(cut.length + " cut vertices");
        for (var i = 0; i < cut.length; i+=2) {
            //addface(icutcentre,cut[i],cut[i+1],color);
        }
    }
    if (0) {
        faceVertexUvs = []
        for (var i = 0; i < faces.length ; i++) {
            faceVertexUvs.push(
                [new THREE.Vector2(0,0),
                 new THREE.Vector2(0,1),
                 new THREE.Vector2(1,1)])
        }
    }
    // Trim down in case the geometry has shrunk
    if (vertices.length > geometry.mlaVertexOffset) vertices.length = geometry.mlaVertexOffset
    console.assert(faces.length == faceVertexUvs.length)
    console.assert(geometry.mlaFaceOffset == geometry.mlaUvOffset)
    if (faces.length > geometry.mlaFaceOffset) {
        faces.length = geometry.mlaFaceOffset
        faceVertexUvs.length = geometry.mlaUvOffset
    }
    if (this.verbose) console.log("Vertices = " + vertices.length + " faces = " + faces.length);
    return geometry;
}

THREE.OFFLoader.edgelinkage = function (off, lambda) {
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
    function addedge(p1,p2,color) {
        //console.log("Face",p1,p2,color)
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
    faces.forEach(function(face) {
        var vlist = face.vlist
        if (vlist.length == 2) {
            var color = face.color
            var p = vertices[vlist[0]];
            var q = vertices[vlist[1]];
            var e = vdiv(vadd(p,q),2)    //edge vector
            if (lambdavertices === undefined) {
                // Should work out a proper approximation for this
                if (lambda > 0.9999) lambda = 0.9999
                if (lambda < -0.9999) lambda = -0.9999
                var len = vdist(p,q)
                var A = vdot(q,q)
                var B = -2*lambda*vdot(p,q)
                var C = lambda*lambda*vdot(p,p) - len*len
                // TBD: check on stable quadratic solutions
                var disc = B*B - 4*A*C
                console.assert(disc >= 0, "No solution to edge link")
                // Either + or - is realizable of course.
                var theta = (-B + Math.sqrt(disc))/(2*A)
                var p1 = vmul(p,lambda)
                var p2 = vmul(p,theta)
                var q1 = vmul(q,lambda)
                var q2 = vmul(q,theta)
                var en = vdiv(vadd(p1,q2),2) // new edge centre
                var n = vdiv(vsub(p1,q2),len) // normal to bisecting plane
                kappa = vdot(en,n)/vdot(e,n) // distance to plane
                //console.log("Plane ", n, n.length(), kappa)
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
    })
    // I guess 'this' here is THREE.OFFLoader
    if (newfaces.length == 0 && !this.alerted) {
        alert("No edges found for edge linkage");
        this.alerted = true;
    }
    return { vertices: newvertices, faces: newfaces }
}

THREE.OFFLoader.dipolygonid0 = function (off) {
    var vadd = THREE.OFFLoader.Utils.vadd
    var vsub = THREE.OFFLoader.Utils.vsub
    var vdiv = THREE.OFFLoader.Utils.vdiv
    var vmul = THREE.OFFLoader.Utils.vmul
    var vdot = THREE.OFFLoader.Utils.vdot
    var vcross = THREE.OFFLoader.Utils.vcross
    var Color = THREE.OFFLoader.Utils.Color
    var vertices = off.vertices
    var faces = off.faces
    var newvertices = []
    var newfaces = []
    var edges = new Map()
    faces.forEach(function(face,findex) {
        var vlist = face.vlist
        var newvlist = []
        vlist.forEach(function(v1,i) {
            var v2 = vlist[(i+1)%vlist.length]
            console.assert(v2 != undefined)
            //console.log("Setting",v1,v2)
            edges.set(v1*1000+v2, { face: findex, index: i })
            newvlist.push(newvertices.length)
            newvertices.push(vertices[v1].clone())
        })
        newfaces.push({ vlist: newvlist, color: Color.red, type: 0 })
    })
    faces.forEach(function(face,findex) {
        var vlist = face.vlist
        var newvlist = []
        vlist.forEach(function(v1,i) {
            var v2 = vlist[(i+1)%vlist.length]
            console.assert(v2 != undefined)
            var e = edges.get(v2*1000+v1) // Other face for edge
            //console.log("Getting",v2,v1,e)
            console.assert(e)
            var f2 = e.face
            var i2 = e.index
            var v3 = newfaces[f2].vlist[i2]
            console.assert(v3 != undefined)
            newvlist.push(v3)
        })
        newfaces.push({ vlist: newvlist, color: Color.green, type: 1 })
    })
    //console.log(newvertices.length, newfaces.length)
    return { vertices: newvertices, faces: newfaces }
}

// Take an ordinary OFF model, with up to 3 set of proper faces,
// and each vertex being part of one of each type of face.
THREE.OFFLoader.dipolygonid = function (off) {
    var vadd = THREE.OFFLoader.Utils.vadd
    var vsub = THREE.OFFLoader.Utils.vsub
    var vdiv = THREE.OFFLoader.Utils.vdiv
    var vmul = THREE.OFFLoader.Utils.vmul
    var vdot = THREE.OFFLoader.Utils.vdot
    var vcross = THREE.OFFLoader.Utils.vcross
    var Color = THREE.OFFLoader.Utils.Color
    var vertices = off.vertices
    var faces = off.faces
    var newvertices = []
    var newfaces = []
    var emap = new Map()
    var vmap = new Map()
    // Represent the type of a face by the dot product of
    // the first two edges. This will distinguish eg. a pentagon
    // from a pentagram, unlike plain edge count.
    function facetype(face) {
        var vlist = face.vlist        
        var type = vdot(vsub(vertices[vlist[1]], vertices[vlist[0]]),
                        vsub(vertices[vlist[2]], vertices[vlist[1]]));
        return type;
    }
    function eqtype(type1,type2) {
        return Math.abs(type1-type2) < 1e-4;
    }
    var primarytype;
    var secondarytype;
    faces.forEach(function(face,findex) {
        var vlist = face.vlist
        if (vlist.length < 3) return;
        var type = facetype(face)
        if (!primarytype) primarytype = type
        if (eqtype(type,primarytype)) {
            // Generate new edge points for each primary face
            var newvlist = vlist.map(function(v) {
                var i = newvertices.length
                newvertices.push(vertices[v].clone())
                // What about vertices shared between several primary faces?
                vmap.set(v,i)
                return i
            })
            // Create a new face
            newfaces.push({ vlist: newvlist, color: Color.red, type: 0 })
            // Remember vertex mapping for each edge
            for (var i = 0; i < vlist.length; i++) {
                var v1 = vlist[i], v2 = vlist[(i+1)%vlist.length]
                var nv1 = newvlist[i], nv2 = newvlist[(i+1)%newvlist.length]
                emap.set(v1*1000+v2, { face: findex, index: i, v1: nv1, v2: nv2 })
            }
        } else {
            if (!secondarytype) secondarytype = type
            var valid = true;
            var newvlist = []
            for (var i = 0; i < vlist.length; i++) {
                var v1 = vlist[i], v2 = vlist[(i+1)%vlist.length]
                var nv1;
                var e = emap.get(v2*1000+v1)
                if (e) nv1 = e.v2 // Get the second vertex in the stored list
                else nv1 = vmap.get(v1) // Otherwise, just use the mapped value
                if (nv1 == undefined) {
                    valid = false;
                    break;
                }
                newvlist.push(nv1)
            }
            if (valid) {
                if (eqtype(type,secondarytype)) {
                    newfaces.push({ vlist: newvlist, color: Color.green, type: 1 })
                } else {
                    newfaces.push({ vlist: newvlist, color: Color.yellow, type: 2 })
                }
            }
        }
    })
    if (secondarytype == undefined) {
        // Didn't find a second face type, duplicate the first set of faces
        faces.forEach(function(face,findex) {
            var vlist = face.vlist
            if (vlist.length >= 3) {
                var newvlist = []
                vlist.forEach(function(v1,i) {
                    var v2 = vlist[(i+1)%vlist.length]
                    console.assert(v2 != undefined)
                    var e = emap.get(v2*1000+v1) // Other face for edge
                    console.assert(e)
                    var v3 = e.v2
                    console.assert(v3 != undefined)
                    newvlist.push(v3)
                })
                newfaces.push({ vlist: newvlist, color: Color.green, type: 1 })
            }
        })
    }
    //console.log(newvertices.length, newfaces.length)
    return { vertices: newvertices, faces: newfaces }
}

// Rotate the "primary" faces of the off model, then translate
// each face enough to make the "secondary" faces have the same
// edge width. Sometimes there isn't a solution. For three faces
// occasionally it works, but usually a face is distorted.
THREE.OFFLoader.twistertransform = function (off, theta) {
    var vadd = THREE.OFFLoader.Utils.vadd
    var vsub = THREE.OFFLoader.Utils.vsub
    var vdiv = THREE.OFFLoader.Utils.vdiv
    var vmul = THREE.OFFLoader.Utils.vmul
    var vdot = THREE.OFFLoader.Utils.vdot
    var vcross = THREE.OFFLoader.Utils.vcross
    var makenormal = THREE.OFFLoader.Utils.makenormal
    var vertices = off.vertices
    var faces = off.faces
    var m = new THREE.Matrix4;
    var face1; // First type 1 face, use this as model
    faces.forEach(function(face) {
        var vlist = face.vlist
        if (vlist.length >= 3) {
            // Save the normal for later
            face.normal = makenormal(vertices[vlist[0]],
                                     vertices[vlist[1]],
                                     vertices[vlist[2]])
            if (face.type == 0) {
                m.makeRotationAxis(face.normal,theta);
                vlist.forEach(function(v) {
                    vertices[v].applyMatrix4(m);
                })
            } else if (face.type == 1 && !face1) {
                face1 = face
            }
        }
    })
    // Now we have to find a factor k such that moving each type 0
    // face out by k will make each type 1 face original size (ie. with
    // side length 1).
    if (!face1) {
        alert("No secondary faces found in twister transform")
        return;
    }
    var normals = face1.vlist.map(function(v) {
        // Every vertex must be a vertex of primary face
        // Find the normal for that face
        for (var i = 0; i < faces.length; i++) {
            if (faces[i].type != 0) break;
            var vlist = faces[i].vlist
            for (var j = 0; j < vlist.length;j++) {
                if (vlist[j] == v) {
                    //console.log(i,faces[i].normal)
                    console.assert(faces[i].normal)
                    return faces[i].normal
                }
            }
        }
        console.assert(false)
    })
    var n0 = normals[0]
    var n1 = normals[1]
    var v0 = vertices[face1.vlist[0]]
    var v1 = vertices[face1.vlist[1]]
    //console.log(vsub(v1,v0).length())
    // Now find k such that |v0 + k*n0 - v1 - k*n1| = 1
    // ie. |(v0-v1) + k(n0-n1)| = 1
    // put V = v0-v1, N = n0-n1
    // then | V + kN| = 1 => V.V + 2kV.N + k^2N.N = 1
    // or: k^2N.N +k(2V.N) + (V.V-1) = 0
    // Quadratic equation with
    // A = N.N, B = 2V.N, C = V.V-1
    var V = vsub(v0,v1)
    var N = vsub(n0,n1)
    //console.log(n0,n1,N)
    var A = vdot(N,N)
    var B = 2*vdot(V,N)
    var C = vdot(V,V)-1
    var d = B*B - 4*A*C
    console.assert(d >= 0)
    if (d < 0 && false) {
        console.log(d,A,B,C)
        console.log("V",V,"N",N)
        console.log(v0,v1)
        console.log(n0,n1)
    }
    if (d < 0) d = 0;
    // This solution gives the right answer for theta = 0
    var k = (-B - Math.sqrt(d))/(2*A)
    //console.log(v0,v1)
    //console.log(n0,n1)
    //console.log(V,N)
    //console.log(A,B,C,d,k,A*k*k + B*k + C)
    //console.log(normals[0],normals[1],normals[2])
    // Now move each primary face out by k*facenormal
    // Require that each vertex belongs to exactly one
    // primary face.
    for (var i = 0; i < faces.length; i++) {
        if (faces[i].type != 0) break;
        var vlist = faces[i].vlist
        var t = vmul(faces[i].normal,k)
        for (var j = 0; j < vlist.length; j++) {
            vertices[vlist[j]].add(t)
        }
    }
}

THREE.OFFLoader.makeStar = function(N,M,theta) {
    var vector = THREE.OFFLoader.Utils.vector;
    N = N || 3
    M = M || 1
    var star = [];
    for (var i = 0; i < N; i++) {
        // Make the star vector length 1 and inclined at 45%
        // so we get a nice sinusoidal polar zono.
        var k = 1/Math.sqrt(2)
        var x = k*Math.cos(2*Math.PI*i*M/N);
        var z = k*Math.sin(2*Math.PI*i*M/N);
        var y = k*theta;
        var v = vector(x,y,z)
        star.push(v);
    }
    return star;
}

// Geo. Hart's excellent zonohedron algorithm
THREE.OFFLoader.starZonohedron = function(star,off,newstar) {
    var vadd = THREE.OFFLoader.Utils.vadd
    var vsub = THREE.OFFLoader.Utils.vsub
    var vdiv = THREE.OFFLoader.Utils.vdiv
    var vmul = THREE.OFFLoader.Utils.vmul
    var vdot = THREE.OFFLoader.Utils.vdot
    var vcross = THREE.OFFLoader.Utils.vcross
    var vnegate = THREE.OFFLoader.Utils.vnegate
    var vector = THREE.OFFLoader.Utils.vector
    var N = star.length;
    function apply(star,component) {
        var result = vector()
        console.assert(star.length == component.length);
        for (var i = 0; i < star.length; i++) {
            result.add(vmul(star[i],component[i]))
        }
        return result;
    }
    function addvector(p) {
        console.assert(newstar)
        // This linear scan doesn't really cut the mustard
        // but will do for now.
        for (var i = 0; i < newstar.length; i++) {
            if (vcross(newstar[i],p).length() < 1e-5) {
                return;
            }
        }
        //console.log("Adding",p)
        newstar.push(p);
    }
    function addvertex(p) {
        for (var i = 0; i < off.vertices.length; i++) {
            if (off.vertices[i].distanceTo(p) < 1e-4) {
                return i;
            }
        }
        i = off.vertices.length;
        off.vertices.push(p)
        return i;
    }
    function addface(fvs) {
        console.assert(off)
        var vlist = []
        for (var i = 0; i < fvs.length;i++) {
            vlist.push(addvertex(fvs[i]))
        }
        off.faces.push({ vlist:vlist })
    }
    // Each component is an N-vector of {-1,0,1}
    for (var sj = 1; sj < N; sj++) {
        for (var si = 0; si < sj; si++) {
            var face = []
            // Normal direction not critical here
            var normal = vcross(star[si],star[sj]);
            var edges = []
            for (var k = 0; k < N; k++) {
                var t = vdot(normal,star[k]);
                var eps = 1e-6;
                if (t < -eps) face[k] = -1;
                else if (t > eps) face[k] = 1;
                else {
                    face[k] = 0;
                    edges.push(k);
                }
            }
            // edges are indices of the face edges
            var K = edges.length;
            console.assert(K >= 2);
            if (edges[0] == si) { // Need to ensure uniqueness
                var npoints = 2*K;
                var p = []
                var centre = apply(star,face); // Also the face normal
                if (K == 2 && false) {
                    // easy case
                    //p[0] = centre - star[edges[0]] - star[edges[1]];
                    //p[1] = centre - star[edges[0]] + star[edges[1]];
                    p[0] = vsub(centre,vadd(star[edges[0]],star[edges[1]]));
                    p[1] = vsub(centre,vsub(star[edges[0]],star[edges[1]]));
                    p[2] = vadd(centre,vadd(star[edges[0]],star[edges[1]]));
                    p[3] = vadd(centre,vsub(star[edges[0]],star[edges[1]]));
                    //p[3] = centre + star[edges[0]] - star[edges[1]];
                    if (newstar) {
                        for (var j = 0; j < 4; j++) {
                            //if (p[j].y > obj.maxy) obj.maxy = p[j].y;
                            addvector(p[j]);
                        }
                    }
                    if (off) {
                        for (var i = 0; i < npoints; i++) {
                            var p0 = p[i];
                            var p1 = p[(i+1)%npoints];
                            addtriangle(centre, p0, p1);
                            addtriangle(vnegate(centre), vnegate(p1), vnegate(p0));
                        }
                    }
                } else {
                    p = [];
                    var fvs = []
                    for (var i = 0; i < K; i++) {
                        p.push(vcross(star[edges[i]], centre));
                    }
                    for (var i = 0; i < K; i++) {
                        var eps = 1e-4;
                        var pvec = vector();
                        for (var j = 0; j < K; j++) {
                            var t = vdot(p[j],star[edges[i]]);
                            if (t > eps) {
                                pvec.add(star[edges[j]]);
                            } else {
                                pvec.sub(star[edges[j]]);
                            }
                        }
                        // We could save a little work here by just adding
                        // one point, but this will do for now.
                        var p0 = vadd(centre,pvec);
                        var p1 = vsub(centre,pvec);
                        fvs.push(p0);
                        fvs.push(p1);
                    }
                    if (newstar) {
                        for (var i = 0; i < fvs.length; i++) {
                            addvector(fvs[i]);
                        }
                    }
                    if (off) {
                        // Need to sort the face vertices so they
                        // are in the right order.
                        var xaxis = vsub(fvs[0],centre);
                        var yaxis = vcross(xaxis,centre);
                        xaxis.normalize();
                        yaxis.normalize();
                        for (var i = 0; i < fvs.length; i++) {
                            // Nice JS trick - I can add an ad-hoc "angle" field in
                            // to the vectors in fvs & use that as a sorting key
                            var fv = vsub(fvs[i],centre)
                            fvs[i].angle = Math.atan2(vdot(fv,yaxis),vdot(fv,xaxis))
                        }
                        fvs.sort(function(v0,v1) { return v0.angle < v1.angle; })
                        addface(fvs);
                        // Opposite face - need to reverse vertex ordering to
                        // keep correct orientation
                        addface(fvs.map(vnegate).reverse())
                    }
                }
            }
        }
    }
    if (off) {
        console.log("StarZono:", off.faces.length, off.vertices.length);
    }
}
