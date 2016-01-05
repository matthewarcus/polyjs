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
        reversefaces: false
    }
};

THREE.OFFLoader.prototype = {
    constructor: THREE.OFFLoader,
    load: function ( url, options, onLoad, onProgress, onError ) {
	var scope = this;
	var loader = new THREE.XHRLoader( scope.manager );
	loader.setCrossOrigin( this.crossOrigin );
	loader.load( url, function ( text, other ) {
	    onLoad( scope.parse( text, options ) );
	}, onProgress, onError );
    },

    setCrossOrigin: function ( value ) {
	this.crossOrigin = value;
    },

    parse: function ( text, options ) {
	console.time( 'OFFLoader' );
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
        
	var geometry = new THREE.Geometry();
        var vertices = geometry.vertices;
        var faces = geometry.faces;
        var faceVertexUvs = geometry.faceVertexUvs;
        faceVertexUvs[0] = [];
        var state = 0, count = 0, nvertices = 0, nfaces = 0;
        function vmul(v0,theta) {
            v = v0.clone();
            v.multiplyScalar(theta);
            return v;
        }
        function vsub(v0, v1) {
            var t = v0.clone();
            t.sub(v1);
            return t;
        }
        function vadd(v0, v1) {
            var t = v0.clone();
            t.add(v1);
            return t;
        }
        function vcross(v0, v1) {
            var t = v0.clone();
            t.cross(v1);
            return t;
        }
        function vdot(v0, v1) {
            return v0.dot(v1);
        }
        function interp(v0,v1,theta) {
            var v = vmul(v0,theta);
            v.add(vmul(v1,1-theta));
            return v;
        }
        function makenormal(v0,v1,v2) {
            var n = vsub(v2,v0);
            n.cross(vsub(v1,v0));
            n.normalize();
            return n;
        }
        function perpdist(v,v0,v1) {
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
        }
        // Remove all occurences of a from s
        function remove(a,s) {
            while(true) {
                var i = s.indexOf(a);
                if (i < 0) break;
                s.splice(i,1);
            }
            return s;
        }
	function addface(a, b, c, uva, uvb, uvc, normal, color) {
	    faces.push(new THREE.Face3(a, b, c, normal, color));
            faceVertexUvs[0].push([uva,uvb,uvc]);
	}
        function addvertex(vertex) {
            var i = vertices.length;
            vertices.push(vertex);
            return i;
        }
        var vertexmodel;
        var edgemodel;
        var defaultvertexcolor = new THREE.Color(0.9,0.9,0.9);
        var defaultedgecolor = defaultvertexcolor;
        var m = new THREE.Matrix4;
        var tt = new THREE.Matrix4;
        var yaxis = new THREE.Vector3(0,1,0);
        function showvertex(v0, color) {
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
            var raxis = vcross(yaxis,v0);
            if (raxis.length() > 1e-4) {
                var theta = Math.acos(vdot(yaxis,v0)/vlen);
                raxis.normalize();
                tt.makeRotationAxis(raxis,theta);
                m.multiply(tt);
            }
            tt.makeTranslation(0,v0.length(),0);
            m.multiply(tt);
            for (var i = 0; i < vertexmodel.faces.length; i++) {
                vertexmodel.faces[i].color = color || defaultvertexcolor;
            }
            geometry.merge(vertexmodel,m);
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
            geometry.merge(edgemodel,m);
        }
        var cut = [];
        var cutlim = -2;
        var cutnormal = new THREE.Vector3(1,2,3);
        cutnormal.normalize();

	var lines = text.split( '\n' );
	for (var lineindex = 0; lineindex < lines.length; lineindex++) {
	    var line = lines[ lineindex ];
	    line = line.trim();
	    if ( line.length === 0 || line.charAt( 0 ) === '#' ) continue;
            var items = remove("",line.split(" "));
            if (state == 0) {
                if (items.length != 1 || items[0] != "OFF") {
                    alert("File load error: invalid OFF line: " + line);
                    return geometry;
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
                    // Mustn't call showvertex until after all vertices
                    // have been set up or the indexing goes awry.
                    if (allvertices) {
                        nvertices = vertices.length
                        for (var i = 0; i < nvertices; i++) {
                            showvertex(vertices[i]); // This extends vertices
                        }
                    }
                    state = 3;
                }
            } else if (state == 3) {
                if (nfaces > 0) {
                    nfaces--;
                    var color;
                    var n = parseInt(items[0]);
                    if (items.length > n+1) {
                        // Reuse color objects?
                        // Can a color just be an index?
                        color = new THREE.Color(parseFloat(items[n+1]),
                                                parseFloat(items[n+2]),
                                                parseFloat(items[n+3]));
                    }
                    var vlist = [];
                    for (var i = 0; i < n; i++) {
                        var v = parseInt(items[i+1]);
                        vlist.push(v);
                    }
                    if (n == 1) {
                        if (!novertices) {
                            var v0 = vertices[vlist[0]];
                            showvertex(v0,color);
                        }
                    } else if (n == 2) {
                        if (!noedges) {
                            var v0 = vertices[vlist[0]];
                            var v1 = vertices[vlist[1]];
                            showedge(v0,v1,color);
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
                            geometry.merge(geom);
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
                                if (reversefaces) addface(icentroid, p1, p0, uvcentroid, uv1, uv0, normal, color);
                                else addface(icentroid, p0, p1, uvcentroid, uv0, uv1, normal, color);
                            }
                        }
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
            faceVertexUvs[0] = []
            for (var i = 0; i < faces.length ; i++) {
                faceVertexUvs[0].push(
                    [new THREE.Vector2(0,0),
                     new THREE.Vector2(0,1),
                     new THREE.Vector2(1,1)])
            }
        }
        console.log("Vertices = " + vertices.length + " faces = " + faces.length);
	console.timeEnd( 'OFFLoader' );
	return geometry;
    }
};
