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

(function() {
    THREE.OFFLoader = function ( manager ) {
        this.manager = ( manager !== undefined ) ? manager : THREE.DefaultLoadingManager;
    };

    THREE.OFFLoader.prototype = {
        constructor: THREE.OFFLoader,
        load: function ( url, onLoad, onProgress, onError ) {
            var scope = this;
            var loader = new THREE.XHRLoader( scope.manager );
            //var loader = new THREE.XHRLoader();
            //loader.setCrossOrigin( this.crossOrigin ); // Not in r74?
            loader.load( url, function ( text, other ) {
                var off = scope.parse(text)
                console.assert(off)
                onLoad(off)
            }, onProgress, onError );
        },

        setCrossOrigin: function ( value ) {
            this.crossOrigin = value;
        },

        verbose: false,
    }

    THREE.OFFLoader.defaultoptions = function () {
        return {
            vertexstyle: "icosahedron",
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
        vtaxi: function (v0, v1) {
            return Math.abs(v0.x-v1.x) + Math.abs(v0.y-v1.y) + Math.abs(v0.z-v1.z)
        },
        vntaxi: function (v0, v1) {
            return Math.abs(v0.x+v1.x) + Math.abs(v0.y+v1.y) + Math.abs(v0.z+v1.z)
        },
        vnormalize: function (v) {
            var t = v.clone();
            t.normalize()
            return t
        },
        // Quaternion multiplication: r = pq.
        // r may be the same as p or q or both.
        qmul: function(p,q,r) {
            var p0 = p.x, p1 = p.y, p2 = p.z, p3 = p.w;
            var q0 = q.x, q1 = q.y, q2 = q.z, q3 = q.w;
            var r0 = p0*q0 - p1*q1 - p2*q2 - p3*q3;
            var r1 = p0*q1 + p2*q3 - p3*q2 + p1*q0;
            var r2 = p0*q2 + p3*q1 - p1*q3 + p2*q0;
            var r3 = p0*q3 + p1*q2 - p2*q1 + p3*q0;
            r.x = r0; r.y = r1; r.z = r2; r.w = r3;
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
            // This isn't a good idea if v0,v1,v2 are small
            // maybe normalize first.
            //if (n.length() < 1e-4) return null;
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
            purple: new THREE.Color(0.5,0.1,0.5),
            straw: new THREE.Color(1,1,0.5),
        }
    }

    THREE.OFFLoader.prototype.parse = function ( text, geometry, options ) {
        // Remove all occurences of a from s
        function remove(a,s) {
            while (true) {
                var i = s.indexOf(a);
                if (i < 0) break;
                s.splice(i,1);
            }
            return s;
        }
        if (this.verbose) console.time( 'OFFLoader' );
        var vertices = [], faces = [];
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
                    // eps can't be too big or texture interpolation goes awry.
                    // We could also treat at 4-space coordinates.
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
                    var n = parseInt(items[0]);
                    var face = {}
                    // Color can also just be an index
                    if (items.length == n+2) {
                        var type = parseInt(items[n+1]);
                        face.type = type;
                    } else if (items.length == n+4) {
                        // Reuse color objects?
                        var color = new THREE.Color(parseFloat(items[n+1]),
                                                    parseFloat(items[n+2]),
                                                    parseFloat(items[n+3]));
                        face.color = color;
                    }
                    var vlist = [];
                    for (var i = 0; i < n; i++) {
                        var v = parseInt(items[i+1]);
                        vlist.push(v);
                    }
                    face.vlist = vlist;
                    faces.push(face);
                }
            }
        }
        if (this.verbose) console.timeEnd( 'OFFLoader' );
        return { vertices: vertices, faces: faces }
    }


    var offMergeNormalMatrix = new THREE.Matrix3;
    function offMerge2(base, geometry, matrix) {
        if ( geometry instanceof THREE.Geometry === false ) {
            console.error( 'THREE.Geometry.offMerge(): geometry not an instance of THREE.Geometry.', geometry );
            return;
        }
        var normalMatrix = offMergeNormalMatrix,
            vertexOffset = base.nvertices,
            vertices1 = base.vertices,
            vertices2 = geometry.vertices,
            faceOffset = base.nfaces,
            faces1 = base.faces,
            faces2 = geometry.faces,
            uvs2 = geometry.faceVertexUvs[ 0 ];

        base.nvertices += vertices2.length
        base.nfaces += faces2.length
        
        if ( matrix !== undefined ) {
            normalMatrix.getNormalMatrix(matrix);
        }
        // vertices
        while (vertices1.length < vertexOffset+vertices2.length ) {
            vertices1.push(new THREE.Vector3);
        }
        for (var i = 0, il = vertices2.length; i < il; i++) {
            var vertex2 = vertices2[ i ];
            var vertex1 = vertices1[vertexOffset + i]
            vertex1.x = vertex2.x
            vertex1.y = vertex2.y
            vertex1.z = vertex2.z
            if ( matrix !== undefined ) vertex1.applyMatrix4( matrix );
        }
        // faces
        while (faces1.length < faceOffset + faces2.length) {
            faces1.push({ vlist: [0,0,0],
                          uvs: [new THREE.Vector2, new THREE.Vector2, new THREE.Vector2],
                          normals: [new THREE.Vector3, new THREE.Vector3, new THREE.Vector3],
                          colors: [0,0,0]
                        });
        }
        for ( i = 0, il = faces2.length; i < il; i ++ ) {
            var face2 = faces2[ i ], face1, color,
                face2VertexNormals = face2.vertexNormals,
                face2VertexColors = face2.vertexColors;
            face1 = faces1[faceOffset+i];
            //faceCopy = new THREE.Face3( face.a + vertexOffset, face.b + vertexOffset, face.c + vertexOffset );
            face1.vlist[0] = face2.a + vertexOffset;
            face1.vlist[1] = face2.b + vertexOffset;
            face1.vlist[2] = face2.c + vertexOffset;

            console.assert(face2VertexNormals.length == 0 || face2VertexNormals.length == 3)
            // If we have both a face normal and vertex normals, which to choose?
            // Probably should be an options, but for now, choose vertex normals.
            if (/*face2.normal ||*/ face2VertexNormals.length == 0) {
                face1.normals[0].copy(face2.normal);
                if ( normalMatrix !== undefined ) {
                    face1.normals[0].applyMatrix3( normalMatrix ).normalize();
                }
                face1.normals[1].copy(face1.normals[0]);
                face1.normals[2].copy(face1.normals[0]);
            } else {
                for ( var j = 0, jl = face2VertexNormals.length; j < jl; j ++ ) {
                    var normal1 = face1.normals[ j ];
                    var normal2 = face2VertexNormals[ j ];
                    normal1.copy(normal2)
                    if ( normalMatrix !== undefined ) {
                        normal1.applyMatrix3( normalMatrix ).normalize();
                    }
                }
            }
            // Color objects generally don't get modified after
            // construction so shallow copy (and if they do get modified,
            // we probably want to see the effect here anyway).
            if (face2VertexColors.length == 0) {
                face1.colors[0] = face1.colors[1] = face1.colors[2] = face2.color;
            } else {
                console.assert(face2VertexColors.length == 3);
                for ( var j = 0; j < 3; j++ ) {
                    face1.colors[j] = face2VertexColors[j];
                }
            }
        }
        for ( i = 0, il = uvs2.length; i < il; i++ ) {
            var uv1 = faces1[faceOffset+i].uvs;
            var uv2 = uvs2[ i ]
            for ( var j = 0, jl = uv2.length; j < jl; j ++ ) {
                uv1[j].x = uv2[j].x
                uv1[j].y = uv2[j].y
            }
        }
    }

    // The "geometry" here is just a vertices/faces pair, eg. could be our context.
    THREE.OFFLoader.display2 = function (off, geometry, options) {
        var vadd = THREE.OFFLoader.Utils.vadd
        var vsub = THREE.OFFLoader.Utils.vsub
        var vdiv = THREE.OFFLoader.Utils.vdiv
        var vdist = THREE.OFFLoader.Utils.vdist
        var vmul = THREE.OFFLoader.Utils.vmul
        var vdot = THREE.OFFLoader.Utils.vdot
        var vcross = THREE.OFFLoader.Utils.vcross
        var makenormal = THREE.OFFLoader.Utils.makenormal
        var Color = THREE.OFFLoader.Utils.Color
        var colors = [ Color.red, Color.yellow, Color.blue,
                       Color.green, Color.purple, Color.orange,
                       Color.cyan, Color.magenta, Color.white ]

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
        
        geometry.nvertices = 0
        geometry.nfaces = 0
        var vertices = geometry.vertices;
        var faces = geometry.faces;
        function addface(a, b, c, uva, uvb, uvc, normal, color) {
            if (geometry.nfaces == faces.length) {
                faces.push({ vlist: [0,0,0],
                             uvs: [new THREE.Vector2,new THREE.Vector2,new THREE.Vector2],
                             normals: [new THREE.Vector3, new THREE.Vector3, new THREE.Vector3],
                             colors: [0,0,0]
                           });
            }
            var face = faces[geometry.nfaces];
            geometry.nfaces++;
            face.vlist[0] = a
            face.vlist[1] = b
            face.vlist[2] = c
            face.normals[0].copy(normal);
            face.normals[1].copy(normal);
            face.normals[2].copy(normal);
            if (!color) color = Color.white;
            face.colors[0] = color;
            face.colors[1] = color;
            face.colors[2] = color;
            face.uvs[0].copy(uva);
            face.uvs[1].copy(uvb);
            face.uvs[2].copy(uvc);
        }
        function addvertex(vertex) {
            if (vertices.length <= geometry.nvertices) {
                vertices.push(new THREE.Vector3);
            }
            vertices[geometry.nvertices].copy(vertex);
            geometry.nvertices++;
            return geometry.nvertices-1;
        }
        var vertexmodel;
        var edgemodel;
        var defaultvertexcolor = new THREE.Color(0.6,0.6,0.6);
        var defaultedgecolor = new THREE.Color(0.9,0.9,0.9);
        var m = new THREE.Matrix4;
        var tt = new THREE.Matrix4;
        var yaxis = new THREE.Vector3(0,1,0);
        var xaxis = new THREE.Vector3(1,0,0);
        function showvertex(v0, color) {
            if (!vertexmodel) {
                if (vertexstyle === "cylinder") {
                    vertexmodel = new THREE.CylinderGeometry( vertexwidth, vertexwidth, vertexwidth );
                } else if (vertexstyle === "sphere") {
                    vertexmodel = new THREE.SphereGeometry( vertexwidth );
                } else {
                    if (vertexstyle !== "icosahedron") {
                        console.log("Unknown vertexstyle: " + vertexstyle);
                    }
                    vertexmodel = new THREE.IcosahedronGeometry( vertexwidth, 1 );
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
            offMerge2(geometry,vertexmodel,m);
        }
        var edge = new THREE.Vector3;  // Probably an absurd microptimization
        var midpoint = new THREE.Vector3; // Ditto.
        function showedge(v0, v1, color) {
            // Make cylinder, stretch to right length, rotate
            // to the correct orientation, then translate into position.
            edge.subVectors(v1,v0);
            var edgelen = edge.length();   // Total height of cylinder
            if (Math.abs(edgelen) < 1e-4) return;

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
            offMerge2(geometry,edgemodel,m);
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
            if (!color && off.faces[faceindex].type != undefined) {
                //console.log("facetype",off.faces[faceindex].type)
                color = colors[off.faces[faceindex].type%colors.length]
            }
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
                    offMerge2(geometry,geom);
                }
                if (!nofaces) {
                    var icentroid = addvertex(centroid);
                    var normal = new THREE.Vector3;
                    // Average the normals over the whole face.
                    // We could chance it and just check one pair of vertices.
                    for (var i = 0; i < n; i++) {
                        var v1 = vertices[vlist[i]];
                        var v2 = vertices[vlist[(i+1)%n]];
                        // If triangle degenerate, don't contribute to normal
                        var normal1 = makenormal(centroid,v2,v1);
                        if (normal1) {
                            normal.add(normal1); // Note order
                        }
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
                addface(icutcentre,cut[i],cut[i+1],color);
            }
        }
        if (this.verbose) {
            console.log("Vertices = " + geometry.nvertices + " faces = " + geometry.nfaces);
        }
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
            if (primarytype == undefined) primarytype = type
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
    THREE.OFFLoader.twistertransform = function (off, options) {
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
        var theta = options.theta || 0
        options.theta = theta + 0.01 //Increment for next step
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
        var vtaxi = THREE.OFFLoader.Utils.vtaxi
        var vntaxi = THREE.OFFLoader.Utils.vntaxi
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
            // This linear scan will do for now.
            for (var i = 0; i < newstar.length; i++) {
                // Cheap test catches most duplicates
                if (vtaxi(newstar[i],p) < 1e-5 || vntaxi(newstar[i],p) < 1e-5) {
                    return;
                }
            }
            for (var i = 0; i < newstar.length; i++) {
                if (vcross(newstar[i],p).length() < 1e-5) {
                    return;
                }
            }
            newstar.push(p);
        }
        var vertices, faces
        if (off) {
            vertices = off.vertices
            faces = off.faces
        }
        function addvertex(p) {
            if (0) {
                for (var i = 0; i < off.vertices.length; i++) {
                    if (off.vertices[i].distanceTo(p) < 1e-4) {
                        return i;
                    }
                }
            }
            i = off.vertices.length;
            off.vertices.push(p)
            return i;
        }
        var facetypes = []
        function getfacetype(vlist,centre) {
            //console.log(vlist)
            var sig = 0;
            // Add the areas of the segments of the face
            // Use Heron's formula
            // Since face is centrally symmetric, only
            // need to do half the edges
            // Could just calculate distance to each point once
            console.assert(vlist.length%2 == 0)
            for (var i = 0; i < vlist.length/2; i++) {
                var p0 = vlist[(i+0)%vlist.length]
                var p1 = vlist[(i+1)%vlist.length]
                //console.log(centre,p0,p1)
                var a = centre.distanceTo(p0)
                var b = centre.distanceTo(p1)
                var c = p0.distanceTo(p1);
                //console.log(a,b,c)
                var s = (a+b+c)/2;
                var area = Math.sqrt(s*(s-a)*(s-b)*(s-c))
                sig += area
            }
            for (var i = 0; i < facetypes.length; i++) {
                if (Math.abs(facetypes[i]-sig) < 1e-5) {
                    return i;
                }
            }
            facetypes.push(sig)
            return facetypes.length-1
        }
        function addface(fvs, facetype) {
            console.assert(off)
            var vlist = []
            for (var i = 0; i < fvs.length;i++) {
                vlist.push(addvertex(fvs[i]))
            }
            off.faces.push({ vlist:vlist, type: facetype })
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
                    var centre = apply(star,face); // Also the face normal
                    var p = []
                    for (var i = 0; i < K; i++) {
                        p.push(vcross(star[edges[i]], centre));
                    }
                    var fvs = []
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
                        var yaxis = vcross(centre,xaxis);
                        xaxis.normalize();
                        yaxis.normalize();
                        for (var i = 0; i < fvs.length; i++) {
                            // Nice JS trick - I can add an ad-hoc "angle" field in
                            // to the vectors in fvs & use that as a sorting key
                            var fv = vsub(fvs[i],centre)
                            fvs[i].angle = Math.atan2(vdot(fv,yaxis),vdot(fv,xaxis))
                        }
                        fvs = fvs.sort(function(v0,v1) {
                            // NB: Not just returning a boolean.
                            if (v0.angle < v1.angle) return -1;
                            else if (v0.angle > v1.angle) return 1;
                            else return 0;
                        })
                        var facetype = getfacetype(fvs,centre)
                        addface(fvs, facetype);
                        // Opposite face - need to reverse vertex ordering to
                        // keep correct orientation
                        addface(fvs.map(vnegate).reverse(), facetype)
                    }
                }
            }
        }
        if (off && false) {
            console.log("Facetypes:", facetypes)
            console.log("StarZono:", off.faces.length, off.vertices.length);
        }
        return off
    }
})()
