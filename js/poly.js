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

// TODO:
// *Snub duals
// *Store separate zrot value
// *Reset zrot (and other values)
// *7 point region drawing
// *Suppress face display
// *Add gradation from centre color mode
// *Sierpinski triangles
// *Reuse geometry when drawing compounds
// *Add check that we have a valid Schwarz triangle
// *Invert though midsphere (minimum of distance to midedge points)
// *Superimpose dual (needs midsphere inversion)
// *Exploding faces
// *Animated inversion
// *Better animated inversion
// *Step to next position
// *Explicitly specify trilinear coords
// *Get rid of global variables
// *Don't draw degenerate sectors of faces
// *Modularize
// *Use proper options object rather than globals
// *Texture coordinates
// *Reading OFF files
// *Invertible snub duals
// *Fix snubcentre for degenerate snub triangles
// *Report errors on texture & off file loading
// *Normal maps
// *Fix flicker from coplanar faces in polyhedra
// *Draw retroflex edges properly in polyhedra
// *Better mechanism to pass parameters to OFF functions
// *Use OFF vertices for compounding
// *BufferGeometry for polyhedra
// *BufferGeometry for OFF files
// *The Fourth Dimension
// !Sort out retroflex snub region edges bug
// 4-space rotations in the viewer
// 4-space OFF files
// Use THREE Vector3 etc. throughout.
// Fix flicker in OFF file geometries
// Generate OFF geometries for polyhedra
// Write out OFF/Ply files
// Sort out the mess of geometry reuse/cloning etc.
// Speed control
// Dimmer switch
// Bump maps
// Better lights
// UI to rotate object rather than scene
// Gesture UI
// Keyboard input for mobile
// Options for eg. drawing snub faces with centre or not
// Separate coplanar faces in explode.
// More efficient use of geometry objects
// Proper uvs for hemi faces & stellations
// Reuse point geometry between compounds
// Table driven options
// Better separation of generation and display
// Separate logical face structure from display structure
// Hollow faces
// Neo filling
// zrotation should be function of time
// Snub stellations
// Dual of stellations
// Stellation of duals
// More generic dual generation
// Stellation parameters
// Make URL for current scene
// More compact URL scheme
// Better colors
// Symmetric colors
// Use homogeneous coordinates

// Our state object.
function PolyContext(options) {
    // Options, modifiable during runtime
    this.verbose = false
    this.dostellate = false;
    this.dosnubify = false;
    this.drawtype = 0;
    this.invertinc = -0.01;
    this.ifact = 0;
    this.donormalize = false;
    this.docompound = false;
    this.dorotate = false;
    this.dozrotate = false;
    this.dochiral = false;
    this.colorstyle = 0;
    this.coloroffset = 0;
    this.hideface = [false,false,false,false]
    this.theta = 0;
    this.tridepth = 0;
    this.regionstyle = 0;
    this.explode = 0;
    this.dohemi = false;
    this.animstep = 10; // Seconds per animation step
    this.texturefile = null;
    this.normalfile = null;
    this.fnames = [];
    this.offfile = null;
    this.offoptions = THREE.OFFLoader.defaultoptions();
        
    // Take a walk around the Schwarz triangle
    this.tours = [
        [], // Can be set from query params
        [ [1,0,0],[1,1,0],[0,1,0],[0,1,1],[0,0,1],[1,0,1],[1,0,0],
          [1,1,1],[1,1,0],[0,1,0],[1,1,1],[0,1,1],[0,0,1],[1,1,1],[1,0,1] ],
        [
            [-1,1,0],[0,1,0],[1,0,0],[1,-1,0],
            [-1,0,1],[0,0,1],[1,0,0],[1,0,-1],
            [0,-1,1],[0,0,1],[0,1,0],[0,1,-1],
            [1,-1,-1],[1,0,0],[0,1,1],[-1,1,1],
            [-1,1,-1],[0,1,0],[1,0,1],[1,-1,1],
            [-1,-1,1],[0,0,1],[1,1,0],[1,1,-1],
            
        ], 
    ];
    this.tournum = 1;

    this.initrunning = false;
    this.initz = 3;
    this.initsym = [2,2,1]; // Single reflection
    this.angles  = [2,3,5];
    this.initt = 0;

    this.processoptions(options);
}

PolyContext.UpdateNone = 0;
PolyContext.UpdateView = 1;
PolyContext.UpdateColors = 2;
PolyContext.UpdateModel = 3;

PolyContext.Stopwatch = function(init,step,running) {
    // When running, reported time is difference between Date.now()
    // and stored time. When stopped, reported time is just the
    // stored time.
    this.running = running;
    this.scale = step*1000;  // Return time in seconds
    this.setTime(init || 0); // Different semantics in running and non-running states
}

PolyContext.Stopwatch.prototype.toggle = function() { 
    this.running = !this.running;
    this.time = Date.now()-this.time;
}
PolyContext.Stopwatch.prototype.setTime = function(t) { 
    if (this.running) this.time = Date.now()-t*this.scale;
    else this.time = t*this.scale;
}
PolyContext.Stopwatch.prototype.incTime = function(t) { 
    if (this.running) {
        this.time -= t*this.scale;
    } else {
        // Apply increment and step to next position
        this.time += t*this.scale;
        this.time = this.scale*(Math.floor(this.time/this.scale));
    }
}

PolyContext.Stopwatch.prototype.getTime = function() {
    if (this.running) return (Date.now()-this.time)/this.scale;
    else return this.time/this.scale;
}

PolyContext.prototype.makeinfostring = function() {
    // TBD: put more here
    var s = [];
    s.push("tri = [" + this.tri + "]");
    s.push("bary = [" + this.bary + "]");
    if (this.docompound) {
        s.push("theta = " + this.theta);
    }
    return s.join("; ");
}

PolyContext.prototype.generateURL = function() {
    // TBD: do this properly
    return this.url +
        (this.theta != 0 ? "&theta=" + this.theta : "") +
        "&colorstyle=" + this.colorstyle +
        // "&z=" + this.camera.position.z + // This won't work
        (this.tridepth ? "&tridepth=" + this.tridepth : "") +
        (this.docompound ? "&docompound" : "") +
        (this.dosnubify ? "&snub" : "") +
        (this.dostellate ? "&stellate" : "") +
        (this.dochiral ? "&chiral" : "") +
        (this.dohemi ? "&hemi" : "") +
        (this.drawtype == 1 ? "&both : this.drawtype == 2 ? "&dual : "");
}

// Given a time t, interpolate between the list of points
PolyContext.prototype.interpolatepoint = function(t,points) {
    t = t%points.length;
    if (t < 0) t += points.length;
    var s = t%1;
    t = Math.floor(t);
    var start = points[t];
    var end = points[(t+1)%points.length];
    var res = [];
    for (var i = 0; i < 3; i++) {
        res[i] = (1-s)*start[i] + s*end[i];
    }
    return res;
}

PolyContext.prototype.makeangles = function(a,b,c,d,e,f) {
    b = b || 1; d = d || 1; f = f || 1;
    a = Number(a); b = Number(b);
    c = Number(c); d = Number(d);
    e = Number(e); f = Number(f);
    // Check b/a + d/c + f/e > 1
    // ie. bce + dea + fac > ace
    if (a < b || c < d || e < f ||
        b*c*e + d*e*a + f*a*c <= a*c*e) {
        alert("Not a Schwarz Triangle");
    }
    return [a/b,c/d,e/f];
}

PolyContext.prototype.processoptions = function(options) {
    // TBD: A more compact parameter format would be good.
    var args = options.split('&');
    var context = this; // Need to see this in the forEach function.
    args.forEach(function (arg) {
        //console.log("Doing parameter '" + arg + "'");
        var matches;
        if (matches = arg.match(/^args=([\d]+)(?:\/([\d]+))?:([\d]+)(?:\/([\d]+))?:([\d]+)(?:\/([\d]+))?$/)) {
            context.angles = context.makeangles(matches[1],matches[2],matches[3],matches[4],matches[5],matches[6]);
            if (!context.url) context.url = "";
            context.url += matches[0];
        } else if (matches = arg.match(/^sym=([\d]+):([\d]+):([\d]+)$/)) {
            context.initsym = context.makeangles(matches[1],1,matches[2],1,matches[3],1);
            context.url += "&" + matches[0];
        } else if (matches = arg.match(/^tri=([\d-.]+):([\d-.]+):([\d-.]+)$/)) {
            context.tours[0].push([Number(matches[1]),Number(matches[2]),Number(matches[3])]);
            context.tournum = 0;
        } else if (matches = arg.match(/^theta=([\d]+)$/)) {
            context.theta = Math.PI/Number(matches[1]);
        } else if (matches = arg.match(/^theta=([\d]+.[\d]+)$/)) {
            context.theta = Number(matches[1]);
        } else if (matches = arg.match(/^colorstyle=([\d]+)$/)) {
            context.colorstyle = Number(matches[1]);
        } else if (matches = arg.match(/^regionstyle=([\d]+)$/)) {
            context.regionstyle = Number(matches[1]);
        } else if (matches = arg.match(/^hide=([\d]+)$/)) {
            context.hideface[Number(matches[1])-1] = true;
        } else if (matches = arg.match(/^t=([\d.]+)$/)) {
            context.initt = Number(matches[1]);
        } else if (matches = arg.match(/^z=([\d.]+)$/)) {
            context.initz = Number(matches[1]);
        } else if (matches = arg.match(/^radius=([\d.]+)$/)) {
            context.radius = Number(matches[1]);
        } else if (matches = arg.match(/^explode=([\d.]+)$/)) {
            context.explode = Number(matches[1]);
        } else if (matches = arg.match(/^tridepth=([\d]+)$/)) {
            context.tridepth = Number(matches[1]);
        } else if (matches = arg.match(/^tour=([\d]+)$/)) {
            context.tournum = Number(matches[1]);
        } else if (matches = arg.match(/^compound$/)) {
            context.docompound = true;
        } else if (matches = arg.match(/^snub$/)) {
            context.dosnubify = true;
        } else if (matches = arg.match(/^stellate$/)) {
            context.dostellate = true;
        } else if (matches = arg.match(/^both$/)) {
            context.drawtype = 1;
        } else if (matches = arg.match(/^dual$/)) {
            context.drawtype = 2;
        } else if (matches = arg.match(/^invert$/)) {
            context.invertinc = -context.invertinc;
            context.ifact = 1;
        } else if (matches = arg.match(/^zrotate$/)) {
            context.dozrotate = true;
        } else if (matches = arg.match(/^chiral$/)) {
            context.dochiral = true;
        } else if (matches = arg.match(/^rotate$/)) {
            context.dorotate = true;
        } else if (matches = arg.match(/^hemi$/)) {
            context.dohemi = true;
        } else if (matches = arg.match(/^run$/)) {
            context.initrunning = true;
        } else if (matches = arg.match(/^texture=(.+)$/)) {
            context.texturefile = matches[1];
        } else if (matches = arg.match(/^normal=(.+)$/)) {
            context.normalfile = matches[1];
        } else if (matches = arg.match(/^function=(.+)$/)) {
            context.fnames.push(matches[1]);
        } else if (matches = arg.match(/^off=(.+)$/)) {
            context.offfile = matches[1];
        } else if (matches = arg.match(/^vertexstyle=(.+)$/)) {
            context.offoptions.vertexstyle = matches[1];
        } else if (matches = arg.match(/^vertexwidth=([\d.]+)$/)) {
            context.offoptions.vertexwidth = Number(matches[1]);
        } else if (matches = arg.match(/^edgewidth=([\d.]+)$/)) {
            context.offoptions.edgewidth = Number(matches[1]);
        } else if (matches = arg.match(/^off.([^=]+)=(.+)$/)) {
            //console.log("Got", matches[1], matches[2])
            context.offoptions[matches[1]] = matches[2];
        } else if (matches = arg.match(/^allvertices$/)) {
            context.offoptions.allvertices = true;
        } else if (matches = arg.match(/^alledges$/)) {
            context.offoptions.alledges = true;
        } else if (matches = arg.match(/^novertices$/)) {
            context.offoptions.novertices = true;
        } else if (matches = arg.match(/^noedges$/)) {
            context.offoptions.noedges = true;
        } else if (matches = arg.match(/^nofaces$/)) {
            context.offoptions.nofaces = true;
        } else if (matches = arg.match(/^n=([\d.]+)$/)) {
            context.offoptions.n = Number(matches[1]);
        } else {
            console.log("Ignoring parameter '" + arg + "'");
        }
    });
    if (this.tours[0].length == 0) {
        this.tours[0].push([1,1,1]);
    }
}

PolyContext.prototype.handleKey = function(key) {
    var context = this;
    var handled = true;
    var update = PolyContext.UpdateModel; // Assume the worst
    switch(key) {
    case ' ':
        this.stopwatch.toggle();
        break;
    case 'z':
        this.dozrotate = !this.dozrotate;
        break;
    case 'a':
        this.theta = 0;
        break;
    case '[':
        this.stopwatch.incTime(-1);
        break;
    case ']':
        this.stopwatch.incTime(1);
        break;
    case 'p':
        this.stopwatch.setTime(0);
        break;
    case ',':
        //this.stopwatch.scale *= 1.05;
        break;
    case ',':
        //this.stopwatch.scale /= 1.05;
        break;
    case 'r':
        this.dorotate = !this.dorotate;
        update = PolyContext.UpdateView;
        break;
    case 'h':
        this.dohemi = !this.dohemi;
        break;
    case 's':
        this.dosnubify = !this.dosnubify;
        break;
    case 'q':
        this.dostellate = !this.dostellate;
        break;
    case 'e':
        this.explode += 1/16;
        break;
    case 'w':
        this.explode -= 1/16;
        break;
    case 'd':
        this.drawtype = (this.drawtype+1)%3
        break;
    case 'n':
        this.donormalize = !this.donormalize;
        break;
    case 'N':
        this.shownormals = !this.shownormals;
        break;
    case 'i':
        this.invertinc = -this.invertinc;
        break;
    case 'f':
        this.colorstyle = (this.colorstyle+1) % this.numcolorstyles
        update = PolyContext.UpdateColors;
        break;
    case 'T':
        if (context.texture) {
            if (context.material.map) context.material.map = null;
            else context.material.map = context.texture;
            context.material.needsUpdate = true;
            context.needclone = true;
        } else {
            update = PolyContext.UpdateNone;
        }
        break;
    case 't':
        this.tournum = (this.tournum+1)%this.tours.length;
        this.stopwatch.setTime(0);
        break;
    case 'U':
        update = PolyContext.UpdateNone;
        alert(this.generateURL());
        break;
    case 'u':
        this.coloroffset += 1;
        update = PolyContext.UpdateColors;
        break;
    case 'V':
        this.verbose = !this.verbose
        if (this.verbose) console.log("Verbose on")
        update = PolyContext.updateNone
        break
    case 'v':
        this.coloroffset-1;
        update = PolyContext.UpdateColors;
        break;
    case 'c':
        this.docompound = !this.docompound;
        this.needclone = true;
        break;
    case 'x':
        this.dochiral = !this.dochiral;
        break;
    case '=':
        this.tridepth = this.tridepth+1;
        break;
    case '-':
        if (this.tridepth > 0) this.tridepth = this.tridepth-1;
        break;
    case 'y':
        this.regionstyle = (this.regionstyle+1)%2;
        break;
    case '?':
        update = PolyContext.UpdateNone;
        alert(this.makeinfostring());
        break;
    case '1':
        this.hideface[0] = !this.hideface[0];
        break;
    case '2':
        this.hideface[1] = !this.hideface[1];
        break;
    case '3':
        this.hideface[2] = !this.hideface[2];
        break;
    case '4':
        this.hideface[3] = !this.hideface[3];
        break;
    case '9':
        this.offoptions.w *= 1.1;
        this.needclone = true;
        break;
    case '0':
        this.offoptions.w /= 1.1;
        this.needclone = true;
        break;
    default:
        handled = false;
    }
    if (handled) context.update(update);
    return handled;
}

// Finally, draw a point (via. the context.drawpoint0 function).
// This is sort of time-critical so inline a lot of stuff.
// Currently, we call this function several times for each actual point
// (unless we are doing an exploded view).
// so there is some significant optimizations possible here.
// Return the index in the geometry's list of points.
PolyContext.prototype.drawpoint = function(p,offset,facefact) {
    var Vector = Geometry.Vector;
    var w = p[3] || 1; // Homogeneous coords
    // For fun, do the explode before the inversion
    if (this.donormalize) w = Vector.length(p);
    var ifact = this.ifact;
    if (ifact != 0) {
        //var q = Vector.invert(p,this.midsphere);
        var k = w*Vector.dot(p,p)/this.midsphere;
        var w0 = (1-ifact)*w+ifact*k;
        //console.log(ifact,k,w,w0);
        w = w0;
    }
    var x = p[0]/w, y = p[1]/w, z = p[2]/w;
    var explode = this.explode;
    if (offset) {
        if (explode != 0) {
            x += explode*offset[0];
            y += explode*offset[1];
            z += explode*offset[2];
        }
    }
    if (facefact) {
        x *= (1+facefact);
        y *= (1+facefact);
        z *= (1+facefact);
    }
    if (this.transmat) {
        var m = this.transmat;
        var x1 = m[0][0]*x + m[0][1]*y + m[0][2]*z;
        var y1 = m[1][0]*x + m[1][1]*y + m[1][2]*z;
        var z1 = m[2][0]*x + m[2][1]*y + m[2][2]*z;
        x = x1; y = y1; z = z1;
    }
    return this.drawpoint0(x,y,z);
}

// TBD: encapsulate the various drawing options
// For n > 0, draw a Sierpinski triangle (or some
// variation thereof).
PolyContext.prototype.drawtriangle = function(p,q,r,uvs,offset,type,i,n) {
    var Vector = Geometry.Vector;
    if (n == 0) {
        var facefact = (this.compound/*+type*/)*0.0001 // Reduce z-fighting
        var index0 = this.drawpoint(p,offset,facefact);
        var index1 = this.drawpoint(q,offset,facefact);
        var index2 = this.drawpoint(r,offset,facefact);
        this.drawtriangle0(index0,index1,index2,uvs,type,i);
    } else {
        // Vary ttt for slightly tacky rotating effect.
        var ttt = 0.5;
        var p1 = Vector.interp(p,q,ttt);
        var q1 = Vector.interp(q,r,ttt);
        var r1 = Vector.interp(r,p,ttt);
        // TBD: Different triangulation styles here
        // k = 0,1,2,3,n%4,random etc.
        //var k = n%4;
        var k = 3;
        // TBD: texture coordinates!
        if (k != 0) this.drawtriangle(p,p1,r1,uvs,offset,type,i,n-1);
        if (k != 1) this.drawtriangle(p1,q,q1,uvs,offset,type,i,n-1);
        if (k != 2) this.drawtriangle(r1,q1,r,uvs,offset,type,i,n-1);
        if (k != 3) this.drawtriangle(p1,q1,r1,uvs,offset,type,i,n-1);
    }
}

PolyContext.prototype.drawface = function(centre,plist,facetype,index,tridepth,parity) {
    function mod(i,n) {
        i %= n;
        if (i < 0) i += n;
        return i;
    }
    var Vector = Geometry.Vector;
    var context = this;
    if (this.drawface0) {
        this.drawface0(plist); // Call alternative function, if defined
    } else {
        var needuvs = context.texturefile;
        if (this.dohemi) {
            // Do this before we do any modifications to plist for retroflex edges.
            // TBD: Proper texture coordinates!
            if (needuvs) {
                var uvs = [new THREE.Vector2(0,0),
                       new THREE.Vector2(1,0),
                       new THREE.Vector2(0,1)];
            }
            for (var i = 0; i < plist.length; i++) {
                var p1 = plist[i];
                var p2 = plist[(i+1)%plist.length];
                if (Vector.taxi(p1,p2) > 1e-4) {
                    this.drawtriangle([0,0,0],p1,p2,uvs,
                                      centre,4,index,tridepth);
                }
            }
        }

        var faceuvs = [];
        if (needuvs) {
            // We really ought to only calculate uvs when we actually want to texture the object.
            // Work out texture coordinates - these should be coherent across the face
            // Firstly, get orthogonal axes in face plane.
            // We use the first point, relative to the centre as the u-axis
            var uaxis;
            for (var i = 0; i < plist.length; i++) {
                // Some points may coincide with the centre so try and find
                // one that doesn't. If we fail, the face is degenerate anyway
                if (Vector.taxi(plist[i],centre) > 1e-4) {
                    uaxis = Vector.normalize(Vector.sub(plist[i],centre));
                    break;
                }
            }
            //console.assert(uaxis,"No suitable u axis for texture coordinate generation");
            if (!uaxis) uaxis = [1,0,0] // eg. if the face is degenerate
            var vaxis = Vector.normalize(Vector.cross(centre,uaxis)); // centre is normal
            if (parity) vaxis = Vector.negate(vaxis); // Mirror image

            // Compute uvs for face points
            for (var i = 0; i < plist.length; i++) {
                var p = plist[i]
                var u = 0.5 + Vector.dot(Vector.sub(p,centre),uaxis);
                var v = 0.5 + Vector.dot(Vector.sub(p,centre),vaxis);
                faceuvs.push(new THREE.Vector2(u,v));
            }
            var uvcentre = new THREE.Vector2(0.5,0.5);
        }
        if (plist.length == 3) {
            // Should do this when some edges are degenerate too.
            this.drawtriangle(plist[0],plist[1],plist[2],
                              faceuvs,centre,facetype,index,tridepth);
        } else {
            // Triangulation phase
            // Now we must work out if we have any retroflex edges
            // Since the face is "semiregular" we just look at first two.
            if (this.doretroflex && facetype < 4 && plist.length > 4) {
                // Only do this for regular poly faces & only if more
                // than 4 sides (otherwise normal drawing out from the
                // centre is OK).
                // We don't need to do this for every vertex - can do
                // calculation once for all the faces of a particular type.
                var newplist = []
                var newfaceuvs = []
                for (var i = 0; i < plist.length; i++) {
                    // Look for a retroflex edge: adjacent edges intersect
                    // edge side of centre, so use that as new vertex and
                    // draw the edge as a sort of ear.
                    var p0 = plist[mod(i-1,plist.length)] // previous vertex
                    var p1 = plist[i]
                    var p2 = plist[mod(i+1,plist.length)] // Edge vertices
                    var c = Vector.mid(p1,p2) // Edge centre
                    var e = Vector.sub(p2,p1) // Edge vector
                    var v = Vector.sub(p1,p0) // Previous edge
                    var tmp = Vector.dot(v,e);
                    if (Math.abs(tmp) < 1e-4) break; // Avoid divide by zero & infinities.
                    var k = -Vector.dot(p0,e)/tmp;
                    var q = Vector.add(p0,Vector.mul(v,k)) // Intersection with centre line
                    // q = kc => q.c = k*c.c => k = q.c/c.c
                    var c0 = Vector.sub(c,centre)
                    var q0 = Vector.sub(q,centre)
                    var efact = Vector.dot(q0,c0)/Vector.dot(c0,c0)
                    if (0 < efact && efact < 1) {
                        if (needuvs) {
                            // q is new point to draw triangle to edge from
                            var u = 0.5 + Vector.dot(Vector.sub(q,centre),uaxis);
                            var v = 0.5 + Vector.dot(Vector.sub(q,centre),vaxis);
                            var newuv = new THREE.Vector2(u,v);
                            newfaceuvs.push(newuv);
                            var uvs = [newuv,faceuvs[i],faceuvs[(i+1)%plist.length]];
                        }
                        newplist.push(q)
                        this.drawtriangle(q,p1,p2,uvs,
                                          centre,facetype,index,tridepth);
                    }
                }
                if (newplist.length > 0) {
                    plist = newplist;
                    faceuvs = newfaceuvs;
                } else {
                    // No retroflex needed for other faces.
                    this.doretroflex = false;
                }
            }
            console.assert(faceuvs.length == 0 || faceuvs.length == plist.length)
            for (var i = 0; i < plist.length; i++) {
                var p1 = plist[i];
                var p2 = plist[(i+1)%plist.length];
                if (needuvs) {
                    var uvs = [uvcentre,faceuvs[i],faceuvs[(i+1)%plist.length]];
                }
                if (Vector.taxi(p1,p2) > 1e-4) { // Don't draw if degenerate
                    this.drawtriangle(centre,p1,p2,uvs,
                                      centre,facetype,index,tridepth);
                }
            }
        }
    }
}

// Draw the Schwarz triangles (fundamental regions) as triangles.
PolyContext.prototype.drawregions = function() {
    var Vector = Geometry.Vector;
    var schwarz = this.schwarz;
    var facedata = this.facedata;
    var points = schwarz.points;
    var regions = schwarz.regions;
    var adjacent = schwarz.adjacent;
    var regionpoints = facedata.regionpoints;
    for (var i = 0; i < regions.length; i++) {
        var region = regions[i];
        var plist = [];
        if (region[3] == 0 || !this.dosnubify) {
            // The method here is to collect certain points on the normal
            // polyhedral figure and invert them in the midsphere (centre point
            // of edges).
            var snubcentre = facedata.snubcentre;
            for (var j = 0; j < 3; j++) {
                plist.push(Vector.mul(points[region[j]], facedata.facedistances[j]));
                if (this.dosnubify) {
                    // This isn't the 'regionpoint' but the 'snubpoint' ie. the
                    // centre of the snub triangle for the face.
                    var adj1 = adjacent[i][(j+2)%3];
                    var p1 = schwarz.applybary(snubcentre,adj1);
                    if (this.regionstyle == 1) {
                        plist.push(p1);
                    } else if (this.regionstyle == 0) {
                        // Include edge points where we cross with the
                        // snub triangle.
                        var adj0 = adjacent[adj1][(j+1)%3];
                        var p0 = regionpoints[adj0];
                        var adj2 = adjacent[adj1][j];
                        var p2 = regionpoints[adj2];
                        plist.push(Vector.mid(regionpoints[i],p0));
                        plist.push(p1);
                        plist.push(Vector.mid(regionpoints[i],p2));
                    } else {
                        // NB: console.assert is a Mozilla feature & non-portable.
                        console.assert(false,"Unknown regionstyle: " + regionstyle);
                    }
                } else if (this.regionstyle == 0) {
                    // Add in edge points for "correct" inversion
                    var edgecentre = schwarz.applybary(facedata.edgecentres[j],i);
                    plist.push(edgecentre);
                }
            }
            for (var j = 0; j < plist.length; j++) {
                plist[j] = Vector.invert(plist[j],this.midsphere);
            }
            var facetype = 4+region[3];
            var centre = Vector.invert(regionpoints[i],this.midsphere);
            var parity = region[3];
            this.drawface(centre,plist,facetype,i,this.tridepth,parity);
        }
    }
}

PolyContext.prototype.drawfaces = function() {
    var Vector = Geometry.Vector;
    var schwarz = this.schwarz;
    var facedata = this.facedata;
    var points = schwarz.points;
    var regions = schwarz.regions;
    var regionpoints = facedata.regionpoints;
    if (this.dostellate && this.texturefile) {
        // TBD: do this properly!
        var uvs = [new THREE.Vector2(0,0),
                   new THREE.Vector2(1,0),
                   new THREE.Vector2(0,1)];
    }
    // Now go through the facepoint lists for each type of face
    for (var type = 0; type < 3; type++) {
        //console.log(facedata.faces[type],type);
        if (!facedata.faces[type] || this.hideface[type]) continue;
        var facepoints = schwarz.faces[type];
        if (this.dostellate) {
            // Draw stellated face.
            for (var i = 0; i < facepoints.length; i++) {
                var plist = facepoints[i];
                if (plist) {
                    // Add a centre vertex, this is a corner of a Schwarz triangle.
                    // facedistances[type] tells us how far it is from the origin (the face
                    // centres, ie. the region points are at distance 1).
                    var centre = Vector.mul(points[i],facedata.facedistances[type])                    
                    var facelines = facedata.facelines[type];
                    if (!facelines) continue;
                    for (var j = 0; j < facelines.length; j++) {
                        // Triangle from centre to two vertices on plist
                        var p1 = schwarz.applybary(facelines[j][0],plist[0]);
                        var p2 = schwarz.applybary(facelines[j][1],plist[0]);
                        this.drawtriangle(centre,p1,p2,uvs,centre,type,i,this.tridepth);
                    }
                }
            }
        } else {
            this.doretroflex = true; // Gets set to false if we don't need it
            for (var i = 0; i < facepoints.length; i++) {
                var fpoints = facepoints[i];
                if (fpoints && (!this.dosnubify || fpoints.length > 4)) {
                    // Add a centre vertex, this is a corner of a Schwarz triangle.
                    // facedistances[type] tells us how far it is from the origin (the face
                    // centres, ie. the region points are at distance 1.
                    var centre = Vector.mul(points[i],facedata.facedistances[type]);
                    // For a snub polyhedron, we only include even (+ve) regions.
                    var inc = this.dosnubify ? 2 : 1;
                    var plist = [];
                    for (var j = 0; j < fpoints.length; j+=inc) {
                        plist.push(regionpoints[fpoints[j]]);
                    }
                    this.drawface(centre,plist,type,i,this.tridepth);
                }
            }
        }
    }
    if (this.dosnubify && !this.dostellate && !this.hideface[3] && facedata.snubcentre) {
        // Now draw a snub triangle for -ve regions
        for (var i = 0; i < regions.length; i++) {
            var t = regions[i];
            if (t[3] == 1) {
                var centre = schwarz.applybary(facedata.snubcentre, i);
                var plist = [regionpoints[schwarz.adjacent[i][0]],
                             regionpoints[schwarz.adjacent[i][1]],
                             regionpoints[schwarz.adjacent[i][2]]];
                this.drawface(centre,plist,3,i,this.tridepth);
            }
        }
    }
}

// Set up context for drawing a polyhedron in the geometry object,
// based on the fundamental regions info in schwarz, using
// trilinear coords tri, and with the colorface coloring function.
PolyContext.prototype.setup = function(tri) {
    var Vector = Geometry.Vector;
    var PointSet = Geometry.PointSet;
    this.npoints = 0;
    this.nfaces = 0;
    console.assert(this.pointset == null);
    console.assert(this.transmat == null);
    var schwarz = this.schwarz;
    var regions = schwarz.regions;
    var points = schwarz.points;
    // Convert trilinear coords to barycentric.
    if (!this.dosnubify) {
        var bary = schwarz.tri2bary(tri);
    } else {
        // For snub polyhedra, we use a bary coord so
        // that 1,1,1 maps through to snub point.
        var bary = [tri[0]*schwarz.snuba[0],
                    tri[1]*schwarz.snuba[1],
                    tri[2]*schwarz.snuba[2]];
        // Scale to unit sphere
        var s = schwarz.applybary(bary,0);
        bary = Vector.div(bary,Vector.length(s));
    }
    var facedata = schwarz.makefacedata(bary);
    if (this.docompound) {
        var regionpoints = facedata.regionpoints;
        // Set up the initial pointset (for compound generation).
        var pointset = [];
        for (var i = 0; i < regionpoints.length; i++) {
            if (!this.dosnubify || regions[i][3] == 0) {
                PointSet.add(regionpoints[i],pointset,1e-6);
            }
        }
        PointSet.sort(pointset);
        this.pointset = pointset;
    }
    this.tri = tri
    this.bary = bary
    this.facedata = facedata;
    if (this.dosnubify) {
        this.midsphere = facedata.snubsphere;
    } else {
        this.midsphere = facedata.midsphere;
    }
    if (this.dostellate) {
        schwarz.stellate(facedata, this.hideface, this.radius);
    }
}

PolyContext.prototype.drawpolyhedron = function() {
    // Do the actual drawing
    if (this.drawtype <= 1) this.drawfaces();
    if (this.drawtype >= 1) this.drawregions();
}

PolyContext.prototype.drawcompound = function(drawpolyhedron) {
    var Vector = Geometry.Vector;
    var PointSet = Geometry.PointSet;
    // We use a group word to describe a particular transformation
    var gens = ["P","Q","R"];
    var rgens = ["PQ","QR","RP"];
    this.compound = 0;
    this.transmat = null;
    this.pointset = null;
    drawpolyhedron(true,"");
    if (this.docompound) {
        console.assert(this.pointset,"No pointset defined for compounding");
        // Set up for compounding
        var ip = Vector.zrot(this.symmetry[0], this.theta);
        var iq = Vector.zrot(this.symmetry[1], this.theta);
        var ir = Vector.zrot(this.symmetry[2], this.theta);
        // Now we do a transitive closure operation. Keep a list of sets of
        // points, pull sets out and apply all basic operations (from the
        // symmetry group), add back in any sets we haven't previously seen.
        // Initially we just know about the basic pointset just generated.
        var pointsets = [{trans: "", pointset: this.pointset}];
        this.compound = 1;
        for (var pindex = 0, count = 0;
             pindex != pointsets.length && count < 400;
             pindex++, count++) {
            var entry = pointsets[pindex];
            var ops = this.dochiral ? rgens : gens; // Rotation only or all operations?
            for (var i = 0; i < ops.length; i++) {
                //console.log(ops[i] + " " + entry.trans);
                //if (entry.trans.length > 0 && entry.trans[entry.trans.length-1] == ops[i][0]) { continue; }
                var newset = Geometry.mapply(ops[i],entry.pointset,ip,iq,ir);
                PointSet.sort(newset);
                // Have we seen this pointset before?
                var seen = false;
                var parity = (entry.trans.length + ops[i].length)%2;
                for (var j = 0; j < pointsets.length; j++) {
                    // If we care about chirality (ie. are generating a snub figure)
                    // then don't consider reflected pointsets the same).
                    if (this.dosnubify && pointsets[j].trans.length%2 != parity) continue;
                    if (PointSet.equal(newset,pointsets[j].pointset,1e-4)) {
                        seen = true;
                        break;
                    }
                }
                //console.log(ops[i] + " " + entry.trans + " " + seen);
                if (!seen) {
                    // A new compound element. Draw it properly and add to queue.
                    // TBD: we could just reuse the geometry with a rotation
                    // and a change of color.
                    var newtrans = entry.trans + ops[i];
                    this.transmat = Geometry.makematrix(newtrans,ip,iq,ir);
                    //console.log(JSON.stringify(this.transmat));
                    //console.log("Compound " + this.compound + " " + newtrans);
                    drawpolyhedron(false);
                    pointsets.push({trans: newtrans, pointset: newset});
                    this.compound++;
                }
            }
        }
        //console.log("Compounds: " + this.compound);
    }
}

function offload(file, context) {
    if (context.verbose) console.log("Loading " + file);
    var manager = new THREE.LoadingManager();
    manager.onProgress = function ( item, loaded, total ) {
        //console.log( item, loaded, total );
    };
    var onProgress = function ( xhr ) {
        if ( xhr.lengthComputable ) {
            var percentComplete = xhr.loaded / xhr.total * 100;
            if (context.verbose) console.log( Math.round(percentComplete, 2) + '% downloaded' );
        }
    };

    var onError = function ( xhr ) {
        alert("Error loading " + file);
        console.log("Error loading " + file);
        console.log(xhr);
    };
    var loader = new THREE.OFFLoader( manager );
    loader.load( file+"?"+new Date().getTime(),
                 function (off) {
                     context.offdata = off
                     context.needclone = true
                     context.update(PolyContext.UpdateModel);
                 },
                 onProgress, onError);
}

    var compoundColors = [
        new THREE.Color(0xff0000),new THREE.Color(0x00ff00),new THREE.Color(0x0000ff),
        new THREE.Color(0xffff00),new THREE.Color(0xffffff),new THREE.Color(0x00ffff),
        new THREE.Color(0x808080),new THREE.Color(0xff00ff)
    ];

PolyContext.prototype.updategeometry = function(geometry,computenormals,colorstyle) {
    var context = this;
    // Pull this out to separate function
    var needclone = !context.vertexarray || context.vertexarray.length < context.nfaces*3*3;
    var needuvs = context.texturefile != undefined;
    if (needclone && context.nfaces > 0) {
        context.vertexarray = new Float32Array(context.nfaces*3*3);
        context.colorarray = new Float32Array(context.nfaces*3*4);
        context.normalarray = new Float32Array(context.nfaces*3*3);
        if (needuvs) context.uvarray = new Float32Array(context.nfaces*3*2);
    }
    var vertexarray = context.vertexarray;
    var colorarray = context.colorarray;
    var normalarray = context.normalarray;
    if (needuvs) var uvarray = context.uvarray;
    var ncolors = compoundColors.length;
    var white = new THREE.Color(0xffffff);
    for (var i = 0; i < context.nfaces; i++) {
        var face = context.faces[i];
        for (var j = 0; j < 3; j++) {
            if (colorstyle === 0) {
                var color = face.colors[j];
            } else if (colorstyle === 1) {
                var color = white;
            } else {
                if (colorstyle == 2) {
                    var compound = Math.floor(i/context.basefaces);
                } else {
                    var compound = i%ncolors;
                }
                //console.log(compound);
                var color = compoundColors[(context.coloroffset + compound)%ncolors];
            }
            var vertex = context.vertices[face.vlist[j]];
            vertexarray[i*3*3+j*3+0] = vertex.x;
            vertexarray[i*3*3+j*3+1] = vertex.y;
            vertexarray[i*3*3+j*3+2] = vertex.z;
            colorarray[i*3*4+j*4+0] = color.r;
            colorarray[i*3*4+j*4+1] = color.g;
            colorarray[i*3*4+j*4+2] = color.b;
            colorarray[i*3*4+j*4+3] = color.a;
            if (needuvs) {
                var uvs = face.uvs[j];
                //console.log(uvs.x, uvs.y);
                uvarray[i*3*2+j*2+0] = uvs.x;
                uvarray[i*3*2+j*2+1] = uvs.y;
            }
            if (!computenormals) {
                var normal = face.normals[j];
                normalarray[i*3*3+j*3+0] = normal.x;
                normalarray[i*3*3+j*3+1] = normal.y;
                normalarray[i*3*3+j*3+2] = normal.z;
            }
        }
    }
    if (needclone && vertexarray) {
        console.log("Cloning geometry",geometry.uuid, context.nfaces);
        // First call dispose() to make sure we don't leak any buffers.
        // This uses the current attributes to determine which buffers to
        // delete so don't eg. try calling removeAttribute on anything.
        geometry.dispose();
        // Now we can add in our newly calculated attribute vectors.
        geometry.addAttribute('position', new THREE.BufferAttribute(vertexarray,3));
        geometry.addAttribute('color', new THREE.BufferAttribute(colorarray,4));
        geometry.addAttribute('normal', new THREE.BufferAttribute(normalarray,3));
        if (needuvs) geometry.addAttribute('uv', new THREE.BufferAttribute(uvarray,2));
        geometry.setDrawRange(0,context.nfaces*3);
    } else if (geometry.attributes.position) {
        // Only display the parts we want
        geometry.attributes.position.needsUpdate = true;
        geometry.attributes.normal.needsUpdate = true;
        geometry.attributes.color.needsUpdate = true;
        if (needuvs) geometry.attributes.uv.needsUpdate = true;
        geometry.setDrawRange(0,context.nfaces*3);
    } else {
        return null;
    }
    if (computenormals) geometry.computeVertexNormals();
    return geometry;
}

PolyContext.prototype.runOnCanvas = function(canvas,width,height,mainwindow) {
    var Vector = Geometry.Vector;
    var PointSet = Geometry.PointSet;
    var context = this; // For benefit of subfunctions
    // http://stackoverflow.com/questions/9899807/three-js-detect-webgl-support-and-fallback-to-regular-canvas
    function webglAvailable(canvas) {
        try {
            return Boolean(window.WebGLRenderingContext && 
                           (canvas.getContext("webgl") || 
                            canvas.getContext("experimental-webgl")));
        } catch(e) { 
            return false;
        } 
    }

    // Set up the scene
    var material = new THREE.MeshPhongMaterial;
    material.side = THREE.DoubleSide;
    material.vertexColors = THREE.FaceColors;
    //material.shininess = 0;
    context.material = material;

    if (context.texturefile) {
        var onLoad = function() {
            console.log("Texture loaded: " + context.texturefile);
            context.update(PolyContext.UpdateModel);
        }
        var onError = function(error) {
            alert("Texture load failed: " + context.texturefile);
            console.log("Texture load failed: " + context.texturefile);
            console.log("Error: ", error);
            context.texture = null;
            material.map = null;
            material.needsUpdate = true;
            context.needclone = true;
            context.update(PolyContext.UpdateModel);
        }
        // Just call function directly, don't seem to need a "new"
        context.texture = new THREE.ImageUtils.loadTexture(context.texturefile,
                                                           THREE.UVMapping,
                                                           onLoad, onError);
        context.texture.wrapS = THREE.RepeatWrapping;
        context.texture.wrapT = THREE.RepeatWrapping;
        //texture.repeat.set( 2, 2 );
        //console.log(texture);
        material.map = context.texture;
        context.material.needsUpdate = true;
    }
    
    if (context.normalfile) {
        // FIXME: cut'n'paste
        var onLoad = function() {
            console.log("Normal map loaded: " + context.normalfile);
            context.update(PolyContext.UpdateModel);
        }
        var onError = function(error) {
            alert("Normal map load failed: " + context.normalfile);
            console.log("Normal map load failed: " + context.normalfile);
            console.log("Error: ", error);
            material.normalMap = null;
            material.needsUpdate = true;
            context.needclone = true;
            context.update(PolyContext.UpdateModel);
        }
        // Just call function directly, don't seem to need a "new"
        var normalmap = new THREE.ImageUtils.loadTexture(context.normalfile,
                                                         THREE.UVMapping,
                                                         onLoad, onError);
        normalmap.wrapS = THREE.RepeatWrapping;
        normalmap.wrapT = THREE.RepeatWrapping;
        //texture.repeat.set( 2, 2 );
        //console.log(texture);
        material.normalMap = normalmap;
        context.material.needsUpdate = true;
    }

    
    var scene = new THREE.Scene();
    // We just have a single main object at the moment.
    var mesh = new THREE.Mesh(new THREE.BufferGeometry(), material);
    mesh.geometry.dynamic = true;
    scene.add(mesh);
    
    var light = new THREE.DirectionalLight(0xffffff, 0.8);
    light.position.set(0,1,1);
    scene.add(light);
    var light = new THREE.DirectionalLight(0xffffff, 0.5);
    light.position.set(0,-1,0);
    scene.add(light);

    var light = new THREE.AmbientLight(0x404040); // soft white light
    scene.add(light);

    var params = { canvas: canvas, antialias: true };

    // Set this if next display loop should recompute geometry
    var needupdate = PolyContext.UpdateNone;
    this.update = function(type) {
        if (type != PolyContext.UpdateNone) {
            if (needupdate == PolyContext.UpdateNone) {
               if (context.verbose) console.log("Requesting animation frame")
                requestAnimationFrame(context.render);
            }
            if (needupdate < type) {
                needupdate = type;
            }
        }
    }
        
    this.renderer =
        (webglAvailable(canvas)) ? 
        new THREE.WebGLRenderer(params) : 
        new THREE.CanvasRenderer(params);
    //this.renderer.setPixelRatio(window.devicePixelRatio);
    this.renderer.setSize(width,height); 

    this.camera = new THREE.PerspectiveCamera(45, width/height, 0.1, 1000);
    this.camera.position.z = this.initz; 

    this.renderscene = function() { context.renderer.render(scene, context.camera); }

    if (mainwindow) {
        var controls = new THREE.OrbitControls( this.camera, canvas );
        var oldupdate = controls.update;
        controls.update = function () {
            oldupdate.call(controls);
            context.update(PolyContext.UpdateView);
            //console.log(context.camera.position);
        }
    }

    if (context.offfile || context.fnames.length > 0) {
        context.numcolorstyles = 4;
    } else {
        this.numcolorstyles = 8;
    }

    this.offcompound = function(off,geometry) {
        var basevertices = geometry.nvertices; // number of faces in the base.
        var basefaces = geometry.nfaces; // number of faces in the base.
        function flip(a) { var tmp = a[1]; a[1] = a[2]; a[2] = tmp; }
        function flipface(face) {
            flip(face.vlist);
            flip(face.normals);
            flip(face.colors);
            flip(face.uvs);
        }
        var matrix = new THREE.Matrix4;
        var normalMatrix = new THREE.Matrix3;
        function drawpolyhedron(init) {
            //console.log("drawpolyhedron " + context.compound);
            if (init) {
                // Use the OFF file vertices here as there will usually be fewer of them
                if (context.verbose) console.log("Generating pointset with " + off.vertices.length + " points")
                var pointset = [];
                for (var i = 0; i < off.vertices.length; i++) {
                    var v = off.vertices[i];
                    // Better to just append && sort
                    // Assumes no duplicates (should check really).
                    //PointSet.add([v.x,v.y,v.z],pointset,1e-4); // Quadratic!
                    pointset.push([v.x,v.y,v.z]); // TBD: PointSet should take a Vector3.
                }
                PointSet.sort(pointset);
                //console.log(pointset);
                context.pointset = pointset;
            } else {
                console.assert(context.transmat);
                var t = context.transmat;
                matrix.set(t[0][0], t[0][1], t[0][2], 0,
                           t[1][0], t[1][1], t[1][2], 0,
                           t[2][0], t[2][1], t[2][2], 0,
                                 0,       0,       0, 1);
                var det = matrix.determinant();
                normalMatrix.getNormalMatrix(matrix);
                var vertexoffset = geometry.nvertices;
                var vertices = geometry.vertices
                geometry.nvertices += basevertices;
                while (vertices.length < vertexoffset+basevertices) {
                    vertices.push(new THREE.Vector3);
                }
                for (var i = 0; i < basevertices; i++) {
                    // Copy over each vertex
                    vertices[vertexoffset+i].copy(vertices[i]);
                    vertices[vertexoffset+i].applyMatrix4(matrix);
                }
                var faceoffset = geometry.nfaces;
                geometry.nfaces += basefaces;
                var faces = geometry.faces;
                while(faces.length < faceoffset+basefaces) {
                    faces.push({ vlist: [0,0,0],
                                 uvs: [new THREE.Vector2, new THREE.Vector2, new THREE.Vector2],
                                 normals: [new THREE.Vector3, new THREE.Vector3, new THREE.Vector3],
                                 colors: [0,0,0]
                               });
                }
                for (var i = 0; i < basefaces; i++) {
                    var face1 = geometry.faces[i];
                    var face2 = geometry.faces[faceoffset+i];
                    face2.vlist[0] = face1.vlist[0]+vertexoffset;
                    face2.vlist[1] = face1.vlist[1]+vertexoffset;
                    face2.vlist[2] = face1.vlist[2]+vertexoffset;
                    face2.uvs[0].copy(face1.uvs[0]);
                    face2.uvs[1].copy(face1.uvs[1]);
                    face2.uvs[2].copy(face1.uvs[2]);
                    face2.normals[0].copy(face1.normals[0]);
                    face2.normals[1].copy(face1.normals[1]);
                    face2.normals[2].copy(face1.normals[2]);
                    face2.normals[0].applyMatrix3(normalMatrix);
                    face2.normals[1].applyMatrix3(normalMatrix);
                    face2.normals[2].applyMatrix3(normalMatrix);
                    face2.colors[0] = face1.colors[0];
                    face2.colors[1] = face1.colors[1];
                    face2.colors[2] = face1.colors[2];
                    if (det < 0) flipface(face2);
                }
            }
        }
        context.drawcompound(drawpolyhedron);
    }

    if (this.offfile) {
        // Load the offfile into the main geometry object
        offload(this.offfile, context);
    }

    // Keep track of time
    this.stopwatch = new PolyContext.Stopwatch(this.initt,
                                               this.animstep,
                                               this.initrunning);

    var basiccolors = ["red","blue","yellow","lime","yellow","aqua","silver"];

    var faceColors = basiccolors.map(function(c){ return new THREE.Color(c); });
    var white = new THREE.Color(0xffffff);
    var black = new THREE.Color(0x000000);

    function getcolor(colors, index) {
        return colors[(index + context.coloroffset)%colors.length];
    }
    this.colorface = function(face,type,index,compound) {
        if (!face.colors) face.colors = [0,0,0];
        var colorstyle = this.colorstyle;
        var color;
        if (colorstyle == 0) {
            color = getcolor(compoundColors,compound);
        } else if (colorstyle == 1) {
            color = getcolor(faceColors,type);
        } else if (colorstyle == 2) {
            color = getcolor(compoundColors,index);
        } else if (colorstyle == 6) {
            color = getcolor(compoundColors,0);
        } else if (colorstyle == 7) {
            color = white;
        }
        if (color) {
            face.colors[0] = face.colors[1] = face.colors[2] = color;
        } else {
            if (colorstyle == 3) {
                face.colors[0] = getcolor(compoundColors,0);
                face.colors[1] = getcolor(compoundColors,1);
                face.colors[2] = getcolor(compoundColors,2);
            } else if (colorstyle == 4) {
                var facecolor = getcolor(faceColors,type);;
                face.colors[0] = black;
                face.colors[1] = face.colors[2] = facecolor;
            } else if (colorstyle == 5) {
                var facecolor = getcolor(compoundColors,compound);
                face.colors[0] = black;
                face.colors[1] = face.colors[2] = facecolor;
            } else {
                console.assert(false,"No colorstyle defined");
            }
        }
    }

    this.drawpoint0 = function (x,y,z) {
        if (this.npoints == this.vertices.length) {
            this.vertices.push(new THREE.Vector3());
        }
        var u = this.vertices[this.npoints];
        u.x = x; u.y = y; u.z = z;
        return this.npoints++;
    }

    this.drawtriangle0 = function(a,b,c,uvs,type,index) {
        var needuvs = context.texturefile;
        if (this.nfaces == this.faces.length) {
            this.faces.push({ vlist: [0,0,0],
                              uvs: [new THREE.Vector2,new THREE.Vector2,new THREE.Vector2] });
            this.needclone = true; // Not sure I really need this now.
        }
        var face = this.faces[this.nfaces];
        face.vlist[0] = a;
        face.vlist[1] = b;
        face.vlist[2] = c;
        if (needuvs) {
            face.uvs[0].copy(uvs[0]);
            face.uvs[1].copy(uvs[1]);
            face.uvs[2].copy(uvs[2]);
        }
        this.colorface(face,type,index,this.compound);
        this.nfaces++;
    }

    // Set up constant aspects of polyhedra generation
    // schwarz has all the details about the particular set of
    // fundamental regions we are working with.
    var Schwarz = Geometry.Schwarz;
    this.schwarz = new Schwarz(this.angles);
    //this.schwarz.describe(false);

    var sym = Schwarz.makesymmetry(this.initsym);
    if (!sym) {
        alert("makesymmetry failed for " + this.initsym);
        return;
    } else {
        this.symmetry = [Vector.normalize(Vector.cross(sym[1],sym[2])),
                         Vector.normalize(Vector.cross(sym[2],sym[0])),
                         Vector.normalize(Vector.cross(sym[0],sym[1]))];
    }

    var rotationx = 0
    var rotationy = 0
    var clock = new THREE.Clock();
    this.render = function (tstamp) {
        var geometry = mesh.geometry
        if (needupdate == PolyContext.UpdateNone) return;
        setTimeout(function() { requestAnimationFrame( context.render ); },
                   1000 / 25 );
        var isoff = context.offfile || context.fnames.length > 0;
        try {
            // We should be able to optimize just updating colors, at least for OFF models.
            if (needupdate == PolyContext.UpdateModel || needupdate == PolyContext.UpdateColors) {
                if (isoff) {
                    context.offdata = context.offdata || {}
                    var off = context.offdata
                    var options = context.offoptions
                    for (var i = 0; i < context.fnames.length; i++) {
                        if (context.verbose) console.log("Calling ", context.fnames[i])
                        off = context[context.fnames[i]](off,context.offoptions,context.stopwatch.running)
                    }
                    if (off && off.vertices) {
                        if (!context.vertices) context.vertices = [];
                        if (!context.faces) context.faces = [];
                        //THREE.OFFLoader.display(off, context.offgeom, options);
                        THREE.OFFLoader.display2(off, context, options);
                        context.basefaces = context.nfaces;
                        if (context.docompound) {
                            context.offcompound(off,context);
                        }
                    }
                    var computenormals = false;
                    var colorstyle = context.colorstyle;
                } else {
                    // Old-style polyhedron. Should merge with OFF-style display
                    // Make sure we have the right color setting in material
                    // We should only do this when material properties actually change
                    var points = context.tours[context.tournum%context.tours.length];
                    var tripoint = context.interpolatepoint(context.stopwatch.getTime(), points);

                    if (!context.vertices) context.vertices = [];
                    if (!context.faces) context.faces = [];
                    context.drawcompound(function (init) {
                        if (init) context.setup(tripoint); // Make basic polyhedron
                        context.drawpolyhedron();
                    });
                    var computenormals = true;
                    var colorstyle = 0;
                }
                if (!context.updategeometry(geometry,computenormals,colorstyle)) return;
                // Display model normals
                // Remove old edges unless we want a persistence of vision effect.
                if (context.normalshelper) {
                    scene.remove(context.normalshelper);
                    context.normalshelper.geometry.dispose();
                }
                if (context.shownormals) {
                    context.normalshelper = new THREE.VertexNormalsHelper(mesh, 2, 0x00ff00, 1);
                    scene.add(context.normalshelper);
                }
            }
            // Render our scene
            context.renderer.render(scene, context.camera);
        } catch(err) {
            console.log("Error: ", err);
        }

        // And update for next time around
        needupdate = PolyContext.UpdateNone;

        // This is handled by THREE.js so we don't need to update model
        if (context.dorotate) {
            //var now = Date.now(); // TBD: use time here
            rotationx += 0.004
            rotationy += 0.002
            mesh.rotation.x = rotationx; 
            mesh.rotation.y = rotationy;
            needupdate = PolyContext.UpdateView;
        }
        // Both z rotation and inversion should be driven
        // from the time rather than however fast we happen
        // to render.
        if (context.dozrotate) {
            // Rotate the symmetry frame
            var thetainc = 0.0025;
            context.theta += thetainc;
            if (context.docompound) needupdate = PolyContext.UpdateModel;
        }
        var ifact = context.ifact;
        ifact += context.invertinc;
        if (ifact > 1) ifact = 1;
        if (ifact < 0) ifact = 0;
        if (ifact != context.ifact) {
            context.ifact = ifact;
            needupdate = PolyContext.UpdateModel;
        }
        // Do we need to update again?
        if (context.stopwatch.running) {
            needupdate = PolyContext.UpdateModel;
        }
    };
    this.update(PolyContext.UpdateModel);
}

// Static function, direct member of PolyContext.
PolyContext.runOnWindow = function(canvas,info,link) {
    var options = window.location.search;
    if (options.length > 0) {
        // Strip off leading '?'
        options = options.slice(1)
    }

    if (!canvas) {
        canvas = document.createElement("canvas");
        document.body.appendChild(canvas);
    }
    var width = window.innerWidth;
    var height = window.innerHeight;
    var context = new PolyContext(options);

    var showinfo = true;
    function setshowinfo() {
        if (showinfo) {
            info.style.display = 'block';
            link.style.display = 'block';
        } else {
            info.style.display = 'none';
            link.style.display = 'none';
        }
    }
    var infostring =
        '&lt;mouse drag&gt;: orbital controls, &lt;up&gt;/&lt;down&gt;: move in/out, r: rotation, f: colorstyle, ?: info display'
    
    info.innerHTML = infostring;
    setshowinfo();
    window.addEventListener("resize", function() {
        if (context.renderer) {
            var w = window.innerWidth;
            var h = window.innerHeight;
            context.renderer.setSize(w,h);
            context.camera.aspect = w/h;
            context.camera.updateProjectionMatrix();
            context.update(PolyContext.UpdateModel);
        }
    });
    window.addEventListener("keypress",
                            function (event) {
                                if (!event.ctrlKey) {
                                    // Ignore event if control key pressed.
                                    var c = String.fromCharCode(event.charCode)
                                    var handled = false;
                                    switch(c) {
                                    case '?':
                                        showinfo = !showinfo;
                                        setshowinfo();
                                        handled = true;
                                        break;
                                    default:
                                        handled = context.handleKey(c);
                                        break;
                                    }
                                    if (handled) event.preventDefault();
                                }
                            },
                            false);
    context.runOnCanvas(canvas,width,height,true);
}
