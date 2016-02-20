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
// !Sort out retroflex snub region edges bug
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
// The Fourth Dimension

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
PolyContext.UpdateModel = 2;

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
            context.fname = matches[1];
        } else if (matches = arg.match(/^off=(.+)$/)) {
            context.offfile = matches[1];
        } else if (matches = arg.match(/^vertexstyle=(.+)$/)) {
            context.offoptions.vertexstyle = matches[1];
        } else if (matches = arg.match(/^vertexwidth=([\d.]+)$/)) {
            context.offoptions.vertexwidth = Number(matches[1]);
        } else if (matches = arg.match(/^edgewidth=([\d.]+)$/)) {
            context.offoptions.edgewidth = Number(matches[1]);
        } else if (matches = arg.match(/^off.([^=]+)=([-\d.]+)$/)) {
            //console.log("Got", matches[1], matches[2])
            context.offoptions[matches[1]] = Number(matches[2]);
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
        } else if (matches = arg.match(/^spline$/)) {
            context.offoptions.dospline = true;
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
        this.updatecolors()
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
        this.updatecolors()
        break;
    case 'V':
        this.verbose = !this.verbose
        if (this.verbose) console.log("Verbose on")
        update = PolyContext.updateNone
        break
    case 'v':
        this.coloroffset-1;
        this.updatecolors()
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
        this.offoptions.C *= 1.1;
        this.needclone = true;
        break;
    case '0':
        this.offoptions.C /= 1.1;
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
// Return the index in the geometries list of points.

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
        } else if (facefact) {
            // Give each facetype a slightly different offset
            x += facefact*offset[0];
            y += facefact*offset[1];
            z += facefact*offset[2];
        }
    }
    if (this.transmat) {
        p = [x,y,z];
        x = Vector.dot(this.transmat[0],p);
        y = Vector.dot(this.transmat[1],p);
        z = Vector.dot(this.transmat[2],p);
    }
    return this.drawpoint0(x,y,z);
}

// TBD: encapsulate the various drawing options
// For n > 0, draw a Sierpinski triangle (or some
// variation thereof).
PolyContext.prototype.drawtriangle = function(p,q,r,tcoords,offset,type,i,n) {
    var Vector = Geometry.Vector;
    if (n == 0) {
        var facefact = (this.compound+type)*0.0001 // Reduce stitching
        var index0 = this.drawpoint(p,offset,facefact);
        var index1 = this.drawpoint(q,offset,facefact);
        var index2 = this.drawpoint(r,offset,facefact);
        this.drawtriangle0(index0,index1,index2,tcoords,type,i);
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
        if (k != 0) this.drawtriangle(p,p1,r1,tcoords,offset,type,i,n-1);
        if (k != 1) this.drawtriangle(p1,q,q1,tcoords,offset,type,i,n-1);
        if (k != 2) this.drawtriangle(r1,q1,r,tcoords,offset,type,i,n-1);
        if (k != 3) this.drawtriangle(p1,q1,r1,tcoords,offset,type,i,n-1);
    }
}

PolyContext.prototype.drawface = function(centre,plist,facetype,index,tridepth,parity) {
    var Vector = Geometry.Vector;
    if (this.drawface0) {
        this.drawface0(plist); // Call alternative function, if defined
    } else {
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

        if (this.dohemi) {
            // Do this first as might modify plist
            // TBD: Proper texture coordinates!
            var uvs = [new THREE.Vector2(0,0),
                       new THREE.Vector2(1,0),
                       new THREE.Vector2(0,1)];
            for (var i = 0; i < plist.length; i++) {
                var p1 = plist[i];
                var p2 = plist[(i+1)%plist.length];
                if (Vector.taxi(p1,p2) > 1e-4) {
                    this.drawtriangle([0,0,0],p1,p2,uvs,
                                      centre,4,index,tridepth);
                }
            }
        }

        // Compute uvs for face points
        var uvs = []; // Try not to clash with hemi uvs
        for (var i = 0; i < plist.length; i++) {
            var u = 0.5 + Vector.dot(Vector.sub(plist[i],centre),uaxis);
            var v = 0.5 + Vector.dot(Vector.sub(plist[i],centre),vaxis);
            uvs.push(new THREE.Vector2(u,v));
        }
        var uvcentre = new THREE.Vector2(0.5,0.5);

        if (plist.length == 3) {
            // Should do this when some edges are degenerate too.
            this.drawtriangle(plist[0],plist[1],plist[2],
                              [uvs[0],uvs[1],uvs[2]],
                              centre,facetype,index,tridepth);
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
                var newuvs = []
                var mod = function(i,n) {
                    i %= n
                    if (i < 0) i += n
                    return i;
                }
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
                        // q is new point to draw triangle to edge from
                        var u = 0.5 + Vector.dot(Vector.sub(q,centre),uaxis);
                        var v = 0.5 + Vector.dot(Vector.sub(q,centre),vaxis);
                        var newuv = new THREE.Vector2(u,v);
                        newuvs.push(newuv)
                        newplist.push(q)
                        this.drawtriangle(q,p1,p2,
                                          [newuv,uvs[i],uvs[(i+1)%plist.length]],
                                          centre,facetype,index,tridepth);
                    }
                }
                if (newplist.length > 0) {
                    plist = newplist;
                    uvs = newuvs;
                } else {
                    // No retroflex needed for other faces.
                    this.doretroflex = false;
                }
            }
            console.assert(uvs.length == plist.length)
            for (var i = 0; i < plist.length; i++) {
                var p1 = plist[i];
                var p2 = plist[(i+1)%plist.length];
                if (Vector.taxi(p1,p2) > 1e-4) { // Don't draw if degenerate
                    this.drawtriangle(centre,p1,p2,
                                      [uvcentre,uvs[i],uvs[(i+1)%plist.length]],
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
            this.drawface(centre,plist,facetype,i,this.tridepth,region[3]);
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
    // TBD: do this properly!
    var uvs = [new THREE.Vector2(0,0),
               new THREE.Vector2(1,0),
               new THREE.Vector2(0,1)];
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
    //this.renderscene();
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
    loader.load( file,
                 function (off) {
                     context.offdata = off
                     context.needclone = true
                     context.update(PolyContext.UpdateModel);
                 },
                 onProgress, onError);
}

PolyContext.prototype.runOnCanvas = function(canvas,width,height) {
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
    //material.color = 0xffffff;
    material.side = THREE.DoubleSide;
    //material.side = THREE.FrontSide;
    material.vertexColors = THREE.FaceColors;
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
    var mesh = new THREE.Mesh(new THREE.Geometry(), material);
    mesh.geometry.dynamic = true;
    context.offgeom = mesh.geometry
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

    var controls = new THREE.OrbitControls( this.camera, canvas );
    var oldupdate = controls.update;
    controls.update = function () {
        oldupdate.call(controls);
        context.update(PolyContext.UpdateView);
        //console.log(context.camera.position);
    }

    this.offcompound = function(off) {
        var offgeom = context.offgeom;
        var offcolors = context.offcolors;
        var geometry = new THREE.Geometry;
        function drawpolyhedron(init) {
            //console.log("drawpolyhedron " + context.compound);
            if (context.docompound && !context.pointset) {
                // Use the OFF file vertices here as there will usually be fewer of them
                if (context.verbose) console.log("Generating pointset with " + off.vertices.length + " points")
                var pointset = [];
                for (var i = 0; i < off.vertices.length; i++) {
                    var v = off.vertices[i];
                    // Better to just append && sort
                    // Assumes no duplicates (should check really).
                    //PointSet.add([v.x,v.y,v.z],pointset,1e-4); // Quadratic!
                    pointset.push([v.x,v.y,v.z]);
                }
                PointSet.sort(pointset);
                context.pointset = pointset;
            }
            var m, det;
            if (context.transmat) {
                console.assert(!init, "init set when matrix supplied");
                var t = context.transmat;
                m = new THREE.Matrix4;
                m.set(t[0][0], t[0][1], t[0][2], 0,
                      t[1][0], t[1][1], t[1][2], 0,
                      t[2][0], t[2][1], t[2][2], 0,
                            0,       0,       0, 1);
                det = m.determinant();
            }
            for (var i = 0; i < offgeom.faces.length; i++) {
                if (context.colorstyle === 0 && offcolors) {
                    offgeom.faces[i].color.copy(offcolors[i]);
                } else if (context.colorstyle === 1) {
                    offgeom.faces[i].color.copy(white);
                } else {
                    var compound = context.compound
                    console.assert(compound != undefined)
                    offgeom.faces[i].color.copy(compoundColors[(context.coloroffset + compound)%compoundColors.length]);
                }
            }
            // This sure is a lousy hack
            // For reflections, need to change order of triangles so they stay front facing,
            // else it seems something flips the normals around.
            // It would be better if Geometry.merge did this itself.
            function flipface(i) {
                var f = offgeom.faces[i];
                var tmp = f.b; f.b = f.c; f.c = tmp;
                if (f.vertexNormals.length > 0) {
                    tmp = f.vertexNormals[1]
                    f.vertexNormals[1] = f.vertexNormals[2]
                    f.vertexNormals[2] = tmp
                }
                if (f.vertexColors.length > 0) {
                    tmp = f.vertexColors[1]
                    f.vertexColors[1] = f.vertexColors[2]
                    f.vertexColors[2] = tmp
                }
                var uvs = offgeom.faceVertexUvs[0][i];
                tmp = uvs[1]; uvs[1] = uvs[2]; uvs[2] = tmp;
                var normals = offgeom.faceVertexUvs[0][i];
                tmp = uvs[1]; uvs[1] = uvs[2]; uvs[2] = tmp;
            }
                
            if (det < 0) {
                console.assert(!init);
                for (var i = 0; i < offgeom.faces.length; i++) flipface(i);
            }
            // Should be able to merge geometries with eg. a face
            // coloring function etc.
            geometry.merge(offgeom,m);
            if (det < 0) {
                console.assert(!init);
                for (var i = 0; i < offgeom.faces.length; i++) flipface(i)
            }
        }
        //var start = Date.now();
        context.drawcompound(drawpolyhedron);
        //console.log("t=" + (Date.now() - start));
        scene.remove(mesh);
        mesh.geometry.dispose();
        mesh.geometry = geometry;
        scene.add(mesh);
    }
    this.updatecolors = function() {
        var geometry = mesh.geometry
        if (context.offcolors) {
            var ncolors = context.offcolors.length;
            // Bzzt! Repetition!
            for (var i = 0; i < geometry.faces.length; i++) {
                if (context.colorstyle === 0) {
                    geometry.faces[i].color.copy(context.offcolors[i%ncolors]);
                } else if (context.colorstyle === 1) {
                    geometry.faces[i].color.copy(white);
                } else {
                    var compound = 0
                    geometry.faces[i].color.copy(compoundColors[(context.coloroffset + compound)%compoundColors.length]);
                }
            }
        }
        geometry.colorsNeedUpdate = true;
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

    var compoundColors = [
        new THREE.Color(0xff0000),new THREE.Color(0x00ff00),new THREE.Color(0x0000ff),
        new THREE.Color(0xffff00),new THREE.Color(0x00ffff),new THREE.Color(0xff00ff),
        new THREE.Color(0x808080),new THREE.Color(0xffffff),
    ];
    var faceColors = basiccolors.map(function(c){ return new THREE.Color(c); });
    var white = new THREE.Color(0xffffff);
    var black = new THREE.Color(0x000000);

    this.numcolorstyles = 8;
    function getcolor(colors, index) {
        return colors[(index + context.coloroffset)%colors.length];
    }
    this.colorface = function(face,type,index,compound) {
        var colorstyle = this.colorstyle;
        if (colorstyle == 0) {
            face.color = getcolor(compoundColors,compound);
        } else if (colorstyle == 1) {
            face.color = getcolor(faceColors,type);
        } else if (colorstyle == 2) {
            face.color = getcolor(compoundColors,index);
        } else if (colorstyle == 3) {
            face.vertexColors = [getcolor(compoundColors,0),
                                 getcolor(compoundColors,1),
                                 getcolor(compoundColors,2)];
        } else if (colorstyle == 4) {
            var facecolor = getcolor(faceColors,type);;
            face.vertexColors = [black,facecolor,facecolor];
        } else if (colorstyle == 5) {
            var facecolor = getcolor(compoundColors,compound);
            face.vertexColors = [black,facecolor,facecolor];
        } else if (colorstyle == 6) {
            face.color = getcolor(compoundColors,0);
        } else if (colorstyle == 7) {
            face.color = white;
        }
    }
    this.colortype = function() {
        var colorstyle = this.colorstyle;
        if (colorstyle == 3 || colorstyle == 4 || colorstyle == 5) {
            return THREE.VertexColors;
        } else {
            return THREE.FaceColors;
        }
    }

    // TBD: Defining a window event handler here seems wrong
    window.addEventListener("keypress",
                            function (event) {
                                if (!event.ctrlKey) {
                                    // Ignore event if control key pressed.
                                    var handled = context.handleKey(String.fromCharCode(event.charCode));
                                    if (handled) event.preventDefault();
                                }
                            },
                            false);

    this.drawpoint0 = function (x,y,z) {
        var geometry = mesh.geometry
        // If we need to add new points to the geometry
        // make sure that we redraw properly.
        if (this.npoints == geometry.vertices.length) {
            geometry.vertices.push(new THREE.Vector3());
            this.needclone = true;
        }
        var u = geometry.vertices[this.npoints];
        u.x = x; u.y = y; u.z = z;
        return this.npoints++;
    }

    this.drawtriangle0 = function(a,b,c,tcoords,type,index) {
        console.assert(tcoords);
        var geometry = mesh.geometry
        if (this.nfaces == geometry.faces.length) {
            geometry.faces.push(new THREE.Face3());
            console.assert(geometry.faceVertexUvs[0].length === this.nfaces);
            geometry.faceVertexUvs[0].push([0,0,0]);
            this.needclone = true;
        }
        var face = geometry.faces[this.nfaces];
        face.a = a; face.b = b; face.c = c;
        if (tcoords) {
            var uv = geometry.faceVertexUvs[0][this.nfaces];
            uv[0] = tcoords[0];
            uv[1] = tcoords[1];
            uv[2] = tcoords[2];
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
    this.render = function (tstamp) {
        var geometry = mesh.geometry
        if (needupdate == PolyContext.UpdateNone) return;
        setTimeout(function() { requestAnimationFrame( context.render ); },
                   1000 / 25 );
        if (needupdate == PolyContext.UpdateModel) {
            if (context.offfile || context.fname) {
                var off = context.offdata
                var options = context.offoptions
                // Maybe this should take a sequence of functions
                if (context.fname) {
                    if (context.verbose) console.log("Calling ", context.fname)
                    off = context[context.fname](off,context.offoptions,context.stopwatch.running)
                    //context.needclone = true;
                }
                if (off) {
                    var nvertices1 = context.offgeom.vertices.length;
                    var nfaces1 = context.offgeom.faces.length;
                    THREE.OFFLoader.display(off, context.offgeom, options);
                    if (nvertices1 != context.offgeom.vertices.length ||
                        nfaces1 != context.offgeom.faces.length) {
                        context.needclone = true;
                    }
                    var geometry = context.offgeom
                    context.numcolorstyles = 3; // Yuk.
                    context.offcolors = [];
                    for (var i = 0; i < geometry.faces.length; i++) {
                        var c = geometry.faces[i].color.clone()
                        context.offcolors.push(c);
                    }
                    if (context.docompound) {
                        context.offcompound(off);
                        if (context.verbose) console.log("Compound of", context.compound);
                    } else if (context.needclone) {
                        context.needclone = false
                        if (context.verbose) console.log("Cloning")
                        scene.remove(mesh);
                        mesh.geometry.dispose()
                        mesh.geometry = context.offgeom.clone()
                        context.offgeom = mesh.geometry
                        context.updatecolors();
                        scene.add(mesh);
                    } else {
                        context.updatecolors();
                        geometry.verticesNeedUpdate = true;
                        geometry.elementsNeedUpdate = true;
                        geometry.normalsNeedUpdate = true;
                        geometry.colorsNeedUpdate = true;
                        geometry.uvsNeedUpdate = true;
                    }
                }
            } else {
                // Old-style polyhedron. Should merge with OFF-style display
                // Make sure we have the right color setting in material
                // We should only do this when material properties actually change
                context.material.vertexColors = context.colortype();
                var points = context.tours[context.tournum%context.tours.length];
                var tripoint = context.interpolatepoint(context.stopwatch.getTime(), points);

                context.drawcompound(function (init) {
                    if (init) context.setup(tripoint); // Make basic polyhedron
                    context.drawpolyhedron();
                });
                // Prepare for rendering...
                console.assert(geometry.faces.length === geometry.faceVertexUvs[0].length);

                // Looks like we need to clone if the number of vertices
                // or faces has changed.
                if (geometry.vertices.length != context.npoints ||
                    geometry.faces.length != context.nfaces) {
                    context.needclone = true;
                    geometry.vertices.length = context.npoints;
                    geometry.faces.length = context.nfaces;
                    geometry.faceVertexUvs[0].length = context.nfaces;
                }
                geometry.computeFaceNormals();
                geometry.verticesNeedUpdate = true;
                geometry.elementsNeedUpdate = true;
                geometry.normalsNeedUpdate = true;
                geometry.colorsNeedUpdate = true;
                geometry.uvsNeedUpdate = true;
                if (context.needclone) {
                    context.needclone = false;
                    //console.log("Cloning geometry: " + context.npoints + " " + context.nfaces);
                    // This is rather horrible, but I haven't found a more economical
                    // way of doing this with the three.js version I'm using (and it's all
                    // changed (for the better I'm sure) in more recent versions).
                    // We don't do this every frame at least.
                    var newgeometry = geometry.clone()
                    scene.remove(mesh);
                    geometry.dispose();
                    geometry = newgeometry;
                    // Reusing the mesh object seems OK
                    mesh.geometry = geometry; 
                    scene.add(mesh);
                }
            }
            // Display model normals
            // Remove old edges unless we want a persistence of vision effect.
            if (context.normalshelper) {
                scene.remove(context.normalshelper);
                context.normalshelper.geometry.dispose();
            }
            if (context.shownormals) {
                context.normalshelper = new THREE.FaceNormalsHelper(mesh, 2, 0x00ff00, 1);
                scene.add(context.normalshelper);
            }
        }
        // Render our scene
        context.renderer.render(scene, context.camera);

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
PolyContext.runOnWindow = function(canvas) {
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
    context.runOnCanvas(canvas,width,height);
}
