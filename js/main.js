"use strict";

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
// Modularize
// Table driven options
// Better separation of generation and display
// Separate logical face structure from display structure
// Hollow faces
// Neo filling
// zrotation should be function of time
// Don't draw degenerate sectors of faces
// Invertible snub duals 
// Snub stellations
// Dual of stellations
// Stellation of duals
// More generic dual generation
// Stellation parameters
// Use proper options object rather than globals
// Make URL for current scene
// More compact URL scheme
// Better colors
// Symmetric colors
// Use homogenous coordinates
// The Fourth Dimension

// Options, modifiable during runtime
var dostellate = false;
var dosnubify = false;
var drawtype = 0;
var invertinc = -0.01;
var ifact = 0;
var donormalize = false;
var docompound = false;
var dorotate = false;
var dozrot = false;
var dorotonly = false;
var colorstyle = 0;
var hideface = [false,false,false,false]
var zrotation = 0;
var tridepth = 0;
var regionstyle = 0;
var explode = 0;
var dohemi = false;

// Constants
var numcolorstyles = 6;
var vertexcolorstyle = 3;

// Take a walk around the Schwarz triangle
var tours = [
    [],
    [ [1,0,0],[1,1,0],[0,1,0],[0,1,1],[0,0,1],[1,0,1],[1,0,0],
      [1,1,1],[1,1,0],[0,1,0],[1,1,1],[0,1,1],[0,0,1],[1,1,1],[1,0,1] ],
    [ [2,-1,0], [1,0,0], [0,1,0], [-1,2,0],
      [0,2,-1], [0,1,0], [0,0,1], [0,-1,2],
      [-1,0,2], [0,0,1], [1,0,0], [2,0,-1] ],
    [ [2,-1,0],[0,-1,2],[-1,0,2], [-1,2,0],[0,2,-1],[2,0,-1] ]
];
var tournum = 0;

function Stopwatch(init,step,running) {
    // When running, reported time is difference between Date.now()
    // and stored time. When stopped, reported time is just the
    // stored time.
    this.toggle = function() { 
        this.running = !this.running;
        this.time = Date.now()-this.time;
    }
    this.setTime = function(t) { 
        if (this.running) this.time = Date.now()-t*this.scale;
        else this.time = t*this.scale;
    }
    this.incTime = function(t) { 
        if (this.running) {
            this.time -= t*this.scale;
        } else {
            // Apply increment and step to next position
            this.time += t*this.scale;
            this.time = this.scale*(Math.floor(this.time/this.scale));
        }
    }
    this.getTime = function() { 
        if (this.running) return (Date.now()-this.time)/this.scale;
        else return this.time/this.scale;
    }
    this.running = running;
    this.scale = step*1000;  // Return time in seconds
    this.setTime(init || 0); // Different semantics in running and non-running states
}

function maketri(t,points) {
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

// Actually generate a THREE point
// We might transform it first.
// Return the point index.
// This is sort of time-critical so inline a lot of stuff
// Maybe Closure can do this for us.
function drawpoint(context,p,offset)
{
    var w = p[3] || 1;
    // For fun, do the explode before the inversion
    if (donormalize) w = Vector.length(p);
    if (ifact != 0) {
        //var q = invert(p,context.midsphere);
        var k = Vector.dot(p,p)/context.midsphere;
        if (ifact == 1) w *= k;
        else w *= (1-ifact)+ifact*k;
    }
    var x = p[0]/w, y = p[1]/w, z = p[2]/w;
    if (offset && explode != 0) {
        x += explode*offset[0];
        y += explode*offset[1];
        z += explode*offset[2];
    }
    if (context.transmat) {
        p = [x,y,z];
        x = Vector.dot(context.transmat[0],p);
        y = Vector.dot(context.transmat[1],p);
        z = Vector.dot(context.transmat[2],p);
    }
    if (context.npoints == context.geometry.vertices.length) {
        context.geometry.vertices.push(new THREE.Vector3());
        context.needclone = true;
    }
    var u = context.geometry.vertices[context.npoints];
    u.x = x; u.y = y; u.z = z;
    return context.npoints++;
}

function drawtriangle0(context,a,b,c,type,index)
{
    if (context.nfaces == context.geometry.faces.length) {
        context.geometry.faces.push(new THREE.Face3());
        context.needclone = true;
    }
    var face = context.geometry.faces[context.nfaces];
    face.a = a; face.b = b; face.c = c;
    context.colorface(face,type,index,context.compound);
    context.nfaces++;
}

// TBD: encapsulate the various drawing options        
function drawtriangle(context,p,q,r,offset,type,i,n)
{
    var geometry = context.geometry
    if (n == 0) {
        var index0 = drawpoint(context,p,offset);
        var index1 = drawpoint(context,q,offset);
        var index2 = drawpoint(context,r,offset);
        drawtriangle0(context,index0,index1,index2,type,i);
    } else {
        var p1 = Vector.mid(p,q);
        var q1 = Vector.mid(q,r);
        var r1 = Vector.mid(r,p);
        // TBD: Different triangulation styles here
        // k = 0,1,2,3,n%4,random etc.
        //var k = n%4;
        var k = 3;
        if (k != 0) drawtriangle(context,p,p1,r1,offset,type,i,n-1);
        if (k != 1) drawtriangle(context,p1,q,q1,offset,type,i,n-1);
        if (k != 2) drawtriangle(context,r1,q1,r,offset,type,i,n-1);
        if (k != 3) drawtriangle(context,p1,q1,r1,offset,type,i,n-1);
    }
}

var origin = [0,0,0]
function drawface(context,centre,plist,facetype,i,tridepth)
{
    if (plist.length == 3) {
        drawtriangle(context,
                     plist[0],plist[1],plist[2],
                     centre,facetype,i,tridepth);
    } else {
        for (var j = 0; j < plist.length; j++) {
            drawtriangle(context,
                         centre,plist[j],plist[(j+1)%plist.length],
                         centre,facetype,i,tridepth);
        }
    }
    if (dohemi) {
        for (var j = 0; j < plist.length; j++) {
            drawtriangle(context,
                         origin,plist[j],plist[(j+1)%plist.length],
                         centre,4,i,tridepth);
        }
    }
}

// Draw the Schwarz triangles (fundamental regions) as triangles.
function drawregions(context)
{
    var schwarz = context.schwarz;
    var facedata = context.facedata;
    var geometry = context.geometry
    var points = schwarz.points;
    var regions = schwarz.regions;
    var adjacent = schwarz.adjacent;
    var regionpoints = facedata.regionpoints;
    for (var i = 0; i < regions.length; i++) {
        var region = regions[i];
        var plist = [];
        if (region[3] == 0 || !dosnubify) {
            var snubcentre = facedata.snubcentre;
            for (var j = 0; j < 3; j++) {
                plist.push(Vector.mul(points[region[j]], facedata.facedistances[j]));
                if (dosnubify && snubcentre) {
                    // This isn't the 'regionpoint' but the 'snubpoint' ie. the
                    // centre of the snub triangle for the face.
                    var adj = adjacent[i][(j+2)%3];
                    plist.push(schwarz.applybary(snubcentre,adj));
                } else if (regionstyle == 1) {
                    // Add in edge points for "correct" inversion
                    var edgecentre = schwarz.applybary(facedata.edgecentres[j],i);
                    plist.push(edgecentre);
                }
            }
            for (var j = 0; j < plist.length; j++) {
                plist[j] = invert(plist[j],context.midsphere);
            }
            var facetype = 4+region[3];
            var centre = invert(regionpoints[i],context.midsphere);
            drawface(context,centre,plist,facetype,i,tridepth)
        }
    }
}

function drawfaces(context)
{
    var schwarz = context.schwarz;
    var facedata = context.facedata;
    var geometry = context.geometry
    var points = schwarz.points;
    var regions = schwarz.regions;
    var regionpoints = facedata.regionpoints;
    // Now go through the facepoint lists for each type of face
    for (var type = 0; type < 3; type++) {
        if (!facedata.faces[type] || hideface[type]) continue;
        if (dostellate) {
            var facepoints = schwarz.faces[type];
            for (var i = 0; i < facepoints.length; i++) {
                var plist = facepoints[i];
                if (plist) {
                    // Add a centre vertex, this is a corner of a Schwarz triangle.
                    // facedistances[type] tells us how far it is from the origin (the face
                    // centres, ie. the region points are at distance 1.
                    var centre = Vector.mul(points[i],facedata.facedistances[type])                    
                    var facelines = facedata.facelines[type];
                    for (var j = 0; j < facelines.length; j++) {
                        // Triangle from centre to two vertices on plist
                        var p1 = schwarz.applybary(facelines[j][0],plist[0]);
                        var p2 = schwarz.applybary(facelines[j][1],plist[0]);
                        drawtriangle(context,centre,p1,p2,centre,type,i,tridepth);
                    }
                }
            }
        } else {
            var facepoints = schwarz.faces[type];
            for (var i = 0; i < facepoints.length; i++) {
                var fpoints = facepoints[i];
                if (fpoints && (!dosnubify || fpoints.length > 4)) {
                    // Add a centre vertex, this is a corner of a Schwarz triangle.
                    // facedistances[type] tells us how far it is from the origin (the face
                    // centres, ie. the region points are at distance 1.
                    var centre = Vector.mul(points[i],facedata.facedistances[type]);
                    // For a snub polyhedron, we only include even (+ve) regions.
                    var inc = dosnubify ? 2 : 1;
                    var plist = [];
                    for (var j = 0; j < fpoints.length; j+=inc) {
                        plist.push(regionpoints[fpoints[j]]);
                    }
                    drawface(context,centre,plist,type,i,tridepth);
                }
            }
        }
    }
    if (dosnubify && !dostellate && !hideface[3] && facedata.snubcentre) {
        // Now draw a snub triangle for -ve regions
        for (var i = 0; i < regions.length; i++) {
            var t = regions[i];
            if (t[3] == 1) {
                var centre = schwarz.applybary(facedata.snubcentre, i);
                var plist = [regionpoints[schwarz.adjacent[i][0]],
                             regionpoints[schwarz.adjacent[i][1]],
                             regionpoints[schwarz.adjacent[i][2]]];
                drawface(context,centre,plist,3,i,tridepth);
            }
        }
    }
}

// Construct a polyhedron in the geometry object,
// based on the fundamental regions info in schwarz, using
// trilinear coords tri, and with the colorface coloring function.
function initcontext(context,schwarz,tri)
{
    var regions = schwarz.regions;
    var points = schwarz.points;
    // Convert trilinear coords to barycentric.
    if (!dosnubify) {
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
    var facedata = makefacedata(schwarz,bary);
    if (docompound) {
        var regionpoints = facedata.regionpoints;
        // Set up the initial pointset (for compound generation).
        var pointset = [];
        for (var i = 0; i < regionpoints.length; i++) {
            if (!dosnubify || regions[i][3] == 0) {
                setadd(regionpoints[i],pointset,1e-6);
            }
        }
        setsort(pointset);
        context.pointset = pointset;
    }
    context.tri = tri
    context.bary = bary
    context.schwarz = schwarz;
    context.facedata = facedata;
    if (dosnubify) {
        context.midsphere = facedata.snubsphere;
    } else {
        context.midsphere = facedata.midsphere;
    }
    if (dostellate) {
        stellate(schwarz,facedata, hideface);
    }
}

function drawpolyhedron(context)
{
    // Do the actual drawing
    if (drawtype <= 1) drawfaces(context);
    if (drawtype >= 1) drawregions(context);
}

function makeangles(a,b,c,d,e,f)
{
    b = b || 1; d = d || 1; f = f || 1;
    //console.log(a,b,c,d,e,f);
    a = Number(a); b = Number(b);
    c = Number(c); d = Number(d);
    e = Number(e); f = Number(f);
    // Check b/a + d/c + f/e > 1
    // ie. bce + dea + fac > ace
    if (a < b || c < d || e < f ||
        b*c*e + d*e*a + f*a*c <= a*c*e) {
        alert("Not a Schwarz Triangle");
    }
    //console.log(b*c*e,d*e*a,f*a*c,a*c*e);
    //console.log(a,b,c,d,e,f);
    return [a/b,c/d,e/f];
}

// http://stackoverflow.com/questions/9899807/three-js-detect-webgl-support-and-fallback-to-regular-canvas
function webglAvailable() {
    try {
        var canvas = document.createElement("canvas");
        return Boolean(window.WebGLRenderingContext && 
                       (canvas.getContext("webgl") || 
                        canvas.getContext("experimental-webgl")));
    } catch(e) { 
        return false;
    } 
}

(function () {
    // Our state object
    function Context(geometry,colorface)
    {
        this.geometry = geometry;
        this.colorface = colorface;
        this.npoints = 0;
        this.nfaces = 0;
        this.needclone = false;
        this.compound = 0;
        this.pointset = [];
        this.transmat = null;
    }
    var context;

    // Initial values
    var initrunning = false;
    var initz = 3;
    var initsym = [2,3,5];
    var angles  = [2,3,5];
    var initt = 0;
    var tri = null;
    var qparams = window.location.search;
    if (qparams.length > 0) {
        // Strip off leading '?'
        processargs(qparams.slice(1));
    }
    if (tours[0].length == 0) tours[0].push([1,1,1]);

    var scene = new THREE.Scene(); 
    var width = window.innerWidth;
    var height = window.innerHeight;
    var camera = new THREE.PerspectiveCamera(45, width/height, 0.1, 1000); 
    var renderer = (webglAvailable()) ? new THREE.WebGLRenderer() : new THREE.CanvasRenderer();
    renderer.setSize(width,height); 
    //console.log("setSize: " + width + " " +  height); 
    document.body.appendChild(renderer.domElement); 

    // Set up the scene
    var material = new THREE.MeshPhongMaterial({color: 0xffffff}); 
    material.side = THREE.DoubleSide;
    var geometry = new THREE.Geometry();
    geometry.dynamic = true;
    var mesh = new THREE.Mesh(geometry, material);
    scene.add(mesh);

    var light = new THREE.DirectionalLight(0xffffff, 0.8);
    light.position.set(0,1,1);
    scene.add(light);
    var light = new THREE.DirectionalLight(0xffffff, 0.5);
    light.position.set(0,-1,0);
    scene.add(light);
    camera.position.z = initz; 

    var light = new THREE.AmbientLight(0x404040); // soft white light
    scene.add(light);

    var controls = new THREE.OrbitControls( camera );
    //var controls = new THREE.TrackballControls( camera );

    // schwarz has all the details about the particular set of
    // fundamental regions we are working with.
    var schwarz = new Schwarz(angles);
    schwarz.describe(false);

    // Miscellaneous state setup
    var animstep = 10; // Seconds per animation step
    var sym = makesymmetry(initsym);
    //var ip0 = sym[0]; var iq0 = sym[1]; var ir0 = sym[2];
    var ip0 = Vector.normalize(Vector.cross(sym[1],sym[2]));
    var iq0 = Vector.normalize(Vector.cross(sym[2],sym[0]));
    var ir0 = Vector.normalize(Vector.cross(sym[0],sym[1]));

    // Keep track of our progression
    var stopwatch = new Stopwatch(initt,animstep,initrunning);

    // Set this if next display loop should recompute geometry
    var needupdate = true;

    var basiccolors = ["red","blue","yellow","lime","yellow","aqua","silver"];

    var compoundColors = [
        new THREE.Color(0xff0000),new THREE.Color(0x00ff00),new THREE.Color(0x0000ff),
        new THREE.Color(0xffff00),new THREE.Color(0x00ffff),new THREE.Color(0xff00ff),
        new THREE.Color(0x808080),new THREE.Color(0xffffff),
    ];
    var faceColors = basiccolors.map(function(c){ return new THREE.Color(c); });
    var white = new THREE.Color(0xffffff);
    var black = new THREE.Color(0x000000);

    // We use a group word to describe a particular transformation
    var gens = ["P","Q","R"];
    var rgens =  ["PQ","QR","RP"];

    function rotatecolors() {
        function rot(a) {
            var t = a[0];
            for (var i = 0; i < a.length-1; i++) {
                a[i] = a[i+1];
            }
            a[a.length-1] = t;
        }
        rot(compoundColors);
        rot(faceColors);
    }
    var colorface = function(face,type,index,compound) {
        //face.color = null;
        //face.vertexColors = null;
        if (colorstyle == 0) {
            face.color = compoundColors[compound%compoundColors.length];
        } else if (colorstyle == 1) {
            face.color = faceColors[type%faceColors.length];
        } else if (colorstyle == 2) {
            face.color = compoundColors[index%compoundColors.length];
        } else if (colorstyle == 3) {
            face.vertexColors = compoundColors.slice(0,3);
        } else if (colorstyle == 4) {
            var facecolor = faceColors[type%faceColors.length];
            face.vertexColors = [black,facecolor,facecolor];
        } else if (colorstyle == 5) {
            var facecolor = compoundColors[compound%compoundColors.length];
            face.vertexColors = [black,facecolor,facecolor];
        }
    }

    function makeinfostring()
    {
        var s = [];
        s.push("tri = [" + context.tri + "]");
        s.push("bary = [" + context.bary + "]");
        if (docompound) {
            s.push("z = " + zrotation);
        }
        return s.join("; ");
    }
    function processargs(params)
    {
        // TBD: A more compact parameter format would be good.
        var args = params.split('&');
        var matches;
        args.forEach(function (arg) {
            //console.log("Doing parameter '" + arg + "'");
            if (matches = arg.match(/^args=([\d]+)(?:\/([\d]+))?:([\d]+)(?:\/([\d]+))?:([\d]+)(?:\/([\d]+))?$/)) {
                angles = makeangles(matches[1],matches[2],matches[3],matches[4],matches[5],matches[6]);
            } else if (matches = arg.match(/^sym=([\d]+):([\d]+):([\d]+)$/)) {
                initsym = makeangles(matches[1],1,matches[2],1,matches[3],1);
            } else if (matches = arg.match(/^tri=([\d.]+):([\d.]+):([\d.]+)$/)) {
                tours[0].push([Number(matches[1]),Number(matches[2]),Number(matches[3])]);
            } else if (matches = arg.match(/^zrot=([\d]+)$/)) {
                zrotation = Math.PI/Number(matches[1]);
            } else if (matches = arg.match(/^zrot=([\d]+.[\d]+)$/)) {
                zrotation = Number(matches[1]);
            } else if (matches = arg.match(/^colorstyle=([\d]+)$/)) {
                colorstyle = Number(matches[1]);
            } else if (matches = arg.match(/^hide=([\d]+)$/)) {
                hideface[Number(matches[1])-1] = true;
            } else if (matches = arg.match(/^t=([\d.]+)$/)) {
                initt = Number(matches[1]);
            } else if (matches = arg.match(/^z=([\d.]+)$/)) {
                initz = Number(matches[1]);
            } else if (matches = arg.match(/^tridepth=([\d]+)$/)) {
                tridepth = Number(matches[1]);
            } else if (matches = arg.match(/^tour=([\d]+)$/)) {
                tournum = Number(matches[1]);
            } else if (matches = arg.match(/^compound$/)) {
                docompound = true;
            } else if (matches = arg.match(/^snub$/)) {
                dosnubify = true;
            } else if (matches = arg.match(/^stellate$/)) {
                dostellate = true;
            } else if (matches = arg.match(/^dual$/)) {
                drawtype = 2;
            } else if (matches = arg.match(/^invert$/)) {
                invertinc = -invertinc;
                ifact = 1;
            } else if (matches = arg.match(/^dozrot$/)) {
                dozrot = true;
            } else if (matches = arg.match(/^dorotonly$/)) {
                dorotonly = true;
            } else if (matches = arg.match(/^dorotate$/)) {
                dorotate = true;
            } else if (matches = arg.match(/^hemi$/)) {
                dohemi = true;
            } else if (matches = arg.match(/^run$/)) {
                initrunning = true;
            } else {
                console.log("Ignoring parameter '" + arg + "'");
            }
        });
    }

    // Keyboard handlers
    function onKeyPress(event) {
        var handled = true;
        switch (String.fromCharCode(event.charCode)) {
        case ' ':
            stopwatch.toggle();
            break;
        case 'z':
            dozrot = !dozrot;
            //if (!dozrot) console.log(zrotation);
            break;
        case 'a':
            zrotation = 0;
            break;
        case '[':
            stopwatch.incTime(-1);
            break;
        case ']':
            stopwatch.incTime(1);
            break;
        case 'p':
            stopwatch.setTime(0);
            break;
        case 'r':
            dorotate = !dorotate;
            break;
        case 'h':
            dohemi = !dohemi;
            break;
        case 's':
            dosnubify = !dosnubify;
            break;
        case 'q':
            dostellate = !dostellate;
            break;
        case 'e':
            explode += 1/16;
            break;
        case 'w':
            explode -= 1/16;
            break;
        case 'd':
            drawtype = (drawtype+1)%3
            break;
        case 'n':
            donormalize = !donormalize;
            break;
        case 'i':
            invertinc = -invertinc;
            break;
        case 'f':
            colorstyle = (colorstyle+1)%numcolorstyles;
            geometry.colorsNeedUpdate = true;
            break;
        case 't':
            tournum = (tournum+1)%tours.length;
            stopwatch.setTime(0);
            break;
        case 'u':
            rotatecolors();
            geometry.colorsNeedUpdate = true;
            break;
        case 'c':
            docompound = !docompound;
            break;
        case 'x':
            dorotonly = !dorotonly;
            break;
        case '=':
            tridepth = tridepth+1;
            break;
        case '-':
            if (tridepth > 0) tridepth = tridepth-1;
            break;
        case 'y':
            regionstyle = (regionstyle+1)%2;
            break;
        case '?':
            alert(makeinfostring());
            break;
        case '1':
            hideface[0] = !hideface[0];
            break;
        case '2':
            hideface[1] = !hideface[1];
            break;
        case '3':
            hideface[2] = !hideface[2];
            break;
        case '4':
            hideface[3] = !hideface[3];
            break;
        default:
            handled = false;
        }
        needupdate = true;
        if (handled) event.preventDefault();
    }
    document.addEventListener("keypress", onKeyPress, false);
    window.addEventListener("resize", function() {
        var w = window.innerWidth;
        var h = window.innerHeight;
        //console.log("resize " + w + " " + h);
        renderer.setSize(w,h);
        camera.aspect = w/h;
        camera.updateProjectionMatrix();
    });
    var render = function () { 
        requestAnimationFrame(render); 
        if (dozrot) {
            // Rotate the symmetry frame
            var zrotinc = 0.001;
            zrotation += zrotinc;
            needupdate = true;
        }
        var oldifact = ifact;
        ifact += invertinc;
        if (ifact > 1) ifact = 1;
        if (ifact < 0) ifact = 0;
        if (ifact != oldifact) needupdate = true;
        if (needupdate) {
            // Make sure we have the right color setting in material
            if (colorstyle >= vertexcolorstyle) {
                material.vertexColors = THREE.VertexColors;
            } else {
                material.vertexColors = THREE.FaceColors;
            }
            // Set up for compounding
            var ip = Vector.zrot(ip0,zrotation);
            var iq = Vector.zrot(iq0,zrotation);
            var ir = Vector.zrot(ir0,zrotation);
            var points = tours[tournum%tours.length];
            tri = maketri(stopwatch.getTime(),points);
            context = new Context(geometry,colorface);
            initcontext(context,schwarz,tri); // Make basic polyhedron
            drawpolyhedron(context);
            if (docompound) {
                // Now we do a transitive closure operation. Keep a list of sets of
                // points, pull sets out and apply all basic operations (from the
                // symmetry group), add back in any sets we haven't previously seen.
                // Initially we just know about the basic pointset just generated.
                var pointsets = [{trans:"", pointset: context.pointset}];
                var pindex = 0; // Where we are in the list of pointsets.
                var count = 0;
                while (pindex != pointsets.length && count < 200) {
                    count++;
                    var entry = pointsets[pindex++];
                    var ops = dorotonly ? rgens : gens; // Rotation only or all operations?
                    for (var i = 0; i < ops.length; i++) {
                        var newset = mapply(ops[i],entry.pointset,ip,iq,ir);
                        setsort(newset);
                        // Have we seen this pointset before?
                        var seen = false;
                        var parity = (entry.trans.length + ops[i].length)%2;
                        for (var j = 0; j < pointsets.length; j++) {
                            // If we care about chirality (ie. are generating a snub figure)
                            // then don't consider reflected pointsets the same).
                            if (dosnubify && pointsets[j].trans.length%2 != parity) continue;
                            if (setequal(newset,pointsets[j].pointset,1e-6)) {
                                seen = true;
                                break;
                            }
                        }
                        if (!seen) {
                            // A new compound element. Draw it properly and add to queue.
                            // TBD: we could just reuse the geometry with a rotation
                            // and a change of color.
                            var newtrans = entry.trans + ops[i];
                            context.transmat = makematrix(newtrans,ip,iq,ir);
                            context.compound++;
                            drawpolyhedron(context);
                            pointsets.push({trans: newtrans, pointset: newset});
                        }
                    }
                }
            }
            // Post process...
            geometry.computeFaceNormals();
            geometry.verticesNeedUpdate = true;
            geometry.elementsNeedUpdate = true;
            geometry.normalsNeedUpdate = true;
            geometry.colorsNeedUpdate = true;
            //geometry.buffersNeedUpdate = true; // If buffer lengths have changed
            if (context.needclone ||
                context.geometry.vertices.length != context.npoints ||
                context.geometry.faces.length != context.nfaces) {
                //console.log("Cloning geometry: " + context.npoints + " " + context.nfaces);
                context.geometry.vertices.length = context.npoints;
                context.geometry.faces.length = context.nfaces;
                var newgeometry = geometry.clone()
                scene.remove(mesh);
                geometry.dispose();
                geometry = newgeometry;
                // Reusing the mesh object seems OK
                mesh.geometry = geometry; 
                scene.add(mesh);
            }
            // Do we need to update again?
            if (!stopwatch.running) needupdate = false;
        }
        if (dorotate) {
            mesh.rotation.x += 0.002; 
            mesh.rotation.y += 0.002; 
        }
        renderer.render(scene, camera);
    };
    render(0); 
})()
