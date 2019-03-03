"use strict";

let initialize = null;

(function() {
    // TODO: get this from the URL
    let shaderID = null
    //shaderID = "4sX3Rn" // iq's menger sponge
    //shaderID = "lttfDH"
    //shaderID = "XdGczw"     // parallepiped
    //shaderID = "4tSBDz"     // inverted spheres


function Renderer(gl,url,oncompletion) {
    // Send off a request for the fragment shader code.
    let renderer = this;
    this.gl = gl;
    this.iMatrix = mat4.create();
    const request = new XMLHttpRequest();
    request.open("GET", url);
    request.onreadystatechange = function() {
        if (request.readyState === 4) {
            renderer.finalize(request.responseText, oncompletion);
        }
    }
    request.send(null); // No body
};

// TODO: better name
Renderer.prototype.finalize = function (fstext,oncompletion) {
    const VS = [
        "#version 300 es",
        "in vec3 aVertexPosition;",
        "out vec2 vTextureCoord;",
        "void main(void) {",
        "  gl_Position = vec4(aVertexPosition,1.0);",
        "  vTextureCoord = gl_Position.xy;",
        "}"
    ].join("\n");

    const FSpreamble = [
        "#version 300 es",
        "precision highp float;",
        "uniform float iTime;",
        "uniform mat4 iView;",
        "uniform mat4 iProjection;",
        "uniform mat4 iMatrix;",
        "uniform vec4 iResolution;",
        "uniform vec4 iMouse;",
        "out vec4 outColor;",
        "in vec2 vTextureCoord;",
        "void mainVR( out vec4 fragColor, in vec2 fragCoord,",
        "             in vec3 fragRayOrigin, in vec3 fragRayDir);",
        "void main() {",
        "  vec4 p = vec4(vTextureCoord,0,1);", // The "screen position", -1 <= z <= 1
        "  vec4 eye = vec4(0,0,1,0);",         // z-infinity
        "  mat4 m = iMatrix;",
        "  p = m*p;",
        "  eye = m*eye;",
        "  p /= p.w;",
        "  eye /= eye.w;",
        "  mainVR(outColor,vec2(0),eye.xyz,normalize(p.xyz-eye.xyz));",
        "}",
        "#define texelFetch(a,b,c) (vec4(0))", // TODO: better way of fixing this
        "#line 1",
        ""
    ].join("\n");
    let gl = this.gl;

    // Compile and link the shaders.
    function makeShader(source, shadertype) {
        const shader = gl.createShader(shadertype);
        gl.shaderSource(shader, source);
        gl.compileShader(shader);
        if (gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
            return shader;
        } else {
            let typestring = "<unknown>";
            if (shadertype == gl.VERTEX_SHADER) typestring = "vertex";
            else if (shadertype == gl.FRAGMENT_SHADER) typestring = "fragment";
            alert("Error for " + shadertype + " shader: " + gl.getShaderInfoLog(shader));
            return null;
        }
    }
    function initShaders(vshader,fshader) {
        let vertexShader = makeShader(vshader,gl.VERTEX_SHADER);
        let fragmentShader = makeShader(fshader,gl.FRAGMENT_SHADER);
        if (!vertexShader || !fragmentShader) return null;
        let program = gl.createProgram();
        gl.attachShader(program, vertexShader);
        gl.attachShader(program, fragmentShader);
        gl.linkProgram(program);
        gl.validateProgram(program); // Check all well
        // If creating the shader program failed, alert
        if (gl.getProgramParameter(program, gl.LINK_STATUS)) {
            return program;
        } else {
            alert("Unable to initialize program: " + gl.getProgramInfoLog(program));
            return null;
        }
    }

    if (shaderID) {
        // Pull data out from Shadertoy JSON response.
        // This is probably a bit fragile.
        let fsobj = JSON.parse(fstext);
        fstext = fsobj['Shader']['renderpass']['0']['code'];
    }
    this.program = initShaders(VS,FSpreamble + fstext);
    if (!this.program) return;

    // Two triangles fill the screen
    const vertices = [
        1.0, 1.0, 0.0,
        -1.0, 1.0, 0.0,
        1.0,-1.0, 0.0,
        -1.0,-1.0, 0.0
    ];
    this.vertBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, this.vertBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(vertices), gl.STATIC_DRAW);
    oncompletion();
}

Renderer.prototype.draw = function (projectionMatrix, viewMatrix, selected) {
    const program = this.program;
    if (!program) {
        console.log("Program not ready!");
        return;
    }
    const gl = this.gl;
    gl.useProgram(this.program);
    const iMatrix = this.iMatrix;
    mat4.mul(iMatrix,projectionMatrix,viewMatrix);
    mat4.invert(iMatrix,iMatrix);
    const index = gl.getAttribLocation(program, "aVertexPosition");
    if (index >= 0) {
        gl.bindBuffer(gl.ARRAY_BUFFER, this.vertBuffer);
        gl.enableVertexAttribArray(index);
        gl.vertexAttribPointer(index, 3, gl.FLOAT, false, 3*4, 0*4);
    }
    gl.uniform1f(gl.getUniformLocation(program, "iTime"), this.time);
    // Parameter 2 is `transpose` - set it to false!
    gl.uniformMatrix4fv(gl.getUniformLocation(program, "iView"), false, viewMatrix);
    gl.uniformMatrix4fv(gl.getUniformLocation(program, "iProjection"), false, projectionMatrix);
    gl.uniformMatrix4fv(gl.getUniformLocation(program, "iMatrix"), false, iMatrix);
    let devicePixelRatio = window.devicePixelRatio || 1; // What to do with this?
    let width = gl.canvas.width, height = gl.canvas.height;
    //let width = gl.canvas.drawingBufferWidth, height = gl.canvas.drawingBufferHeight;
    gl.uniform4f(gl.getUniformLocation(program, "iResolution"),
                 width*devicePixelRatio, height*devicePixelRatio,0,0);
    if (selected) {
        // TBD: think of a better way of doing this
        gl.uniform4f(gl.getUniformLocation(program, "iMouse"),1,1,0,0);
    } else {
        gl.uniform4f(gl.getUniformLocation(program, "iMouse"),0,0,0,0);
    }
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
    let err = gl.getError();
    if (err) console.log("GL Error: ", err);
};

Renderer.prototype.startFrame = function () {
    //console.log("startFrame");
    let t = Date.now();
    if (!this.startTime) this.startTime = t;
    this.time = (t-this.startTime)/1000;
}

Renderer.prototype.endFrame = function () {
    //console.log("endFrame");
}

// Code from cottontail, cut down to just WebGL2
// Creates a WebGL context and initializes it with some common default state.

function createWebGLContext(glAttribs) {
  glAttribs = glAttribs || {alpha: false}; // Why no alpha?

  let webglCanvas = document.createElement('canvas');
  let contextTypes = ['webgl2'];
  let context = null;

  for (let contextType of contextTypes) {
    context = webglCanvas.getContext(contextType, glAttribs);
    if (context) {
      break;
    }
  }

  if (!context) {
    let webglType = 'WebGL 2';
    console.error('This browser does not support ' + webglType + '.');
    return null;
  }

  return context;
}

initialize = function() {
    // Apply the version shim.
    var versionShim = new WebXRVersionShim();

    // XR globals.
    let xrButton = null;
    let xrImmersiveRefSpace = null;
    let xrNonImmersiveRefSpace = null;

    // WebGL scene globals.
    let gl = null;
    let renderer = null;

    let projectionMatrix = mat4.create();
    let viewMatrix = mat4.create();

    let running = true;
    let selected = false;

    let shaderurl;

    if (shaderID) {
        shaderurl = "https://www.shadertoy.com/api/v1/shaders/" + shaderID + "?key=fdntwh";
    } else {
        shaderurl = "goursat.glsl";
        shaderurl += "?" + new Date().getTime(); // Skip caching
    }

    // When we hit the fallback path, we'll need to initialize a few extra
    // variables in order to render correctly.
    function initFallback() {
        initGL(false, function() {
            document.body.appendChild(gl.canvas);

            // Using a simple identity matrix for the view.
            mat4.identity(viewMatrix);

            // Space bar is stop/start animation
            function keypressHandler(event) {
                if (!event.ctrlKey) {
                    // Ignore event if control key pressed.
                    var c = String.fromCharCode(event.charCode)
                    switch(c) {
                    case ' ':
                        running = !running;
                        if (running) {
                            // If we are now running, start animating.
                            window.requestAnimationFrame(onWindowFrame);
                        }
                        break;
                    }
                }
            }
            // We need to track the canvas size in order to resize the WebGL
            // backbuffer width and height, as well as update the projection matrix
            // and adjust the viewport.
            function onResize() {
                gl.canvas.width = gl.canvas.offsetWidth * window.devicePixelRatio;
                gl.canvas.height = gl.canvas.offsetHeight * window.devicePixelRatio;
                mat4.perspective(projectionMatrix, Math.PI*0.4,
                                 gl.canvas.width/gl.canvas.height,
                                 0.1, 1000.0);
                gl.viewport(0, 0, gl.drawingBufferWidth, gl.drawingBufferHeight);
            }
            window.addEventListener('resize', onResize);
            window.addEventListener("keypress",keypressHandler,false);
            onResize();

            // We'll kick off the render loop using the window's
            // requestAnimationFrame function.
            window.requestAnimationFrame(onWindowFrame);
        });
    }
              
    // Since both the XR and fallback paths need to do the same WebGL
    // initialization code, we've moved that code out to it's own function.
    function initGL(needXR, callback) {
        console.assert(!gl);

        gl = createWebGLContext({
            xrCompatible: needXR
        });

        renderer = new Renderer(gl,shaderurl,callback);
    }

    function onRequestSession() {
        let mirrorCanvas = document.createElement('canvas');
        let ctx = mirrorCanvas.getContext('xrpresent');
        mirrorCanvas.setAttribute('id', 'mirror-canvas');
        document.body.appendChild(mirrorCanvas);
        navigator.xr.requestSession({ mode: 'immersive-vr', outputContext: ctx }).then((session) => {
            xrButton.setSession(session);
            onSessionStarted(session);
        });
    }

    function onSelect() {
        selected = !selected;
    }
    
    function onSessionStarted(session) {
        session.addEventListener('end', onSessionEnded);
        session.addEventListener('select', onSelect);

        initGL(true, function() {
            session.baseLayer = new XRWebGLLayer(session, gl);

            session.requestReferenceSpace({ type: 'stationary', subtype: 'eye-level' }).then((refSpace) => {
                if (session.immersive) {
                    xrImmersiveRefSpace = refSpace;
                } else {
                    xrNonImmersiveRefSpace = refSpace;
                }
                session.requestAnimationFrame(onXRFrame);
            });
        });
    }

    function onEndSession(session) {
        session.end();
    }

    function onSessionEnded(event) {
        if (event.session.immersive) {
            document.body.removeChild(document.querySelector('#mirror-canvas'));
            xrButton.setSession(null);
        }
    }

    function onXRFrame(timestamp, frame) {
        let session = frame.session;
        let refSpace = session.immersive ?
            xrImmersiveRefSpace :
            xrNonImmersiveRefSpace;
        let pose = frame.getViewerPose(refSpace);

        session.requestAnimationFrame(onXRFrame);

        renderer.startFrame();

        if (pose) {
            gl.bindFramebuffer(gl.FRAMEBUFFER, session.baseLayer.framebuffer);
            gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

            for (let view of pose.views) {
                let viewport = session.baseLayer.getViewport(view);
                gl.viewport(viewport.x, viewport.y,
                            viewport.width, viewport.height);
                // Pass timestamp?
                renderer.draw(view.projectionMatrix, view.viewMatrix, selected);
            }
        }

        renderer.endFrame();
    }

    // This is the bulk of the fallback rendering loop. Notice that it looks
    // a lot like a simplified version of the XR frame loop. Samples after
    // this one will do some work to hide this for code readability purposes.
    function onWindowFrame(t) {
        if (running) window.requestAnimationFrame(onWindowFrame);

        renderer.startFrame();

        // We can skip setting the framebuffer and viewport every frame, because
        // it won't change from frame to frame and we're updating the viewport
        // only when we resize for efficency.
        gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

        // We're drawing with our own projection and view matrix now, and we 
        // don't have a list of view to loop through, but otherwise all of the
        // WebGL drawing logic is exactly the same.

        renderer.draw(projectionMatrix, viewMatrix, false);
        renderer.endFrame();
    }
    xrButton = new XRDeviceButton({
        onRequestSession: onRequestSession,
        onEndSession: onEndSession
    });
    document.querySelector('header').appendChild(xrButton.domElement);

    // Find browser capabilities
    if (navigator.xr) {
        navigator.xr.supportsSessionMode('immersive-vr').then(() => {
            xrButton.enabled = true;
        });

        let outputCanvas = document.createElement('canvas');
        let ctx = outputCanvas.getContext('xrpresent');

        navigator.xr.requestSession({ outputContext: ctx })
            .then((session) => {
                document.body.appendChild(outputCanvas);
                onSessionStarted(session);
            })
            .catch((err) => {
                initFallback();
            });
    } else {
        // If navigator.xr isn't present in the browser then we need to use
        // the fallback rendering path.
        initFallback();
    }
}
})();
