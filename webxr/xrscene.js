"use strict"

function Renderer(gl,callback) {
    // Fire off a request for the fragment shader code.
    let renderer = this;
    this.gl = gl;
    let url = "goursat.frag";
    url += "?" + new Date().getTime();
    const request = new XMLHttpRequest();
    request.open("GET", url);
    request.onreadystatechange = function() {
        if (request.readyState === 4) {
            if (renderer.finalize(request.responseText)) {
                callback();
            }
        }
    }
    request.send(null); // No body
};

Renderer.prototype.finalize = function (fstext) {
    const VS = [
        "#version 300 es",
        "in vec3 aVertexPosition;",
        "out vec2 vTextureCoord;",
        "void main(void) {",
        "  gl_Position = vec4(aVertexPosition,1.0);",
        "  vTextureCoord = gl_Position.xy;",
        "}"
    ].join("\n");

    // Not using this
    const FS = [
        "#version 300 es",
        "precision highp float;",
        "uniform float iTime;",
        "uniform mat4 iView;",
        "out vec4 outColor;",
        //"uniform sampler2D diffuse;",
        "in vec2 vTextureCoord;",

        "void main() {",
        "  vec4 p = vec4(vTextureCoord,0,1);",
        "  p = inverse(iView)*p;",
        "  if (length(p.xy) < 0.5) outColor = vec4(1,0.5+0.5*sin(iTime),0,1);",
        "  else outColor = vec4(0,1,0,1);",
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
        //"  vec4 p = vec4(gl_FragCoord.xy/iResolution.xy*2.0-1.0,0,1);",
        "  vec4 p = vec4(vTextureCoord,0,1);",
        "  vec4 eye = vec4(0,0,1,0);",
        "  mat4 m = iMatrix; //inverse(iProjection*iView);",
        "  p = m*p;",
        "  eye = m*eye;",
        "  p /= p.w;",  // Don't know we really need this
        "  eye /= eye.w;", // Or this
        "  mainVR(outColor,vec2(0),eye.xyz,normalize(p.xyz-eye.xyz));",
        "}",
        "#define texelFetch(a,b,c) (vec4(0))",
        "#line 1",
        ""
    ].join("\n");
    let gl = this.gl;

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
    // Initialize the shaders.
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

    this.program = initShaders(VS,FSpreamble + fstext);
    if (!this.program) return false;
    
    const vertices = [
        1.0, 1.0, 0.0,
        -1.0, 1.0, 0.0,
        1.0,-1.0, 0.0,
        -1.0,-1.0, 0.0
    ];
    this.vertBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, this.vertBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(vertices), gl.STATIC_DRAW);
    this.startTime = Date.now();
    return true;
}

Renderer.prototype.draw = function (projectionMatrix, viewMatrix) {
    const gl = this.gl;
    const program = this.program;
    if (!program) {
        console.log("Not ready!");
        return;
    }
    let iMatrix = mat4.clone(projectionMatrix);
    mat4.mul(iMatrix,iMatrix,viewMatrix);
    mat4.invert(iMatrix,iMatrix);
    gl.useProgram(this.program);
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
    gl.uniform4f(gl.getUniformLocation(program, "iResolution"), gl.canvas.width, gl.canvas.height,0,0);
    gl.uniform4f(gl.getUniformLocation(program, "iMouse"),0,0,0,0);
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
    let err = gl.getError();
    if (err) console.log("GL Error: ", err);
};

Renderer.prototype.startFrame = function () {
    //console.log("startFrame");
    this.time = (Date.now()-this.startTime)/1000;
}

Renderer.prototype.endFrame = function () {
    //console.log("endFrame");
}

// Nicked from cottontail - cut down to just WebGL2?
// Creates a WebGL context and initializes it with some common default state.
function createWebGLContext(glAttribs) {
  glAttribs = glAttribs || {alpha: false};

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
    let webglType = (glAttribs.webgl2 ? 'WebGL 2' : 'WebGL');
    console.error('This browser does not support ' + webglType + '.');
    return null;
  }

  return context;
}

function initialize() {
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
    
    function initXR() {
        xrButton = new XRDeviceButton({
            onRequestSession: onRequestSession,
            onEndSession: onEndSession
        });
        document.querySelector('header').appendChild(xrButton.domElement);

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

    // When we hit the fallback path, we'll need to initialize a few extra
    // variables in order to render correctly.
    function initFallback() {
        initGL(false, function() {
            document.body.appendChild(gl.canvas);

            // Using a simple identity matrix for the view.
            mat4.identity(viewMatrix);

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

        renderer = new Renderer(gl,callback);
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

    function onSessionStarted(session) {
        session.addEventListener('end', onSessionEnded);

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

    function onXRFrame(t, frame) {
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

                renderer.draw(view.projectionMatrix, view.viewMatrix);
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

        renderer.draw(projectionMatrix, viewMatrix);
        renderer.endFrame();
    }

    // Start the XR application.
    initXR();
}
