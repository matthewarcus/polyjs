"use strict"

const VR_MODE = 0;
const AR_MODE = 1;
function xrscene(mode) {
    if (mode != VR_MODE && mode != AR_MODE) {
        throw new Error('Invalid XR mode')
    }
    let sessiontype =
        mode == VR_MODE ? 'immersive-vr':
        mode == AR_MODE ? 'immersive-ar' :
        null

    let referencespace = 'local';
    
    let shaderfile = null
    shaderfile = "goursat.glsl";
    let shaderID = "llKfDc";
    shaderID = "4sX3Rn" // iq's menger sponge
    //shaderID = "lttfDH"
    //shaderID = "XdGczw"     // parallepiped
    //shaderID = "4tSBDz"     // inverted spheres

    const xrButton = document.getElementById('xr-button');
    let xrSession = null;
    let xrRefSpace = null;
    let gl = null;
    let renderer = null;

    function makeShader(source, shadertype) {
        const shader = gl.createShader(shadertype);
        gl.shaderSource(shader, source);
        gl.compileShader(shader);
        const infolog = gl.getShaderInfoLog(shader);
        if (gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
            if (infolog) alert("Info for " + shadertype + " shader: " + infolog);
            return shader;
        } else {
            let typestring = "<unknown>";
            if (shadertype == gl.VERTEX_SHADER) typestring = "vertex";
            else if (shadertype == gl.FRAGMENT_SHADER) typestring = "fragment";
            alert("Error for " + shadertype + " shader: " + infolog);
        }
    }

    function Renderer(fsource) {
        this.iMatrix = mat4.create();

        // Compile and link the shaders.
        const vsource = [
            "#version 300 es",
            "in vec3 aVertexPosition;",
            "out vec2 vTextureCoord;",
            "void main(void) {",
            "  gl_Position = vec4(aVertexPosition,1.0);",
            "  vTextureCoord = gl_Position.xy;",
            "}"
        ].join("\n");

        const preamble = [
            "#version 300 es",
            "precision highp float;",
            // "uniform mat4 iView;", // Define in shader, if needed
            // "uniform mat4 iProjection;", 
            "uniform mat4 iMatrix;",
            "uniform float iTime;",
            "uniform int iFrame;",
            "uniform vec4 iResolution;",
            "uniform vec4 iMouse;",
            "out vec4 outColor;",
            "in vec2 vTextureCoord;",
            "void mainVR( out vec4 fragColor, vec4 eye, vec4 screenpos);",
            "void main() {",
            "  vec4 screenpos = vec4(vTextureCoord,0,1);", // The "screen position", -1 <= z <= 1
            "  vec4 eye = vec4(0,0,1,0);",         // z-infinity
            "  mat4 m = iMatrix;",
            "  mainVR(outColor,m*eye,m*screenpos);",
            "}",
            "#define HW_PERFORMANCE 0", // Some Shadertoys use this to set AA
            "#line 1",
            ""
        ].join("\n");

        if (!shaderfile) {
            // Pull data out from Shadertoy JSON response.
            // This is probably a bit fragile.
            const json = JSON.parse(fsource);
            let error = json['Error'];
            if (error) throw new Error(error);
            const shader = json['Shader']
            console.log(shader.info);
            fsource = "";
            for (const pass of shader.renderpass) {
                if (pass.type == 'common') fsource = pass.code+fsource;
                if (pass.type == 'image') fsource = fsource+pass.code;
            }
            if (!fsource) throw new Error('No shader source found');
        }
        fsource = preamble + fsource;
        const vshader = makeShader(vsource,gl.VERTEX_SHADER);
        const fshader = makeShader(fsource,gl.FRAGMENT_SHADER);
        if (!vshader || !fshader) throw new Error("compilation failed");
        const program = gl.createProgram();
        if (!program) throw new Error("program creation failed");
        gl.attachShader(program, vshader);
        gl.attachShader(program, fshader);
        gl.linkProgram(program);
        gl.validateProgram(program); // Check all well
        if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
            throw(new Error("Unable to initialize program: " + gl.getProgramInfoLog(program)));
        }

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
        this.program = program;
    }

    const devicePixelRatio = window.devicePixelRatio || 1; // This is 3 on my phone
    //if (devicePixelRatio != 1) alert("devicePixelRatio: " + devicePixelRatio);

    let framecounter = 0;
    Renderer.prototype.draw = function (projectionMatrix, viewMatrix, timestamp) {
        //console.log(projectionMatrix);
        //console.log(viewMatrix);
        if (this.start == null) this.start = timestamp;
        const iMatrix = this.iMatrix;
        const program = this.program;
        console.assert(iMatrix);
        console.assert(program);
        gl.useProgram(program);
        mat4.copy(iMatrix,projectionMatrix);
        mat4.mul(iMatrix,iMatrix,viewMatrix);
        mat4.invert(iMatrix,iMatrix);
        const index = gl.getAttribLocation(program, "aVertexPosition");
        if (index >= 0) {
            gl.bindBuffer(gl.ARRAY_BUFFER, this.vertBuffer);
            gl.enableVertexAttribArray(index);
            gl.vertexAttribPointer(index, 3, gl.FLOAT, false, 3*4, 0*4);
        }
        gl.uniform1i(gl.getUniformLocation(program, "iFrame"), framecounter++);
        gl.uniform1f(gl.getUniformLocation(program, "iTime"), (timestamp-this.start)/1000);
        // Parameter 2 is `transpose` - set it to false!
        gl.uniformMatrix4fv(gl.getUniformLocation(program, "iView"), false, viewMatrix);
        gl.uniformMatrix4fv(gl.getUniformLocation(program, "iProjection"), false, projectionMatrix);
        gl.uniformMatrix4fv(gl.getUniformLocation(program, "iMatrix"), false, iMatrix);
        // drawingBufferWidth/Height should be the right thing
        const width = gl.canvas.drawingBufferWidth, height = gl.canvas.drawingBufferHeight;
        gl.uniform4f(gl.getUniformLocation(program, "iResolution"),
                     width*devicePixelRatio, height*devicePixelRatio,0,0);
        gl.uniform4f(gl.getUniformLocation(program, "iMouse"),0,0,0,0);
        gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
        const err = gl.getError();
        if (err) console.log("GL Error: ", err);
    };

    // Called either when the user has explicitly ended the session by calling
    // session.end() or when the UA has ended the session for any reason.
    // At this point the session object is no longer usable and should be
    // discarded.
    function onSessionEnded(event) {
        console.log("WebXR session ended");
        xrSession = null;
        xrButton.textContent = 'Enter VR';

        // In this simple case discard the WebGL context too, since we're not
        // rendering anything else to the screen with it.
        gl = null;
    }
    
    function onXRFrame(timestamp, frame) {
        const session = frame.session;
        const pose = frame.getViewerPose(xrRefSpace);
        session.requestAnimationFrame(onXRFrame);
        if (pose) {
            const baseLayer = session.renderState.baseLayer;
            gl.bindFramebuffer(gl.FRAMEBUFFER, baseLayer.framebuffer);
            gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

            for (const view of pose.views) {
                const viewport = baseLayer.getViewport(view);
                gl.viewport(viewport.x, viewport.y,
                            viewport.width, viewport.height);
                // Pass timestamp?
                //console.log(view.transform);
                renderer.draw(view.projectionMatrix, view.transform.inverse.matrix, timestamp);
            }
        }
    }

    function makeshaderurl(shaderID) {
        const qparam = new Date().getTime();  // Skip caching
        if (shaderfile) return shaderfile + "?" + qparam;
        else if (shaderID) return "https://www.shadertoy.com/api/v1/shaders/" + shaderID + "?key=fdntwh&" + qparam;
        else throw new Error('No shader specified');
    }
        
    function onSessionStarted (session) {
        xrSession = session;
        xrButton.textContent = 'Exit VR';
        session.addEventListener('end', onSessionEnded);
        const webglCanvas = document.createElement('canvas');
        gl = webglCanvas.getContext('webgl2', { xrCompatible: true })

        // Send off a request for the fragment shader code.
        const request = new XMLHttpRequest();
        request.open("GET", makeshaderurl(shaderID));
        request.onreadystatechange = function() {
            //console.log(request.readyState);
            if (request.readyState === 4) {
                function setRefSpace(refSpace) {
                    xrRefSpace = refSpace;
                    session.requestAnimationFrame(onXRFrame);
                }
                renderer = new Renderer(request.responseText);
                session.updateRenderState({ baseLayer: new XRWebGLLayer(session, gl) });
                session.requestReferenceSpace(referencespace).
                    then(setRefSpace,
                         error => {
                             session.requestReferenceSpace('local').then(setRefSpace);
                         });
            }
        }
        request.send(null); // No body
    }

    // Called when the user clicks the button to enter XR. If we don't have a
    // session we'll request one, and if we do have a session we'll end it.
    function onButtonClicked() {
        if (xrSession) xrSession.end();
        else navigator.xr.requestSession(sessiontype).then(onSessionStarted);
    }

    // Find browser capabilities
    if (navigator.xr) {
        navigator.xr.isSessionSupported(sessiontype).then((supported) => {
            if (supported) {
                xrButton.addEventListener('click', onButtonClicked);
                xrButton.textContent = 'Enter VR';
                xrButton.disabled = false;
            }
        });
    }
}
