
/**
 * @file A simple WebGL example drawing central Illinois style terrain
 * @author Jerry Chang <jerryc2@illinois.edu> 
 */

/** @global The WebGL context */
var gl;

/** @global The HTML5 canvas we draw on */
var canvas;

/** @global A simple GLSL shader program */
var shaderProgram;

/** @global The Modelview matrix */
var mvMatrix = glMatrix.mat4.create();

/** @global The Projection matrix */
var pMatrix = glMatrix.mat4.create();

/** @global The Normal matrix */
var nMatrix = glMatrix.mat3.create();

/** @global The matrix stack for hierarchical modeling */
var mvMatrixStack = [];

/** @global The angle of rotation around the y axis */
var viewRot = 10;

/** @global A glmatrix vector to use for transformations */
var transformVec = glMatrix.vec3.create();    

// Initialize the vector....
glMatrix.vec3.set(transformVec,0.0,0.0,-2.0);

/** @global An object holding the geometry for a 3D terrain */
var myTerrain;


// View parameters
/** @global Location of the camera in world coordinates */
var eyePt = glMatrix.vec3.fromValues(0.0,0.0,0.0);
/** @global Direction of the view in world coordinates */
var viewDir = glMatrix.vec3.fromValues(0.0,0.0,-1.0);
/** @global Up vector for view matrix creation, in world coordinates */
var up = glMatrix.vec3.fromValues(0.0,1.0,0.0);
/** @global Location of a point along viewDir in world coordinates */
var viewPt = glMatrix.vec3.fromValues(0.0,0.0,0.0);

//Light parameters
/** @global Light position in VIEW coordinates */
var lightPosition = [0,3,3];
/** @global Ambient light color/intensity for Phong reflection */
var lAmbient = [0,0,0];
/** @global Diffuse light color/intensity for Phong reflection */
var lDiffuse = [1,1,1];
/** @global Specular light color/intensity for Phong reflection */
var lSpecular =[1,1,1];

//Material parameters
/** @global Ambient material color/intensity for Phong reflection */
var kAmbient = [1.0,1.0,1.0];
/** @global Diffuse material color/intensity for Phong reflection */
var kTerrainDiffuse = [205.0/255.0,163.0/255.0,63.0/255.0];
/** @global Specular material color/intensity for Phong reflection */
var kSpecular = [0.5,0.5, 0.5];
/** @global Shininess exponent for Phong reflection */
var shininess = 23;
/** @global Edge color fpr wireframeish rendering */
var kEdgeBlack = [0.0,0.0,0.0];
/** @global Edge color for wireframe rendering */
var kEdgeWhite = [1.0,1.0,1.0];



//-------------------------------------------------------------------------
/**
 * Sends Modelview matrix to shader
 */
function uploadModelViewMatrixToShader() {
  gl.uniformMatrix4fv(shaderProgram.mvMatrixUniform, false, mvMatrix);
}

//-------------------------------------------------------------------------
/**
 * Sends projection matrix to shader
 */
function uploadProjectionMatrixToShader() {
  gl.uniformMatrix4fv(shaderProgram.pMatrixUniform, 
                      false, pMatrix);
}

//-------------------------------------------------------------------------
/**
 * Generates and sends the normal matrix to the shader
 */
function uploadNormalMatrixToShader() {
  glMatrix.mat3.fromMat4(nMatrix,mvMatrix);
  glMatrix.mat3.transpose(nMatrix,nMatrix);
  glMatrix.mat3.invert(nMatrix,nMatrix);
  gl.uniformMatrix3fv(shaderProgram.nMatrixUniform, false, nMatrix);
}

//----------------------------------------------------------------------------------
/**
 * Pushes matrix onto modelview matrix stack
 */
function mvPushMatrix() {
    var copy = glMatrix.mat4.clone(mvMatrix);
    mvMatrixStack.push(copy);
}


//----------------------------------------------------------------------------------
/**
 * Pops matrix off of modelview matrix stack
 */
function mvPopMatrix() {
    if (mvMatrixStack.length == 0) {
      throw "Invalid popMatrix!";
    }
    mvMatrix = mvMatrixStack.pop();
}

//----------------------------------------------------------------------------------
/**
 * Sends projection/modelview matrices to shader
 */
function setMatrixUniforms() {
    uploadModelViewMatrixToShader();
    uploadNormalMatrixToShader();
    uploadProjectionMatrixToShader();
}

//----------------------------------------------------------------------------------
/**
 * Translates degrees to radians
 * @param {Number} degrees Degree input to function
 * @return {Number} The radians that correspond to the degree input
 */
function degToRad(degrees) {
        return degrees * Math.PI / 180;
}

//----------------------------------------------------------------------------------
/**
 * Creates a context for WebGL
 * @param {element} canvas WebGL canvas
 * @return {Object} WebGL context
 */
function createGLContext(canvas) {
  var names = ["webgl", "experimental-webgl"];
  var context = null;
  for (var i=0; i < names.length; i++) {
    try {
      context = canvas.getContext(names[i]);
    } catch(e) {}
    if (context) {
      break;
    }
  }
  if (context) {
    context.viewportWidth = canvas.width;
    context.viewportHeight = canvas.height;
  } else {
    alert("Failed to create WebGL context!");
  }
  return context;
}

//----------------------------------------------------------------------------------
/**
 * Loads Shaders
 * @param {string} id ID string for shader to load. Either vertex shader/fragment shader
 */
function loadShaderFromDOM(id) {
  var shaderScript = document.getElementById(id);
  
  // If we don't find an element with the specified id
  // we do an early exit 
  if (!shaderScript) {
    return null;
  }
  
  // Loop through the children for the found DOM element and
  // build up the shader source code as a string
  var shaderSource = "";
  var currentChild = shaderScript.firstChild;
  while (currentChild) {
    if (currentChild.nodeType == 3) { // 3 corresponds to TEXT_NODE
      shaderSource += currentChild.textContent;
    }
    currentChild = currentChild.nextSibling;
  }
 
  var shader;
  if (shaderScript.type == "x-shader/x-fragment") {
    shader = gl.createShader(gl.FRAGMENT_SHADER);
  } else if (shaderScript.type == "x-shader/x-vertex") {
    shader = gl.createShader(gl.VERTEX_SHADER);
  } else {
    return null;
  }
 
  gl.shaderSource(shader, shaderSource);
  gl.compileShader(shader);
 
  if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
    alert(gl.getShaderInfoLog(shader));
    return null;
  } 
  return shader;
}

//----------------------------------------------------------------------------------
/**
 * Setup the fragment and vertex shaders
 */
function setupShaders() {
  vertexShader = loadShaderFromDOM("shader-vs");
  fragmentShader = loadShaderFromDOM("shader-fs");
  
  shaderProgram = gl.createProgram();
  gl.attachShader(shaderProgram, vertexShader);
  gl.attachShader(shaderProgram, fragmentShader);
  gl.linkProgram(shaderProgram);

  if (!gl.getProgramParameter(shaderProgram, gl.LINK_STATUS)) {
    alert("Failed to setup shaders");
  }

  gl.useProgram(shaderProgram);

  shaderProgram.vertexPositionAttribute = gl.getAttribLocation(shaderProgram, "aVertexPosition");
  gl.enableVertexAttribArray(shaderProgram.vertexPositionAttribute);

  shaderProgram.vertexNormalAttribute = gl.getAttribLocation(shaderProgram, "aVertexNormal");
  gl.enableVertexAttribArray(shaderProgram.vertexNormalAttribute);

  shaderProgram.mvMatrixUniform = gl.getUniformLocation(shaderProgram, "uMVMatrix");
  shaderProgram.pMatrixUniform = gl.getUniformLocation(shaderProgram, "uPMatrix");
  shaderProgram.nMatrixUniform = gl.getUniformLocation(shaderProgram, "uNMatrix");
  shaderProgram.uniformLightPositionLoc = gl.getUniformLocation(shaderProgram, "uLightPosition");    
  shaderProgram.uniformAmbientLightColorLoc = gl.getUniformLocation(shaderProgram, "uAmbientLightColor");  
  shaderProgram.uniformDiffuseLightColorLoc = gl.getUniformLocation(shaderProgram, "uDiffuseLightColor");
  shaderProgram.uniformSpecularLightColorLoc = gl.getUniformLocation(shaderProgram, "uSpecularLightColor");
  shaderProgram.uniformShininessLoc = gl.getUniformLocation(shaderProgram, "uShininess");    
  shaderProgram.uniformAmbientMaterialColorLoc = gl.getUniformLocation(shaderProgram, "uAmbientMaterialColor");  
  shaderProgram.uniformDiffuseMaterialColorLoc = gl.getUniformLocation(shaderProgram, "uDiffuseMaterialColor");
  shaderProgram.uniformSpecularMaterialColorLoc = gl.getUniformLocation(shaderProgram, "uSpecularMaterialColor");
  shaderProgram.uFog = gl.getUniformLocation(shaderProgram, "uFog");    //uFog is a html variable in shader program
}

//-------------------------------------------------------------------------
/**
 * Sends material information to the shader
 * @param {Float32} alpha shininess coefficient
 * @param {Float32Array} a Ambient material color
 * @param {Float32Array} d Diffuse material color
 * @param {Float32Array} s Specular material color
 */
function setMaterialUniforms(alpha,a,d,s) {
  gl.uniform1f(shaderProgram.uniformShininessLoc, alpha);
  gl.uniform3fv(shaderProgram.uniformAmbientMaterialColorLoc, a);
  gl.uniform3fv(shaderProgram.uniformDiffuseMaterialColorLoc, d);
  gl.uniform3fv(shaderProgram.uniformSpecularMaterialColorLoc, s);
}

//-------------------------------------------------------------------------
/**
 * Sends light information to the shader
 * @param {Float32Array} loc Location of light source
 * @param {Float32Array} a Ambient light strength
 * @param {Float32Array} d Diffuse light strength
 * @param {Float32Array} s Specular light strength
 */
function setLightUniforms(loc,a,d,s) {
  gl.uniform3fv(shaderProgram.uniformLightPositionLoc, loc);
  gl.uniform3fv(shaderProgram.uniformAmbientLightColorLoc, a);
  gl.uniform3fv(shaderProgram.uniformDiffuseLightColorLoc, d);
  gl.uniform3fv(shaderProgram.uniformSpecularLightColorLoc, s);
}

//----------------------------------------------------------------------------------
/**
 * Populate buffers with data
 */
function setupBuffers() {
    myTerrain = new Terrain(64,-0.5,0.5,-0.5,0.5);
    myTerrain.loadBuffers();
}

//----------------------------------------------------------------------------------
/**
 * Draw call that applies matrix transformations to model and draws model in frame
 */
function draw() { 
    //console.log("function draw()")
    var transformVec = glMatrix.vec3.create();
  
    gl.viewport(0, 0, gl.viewportWidth, gl.viewportHeight);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

    // We'll use perspective 
    glMatrix.mat4.perspective(pMatrix,degToRad(45), 
                     gl.viewportWidth / gl.viewportHeight,
                     0.1, 200.0);

    glMatrix.vec3.transformQuat(r_up, up, orient);
    glMatrix.vec3.transformQuat(r_viewDir, viewDir, orient);

    // We want to look down -z, so create a lookat point in that direction    
    glMatrix.vec3.add(viewPt, eyePt, r_viewDir);
    // Then generate the lookat matrix and initialize the MV matrix to that view
    glMatrix.mat4.lookAt(mvMatrix,eyePt,viewPt,r_up);    
 
    //Draw Terrain
    mvPushMatrix();
    glMatrix.vec3.set(transformVec,0.0,-0.25,-2.0);
    glMatrix.mat4.translate(mvMatrix, mvMatrix,transformVec);
    glMatrix.mat4.rotateY(mvMatrix, mvMatrix, degToRad(viewRot));
    glMatrix.mat4.rotateX(mvMatrix, mvMatrix, degToRad(-75));
    setMatrixUniforms();
    setLightUniforms(lightPosition,lAmbient,lDiffuse,lSpecular);
    
    if ((document.getElementById("polygon").checked) || (document.getElementById("wirepoly").checked))
    { 
      setMaterialUniforms(shininess,kAmbient,kTerrainDiffuse,kSpecular); 
      myTerrain.drawTriangles();
    }
    
    if(document.getElementById("wirepoly").checked)
    {
      setMaterialUniforms(shininess,kAmbient,kEdgeBlack,kSpecular);
      myTerrain.drawEdges();
    }

    if(document.getElementById("wireframe").checked)
    {
      setMaterialUniforms(shininess,kAmbient,kEdgeWhite,kSpecular);
      myTerrain.drawEdges();
    }
    if(document.getElementById("fog").checked)      //set uFog to true if checked on screen, else false
    {
      gl.uniform1i(shaderProgram.uFog, true);
    }else{
      gl.uniform1i(shaderProgram.uFog, false);
    }
    mvPopMatrix();

  
}
var currentlyPressedKeys = {};
let rotation = glMatrix.quat.create();    //roation matrix
let orient = glMatrix.quat.create();

let r_up = glMatrix.vec3.clone(up);      //create rotated_up matrix
let r_viewDir = glMatrix.vec3.clone(viewDir);               //create rotated_viewDir Matrix
//----------------------------------------------------------------------------------
/**
 * sense a keypress
 */
function handlekeydown(event) {
  console.log("Key down ", event.key, " code ", event.code);

  if(event.key == "ArrowDown" || event.key== "ArrowUp" || event.key== "ArrowRight" || event.key== "ArrowLeft"){
    event.preventDefault();
  }
  currentlyPressedKeys[event.key] = true;
}

//----------------------------------------------------------------------------------
/**
 * sense a key release
 */
function handlekeyup() {
  console.log("Key up ", event.key, " code ", event.code);
  currentlyPressedKeys[event.key] = false;
}

let speed = 0.001;
//----------------------------------------------------------------------------------
/**
 * Animation function
 */
function animate() {
  let change = glMatrix.vec3.create();
  change = glMatrix.vec3.fromValues(speed*r_viewDir[0],speed*r_viewDir[1],speed*r_viewDir[2]);
  glMatrix.vec3.add(eyePt, eyePt, change);

  if(currentlyPressedKeys["ArrowRight"]){             //rotate right
    let right = glMatrix.vec3.create();
    glMatrix.vec3.cross(right, r_viewDir, r_up);
    r_up = glMatrix.vec3.fromValues(0.0,1.0,0.0);
    r_viewDir = glMatrix.vec3.fromValues(0.0,0.0,-1.0);
    glMatrix.vec3.transformQuat(r_up, r_up, orient);
    glMatrix.vec3.transformQuat(r_viewDir, r_viewDir, orient);
    
    glMatrix.quat.rotationTo(rotation, r_up, right);         //rotation quaternion
    glMatrix.quat.pow(rotation, rotation, 0.007);    //smaller rotation quaternion
    glMatrix.quat.multiply(orient, orient, rotation);
  }
  if(currentlyPressedKeys["ArrowLeft"]){        //rotate left
    let left = glMatrix.vec3.create();
    glMatrix.vec3.cross(left, r_viewDir, r_up);
    r_up = glMatrix.vec3.fromValues(0.0,1.0,0.0);
    r_viewDir = glMatrix.vec3.fromValues(0.0,0.0,-1.0);
    glMatrix.vec3.transformQuat(r_up, r_up, orient);
    glMatrix.vec3.transformQuat(r_viewDir, r_viewDir, orient);
    
    glMatrix.quat.rotationTo(rotation, left, r_up);         //rotation quaternion
    glMatrix.quat.pow(rotation, rotation, 0.007);    //smaller rotation quaternion
    glMatrix.quat.multiply(orient, orient, rotation);
  }
  if(currentlyPressedKeys["ArrowUp"]){                      //rotate up
    r_up = glMatrix.vec3.fromValues(0.0,1.0,0.0);
    r_viewDir = glMatrix.vec3.fromValues(0.0,0.0,-1.0);
    glMatrix.vec3.transformQuat(r_up, r_up, orient);
    glMatrix.vec3.transformQuat(r_viewDir, r_viewDir, orient);
    
    glMatrix.quat.rotationTo(rotation, r_viewDir, r_up);         //rotation quaternion
    glMatrix.quat.pow(rotation, rotation, 0.007);    //smaller rotation quaternion
    glMatrix.quat.multiply(orient, orient, rotation);
  }
  if(currentlyPressedKeys["ArrowDown"]){                  //rotate down
    r_up = glMatrix.vec3.fromValues(0.0,1.0,0.0);
    r_viewDir = glMatrix.vec3.fromValues(0.0,0.0,-1.0);
    glMatrix.vec3.transformQuat(r_up, r_up, orient);
    glMatrix.vec3.transformQuat(r_viewDir, r_viewDir, orient);
    
    glMatrix.quat.rotationTo(rotation, r_up, r_viewDir);         //rotation quaternion
    glMatrix.quat.pow(rotation, rotation, 0.007);    //smaller rotation quaternion
    glMatrix.quat.multiply(orient, orient, rotation);
  }
  if(currentlyPressedKeys["-"]){                                    //decrease speed
      speed -= 0.0002;
  }
  if(currentlyPressedKeys["+"]){                                  //increase speed
    speed +=0.0002;
  }

}

//----------------------------------------------------------------------------------
/**
 * Startup function called from html code to start program.
 */
 function startup() {
  canvas = document.getElementById("myGLCanvas");
  gl = createGLContext(canvas);
  setupShaders();
  setupBuffers();
  gl.clearColor(1.0, 1.0, 1.0, 1.0);
  gl.enable(gl.DEPTH_TEST);
  document.onkeydown = handlekeydown;
  document.onkeyup = handlekeyup;
  tick();
}

//----------------------------------------------------------------------------------
/**
 * Keeping drawing frames....
 */
function tick() {
    animate();
    requestAnimFrame(tick);
    draw();
}


