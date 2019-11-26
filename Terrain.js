/**
 * @fileoverview Terrain - A simple 3D terrain using WebGL
 * Author: Jerry Chang <jerryc2.illinois.edu>
 */

/** Class implementing 3D terrain. */
class Terrain{   
/**
 * Initialize members of a Terrain object
 * @param {number} div Number of triangles along x axis and y axis
 * @param {number} minX Minimum X coordinate value
 * @param {number} maxX Maximum X coordinate value
 * @param {number} minY Minimum Y coordinate value
 * @param {number} maxY Maximum Y coordinate value
 */
    constructor(div,minX,maxX,minY,maxY){
        this.div = div;
        this.minX=minX;
        this.minY=minY;
        this.maxX=maxX;
        this.maxY=maxY;
        // Allocate vertex array
        this.vBuffer = [];
        // Allocate triangle array
        this.fBuffer = [];
        // Allocate normal array
        this.nBuffer = [];
        this.count = [];
        //Allocate x-buffer
        //this.zBuffer = [];
        // Allocate array for edges so we can draw wireframe
        this.eBuffer = [];
        console.log("Terrain: Allocated buffers");

        this.generateTriangles();
        console.log("Terrain: Generated triangles");

        this.z_generation();

        this.avg_normals();

        this.generateLines();
        console.log("Terrain: Generated lines");
        
        // Get extension for 4 byte integer indices for drwElements
        var ext = gl.getExtension('OES_element_index_uint');
        if (ext ==null){
            alert("OES_element_index_uint is unsupported by your browser and terrain generation cannot proceed.");
        }
    }
    
    /**
    * Set the x,y,z coords of a vertex at location(i,j)
    * @param {Object} v an an array of length 3 holding x,y,z coordinates
    * @param {number} i the ith row of vertices
    * @param {number} j the jth column of vertices
    */
    setVertex(v,i,j)
    {
        //Your code here
        var vertex_id = 3 * (i*(this.div+1) + j);   //populate vBuffers
        this.vBuffer[vertedx_id] = v[0];
        this.vBuffer[vertedx_id + 1] = v[1];
        this.vBuffer[vertedx_id + 2] = v[2];
    }
    
    /**
    * Return the x,y,z coordinates of a vertex at location (i,j)
    * @param {Object} v an an array of length 3 holding x,y,z coordinates
    * @param {number} i the ith row of vertices
    * @param {number} j the jth column of vertices
    */
    getVertex(v,i,j)
    {
        //Your code here
        var vertex_id = 3 * (i*(this.div+1) + j);
        v[0] = this.vBuffer[vertex_id];
        v[1] = this.vBuffer[vertex_id + 1];
        v[2] = this.vBuffer[vertex_id + 2];

    }
    
    /**
    * Send the buffer objects to WebGL for rendering 
    */
    loadBuffers()
    {
        // Specify the vertex coordinates
        this.VertexPositionBuffer = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, this.VertexPositionBuffer);      
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(this.vBuffer), gl.STATIC_DRAW);
        this.VertexPositionBuffer.itemSize = 3;
        this.VertexPositionBuffer.numItems = this.numVertices;
        console.log("Loaded ", this.VertexPositionBuffer.numItems, " vertices");
    
        // Specify normals to be able to do lighting calculations
        this.VertexNormalBuffer = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, this.VertexNormalBuffer);
        gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(this.nBuffer),
                  gl.STATIC_DRAW);
        this.VertexNormalBuffer.itemSize = 3;
        this.VertexNormalBuffer.numItems = this.numVertices;
        console.log("Loaded ", this.VertexNormalBuffer.numItems, " normals");
    
        // Specify faces of the terrain 
        this.IndexTriBuffer = gl.createBuffer();
        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.IndexTriBuffer);
        gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(this.fBuffer),
                  gl.STATIC_DRAW);
        this.IndexTriBuffer.itemSize = 1;
        this.IndexTriBuffer.numItems = this.fBuffer.length;
        console.log("Loaded ", this.IndexTriBuffer.numItems, " triangles");
    
        //Setup Edges  
        this.IndexEdgeBuffer = gl.createBuffer();
        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.IndexEdgeBuffer);
        gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint32Array(this.eBuffer),
                  gl.STATIC_DRAW);
        this.IndexEdgeBuffer.itemSize = 1;
        this.IndexEdgeBuffer.numItems = this.eBuffer.length;
        
        console.log("triangulatedPlane: loadBuffers");
    }
    
    /**
    * Render the triangles 
    */
    drawTriangles(){
        gl.bindBuffer(gl.ARRAY_BUFFER, this.VertexPositionBuffer);
        gl.vertexAttribPointer(shaderProgram.vertexPositionAttribute, this.VertexPositionBuffer.itemSize, 
                         gl.FLOAT, false, 0, 0);

        // Bind normal buffer
        gl.bindBuffer(gl.ARRAY_BUFFER, this.VertexNormalBuffer);
        gl.vertexAttribPointer(shaderProgram.vertexNormalAttribute, 
                           this.VertexNormalBuffer.itemSize,
                           gl.FLOAT, false, 0, 0);   
    
        //Draw 
        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.IndexTriBuffer);
        gl.drawElements(gl.TRIANGLES, this.IndexTriBuffer.numItems, gl.UNSIGNED_INT,0);
    }
    
    /**
    * Render the triangle edges wireframe style 
    */
    drawEdges(){
    
        gl.bindBuffer(gl.ARRAY_BUFFER, this.VertexPositionBuffer);
        gl.vertexAttribPointer(shaderProgram.vertexPositionAttribute, this.VertexPositionBuffer.itemSize, 
                         gl.FLOAT, false, 0, 0);

        // Bind normal buffer
        gl.bindBuffer(gl.ARRAY_BUFFER, this.VertexNormalBuffer);
        gl.vertexAttribPointer(shaderProgram.vertexNormalAttribute, 
                           this.VertexNormalBuffer.itemSize,
                           gl.FLOAT, false, 0, 0);   
    
        //Draw 
        gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.IndexEdgeBuffer);
        gl.drawElements(gl.LINES, this.IndexEdgeBuffer.numItems, gl.UNSIGNED_INT,0);   
    }
/**
 * Fill the vertex and buffer arrays 
 * Inputs: none
 * Outputs: Terrain map generated, populated buffers.
 */    
generateTriangles()
{

    //Your code here
    var delta_x = (this.maxX-this.minX)/this.div;
    var delta_y = (this.maxY-this.minY)/this.div;

    for(var i=0; i<=this.div; i++){
        for(var j=0; j<=this.div; j++){
            this.vBuffer.push(this.minX+delta_x*j);             //vertices (x,y,z)
            this.vBuffer.push(this.minY+delta_y*i);
            this.vBuffer.push(0);         //flat generation for value of z

            this.nBuffer.push(0);       //normal vector
            this.nBuffer.push(0);
            this.nBuffer.push(0);
        }
    }
    for(var i=0; i<this.div; i++){          //V_id by point in an array
        for(var j=0; j<this.div; j++){
            var vid = i*(this.div+1) + j;       //row-major order calculation for the triangle-division
            this.fBuffer.push(vid);             //trangle FACE buffer
            this.fBuffer.push(vid+1);
            this.fBuffer.push(vid+this.div+1);

            this.fBuffer.push(vid + 1);         //next triange-face neighbor
            this.fBuffer.push(vid+1+this.div+1);
            this.fBuffer.push(vid+this.div+1);
            
        }
    }
    
    console.log("vBUffer", this.vBuffer.length);
    this.numVertices = this.vBuffer.length/3;
    this.numFaces = this.fBuffer.length/3;
}

/**
 * Make z generation for terrain, raise the elevation using random generation of nromal vectors
 * Inputs: none
 * Outputs: raised z-coordinate, modeled terrain
 */
z_generation(){
    
    var delta_change = 0.0055;
    for(var k=0; k<175; k++){
    var normal_x = Math.random() - this.maxX;      //generate a random x-coordinate in the terrain map
    var normal_y = Math.random() - this.maxY;      //generate a random y-coordinate in the terrain map
    var two_x =  Math.random() - this.maxX;      //generate a random x-coordinate in the terrain map
    var two_y =  Math.random() - this.maxY;      //generate a random y-coordinate in the terrain map
    // var normal = glMatrix.vec3.fromValues(normal_x-two_x, normal_y-two_y, 0);  
    let last_normal = glMatrix.vec3.create();
    last_normal = glMatrix.vec3.fromValues(normal_x, normal_y, 0);

    for(var i=2; i<=(this.div)*(this.div)*4; i+=3){
            var varx = this.vBuffer[i-2] - two_x;
            var vary = this.vBuffer[i-1] - two_y;
            var normal = Math.sqrt(Math.pow(varx, 2) + Math.pow(vary, 2));      //normalize the vector
            varx = varx/normal;
            vary = vary/normal;
            let second_vec = glMatrix.vec3.create();
            second_vec = glMatrix.vec3.fromValues(varx, vary, 0);

            let result = glMatrix.vec3.create();
            result = glMatrix.vec3.dot(second_vec, last_normal);
            if(result < 0){                         //subtract delta_change
                this.vBuffer[i] -= delta_change;
            } else{                                    //add delta_change
                this.vBuffer[i] += delta_change;
            }
    }
}
}

/*calculate the average normals for each plane in the terrain
Inputs: none
Outputs: shaded terrain*/ 
avg_normals(){
    for(var j=0; j<this.nBuffer.length; j++){
        this.count.push(0);
    }
    for(var i=0; i<this.fBuffer.length; i+=3){
        var face1 = this.fBuffer[i]*3;      //get vertex indices from fBuffer
        var face2 = this.fBuffer[i+1]*3;
        var face3 = this.fBuffer[i+2]*3;
/**************************************start creating normals, populating buffers***************************************** */
        var x_one = this.vBuffer[face1];         //get vertex information from vBuffer, corresponding to faces
        var y_one = this.vBuffer[face1+1];      //set one
        var z_one = this.vBuffer[face1+2];

        var x_two = this.vBuffer[face2];        //set two
        var y_two = this.vBuffer[face2+1];
        var z_two = this.vBuffer[face2+2];

        var x_three = this.vBuffer[face3];      //set three
        var y_three = this.vBuffer[face3+1];
        var z_three = this.vBuffer[face3+2];

        var vertex_one = [x_two - x_one, y_two - y_one, z_two - z_one];         //calculate vector for cross products
        var vertex_two = [x_three - x_one, y_three - y_one, z_three - z_one];  
        var norm = glMatrix.vec3.create();
        glMatrix.vec3.cross(norm, vertex_one, vertex_two)

        this.nBuffer[face1] += norm[0];     //add the normalized components to all 3 vertices
        this.nBuffer[face1+1] += norm[1];
        this.nBuffer[face1+2] += norm[2];
        this.nBuffer[face2] += norm[0];     //vertex 2
        this.nBuffer[face2+1] += norm[1];
        this.nBuffer[face2+2] += norm[2];
        this.nBuffer[face3] += norm[0];     //vertex 3
        this.nBuffer[face3+1] += norm[1];
        this.nBuffer[face3+2] += norm[2];
        
        
        this.count[face1] +=1;      //count for vertex 1 
        this.count[face1+1] +=1;
        this.count[face1+2] +=1;
        this.count[face2] +=1;      //count for vertex 2
        this.count[face2+1] +=1;
        this.count[face2+2] +=1;
        this.count[face3] +=1;      //count for vertex 3
        this.count[face3+1] +=1;
        this.count[face3+2] +=1;
    }
    for(var k=0; k<this.nBuffer.length; k+= 3){                   //take the average of the values at the end 
        // this.nBuffer[k] = this.nBuffer[k] / this.count[k];
        let norm = [this.nBuffer[k], this.nBuffer[k+1], this.nBuffer[k+2]];
        glMatrix.vec3.normalize(norm, norm);
        this.nBuffer[k] = norm[0];
        this.nBuffer[k+1] = norm[1];
        this.nBuffer[k+2] = norm[2];
    }
}

/**
 * Print vertices and triangles to console for debugging
 */
printBuffers()
    {
        
    for(var i=0;i<this.numVertices;i++)
          {
           console.log("v ", this.vBuffer[i*3], " ", 
                             this.vBuffer[i*3 + 1], " ",
                             this.vBuffer[i*3 + 2], " ");
                       
          }
    
      for(var i=0;i<this.numFaces;i++)
          {
           console.log("f ", this.fBuffer[i*3], " ", 
                             this.fBuffer[i*3 + 1], " ",
                             this.fBuffer[i*3 + 2], " ");
                       
          }
        
    }

/**
 * Generates line values from faces in faceArray
 * to enable wireframe rendering
 */
generateLines()
{
    var numTris=this.fBuffer.length/3;
    for(var f=0;f<numTris;f++)
    {
        var fid=f*3;
        this.eBuffer.push(this.fBuffer[fid]);
        this.eBuffer.push(this.fBuffer[fid+1]);
        
        this.eBuffer.push(this.fBuffer[fid+1]);
        this.eBuffer.push(this.fBuffer[fid+2]);
        
        this.eBuffer.push(this.fBuffer[fid+2]);
        this.eBuffer.push(this.fBuffer[fid]);
    }
    
}
    
}
