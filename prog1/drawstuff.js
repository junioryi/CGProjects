/* classes */ 

// Color class
class Color {
    constructor(r,g,b,a) {
        try {
            if ((typeof(r) !== "number") || (typeof(g) !== "number") || (typeof(b) !== "number") || (typeof(a) !== "number"))
                throw "color component not a number";
            else if ((r<0) || (g<0) || (b<0) || (a<0)) 
                throw "color component less than 0";
            else if ((r>255) || (g>255) || (b>255) || (a>255)) 
                throw "color component bigger than 255";
            else {
                this[0] = r; this[1] = g; this[2] = b; this[3] = a; 
            }
        } // end try
        
        catch (e) {
            console.log(e);
        }
    } // end Color constructor

        // Color change method
    change(r,g,b,a) {
        try {
            if ((typeof(r) !== "number") || (typeof(g) !== "number") || (typeof(b) !== "number") || (typeof(a) !== "number"))
                throw "color component not a number";
            else if ((r<0) || (g<0) || (b<0) || (a<0)) 
                throw "color component less than 0";
            else if ((r>255) || (g>255) || (b>255) || (a>255)) 
                throw "color component bigger than 255";
            else {
                this[0] = r; this[1] = g; this[2] = b; this[3] = a; 
            }
        } // end try
        
        catch (e) {
            console.log(e);
        }
    } // end Color change method 

    add(c) {
        try {
            if (!(c instanceof Color)) 
                throw "Adding non-color to color.";
            else {
                for (var i=0; i<3; ++i) 
                    this[i] += c[i];
            }
        }
        catch(e) {
            console.log(e);
        }
    }

    toConsole(prefix) {
        console.log(prefix+"["+this[0]+","+this[1]+","+this[2]+"]");
    } // end to console

} // end color class

// Vector class
class Vector { 
    constructor(x,y,z) {
        this.set(x,y,z);
    } // end constructor
    
    // sets the components of a vector
    set(x,y,z) {
        try {
            if ((typeof(x) !== "number") || (typeof(y) !== "number") || (typeof(z) !== "number"))
                throw "vector component not a number";
            else
                this.x = x; this.y = y; this.z = z; 
        } // end try
        
        catch(e) {
            console.log(e);
        }
    } // end vector set
    
    // copy the passed vector into this one
    copy(v) {
        try {
            if (!(v instanceof Vector))
                throw "Vector.copy: non-vector parameter";
            else
                this.x = v.x; this.y = v.y; this.z = v.z;
        } // end try
        
        catch(e) {
            console.log(e);
        }
    }
    
    toConsole(prefix) {
        console.log(prefix+"["+this.x+","+this.y+","+this.z+"]");
    } // end to console
    
    // static dot method
    static dot(v1,v2) {
        try {
            if (!(v1 instanceof Vector) || !(v2 instanceof Vector))
                throw "Vector.dot: non-vector parameter";
            else
                return(v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
        } // end try
        
        catch(e) {
            console.log(e);
            return(NaN);
        }
    } // end dot static method
    
    // static add method
    static add(v1,v2) {
        try {
            if (!(v1 instanceof Vector) || !(v2 instanceof Vector))
                throw "Vector.add: non-vector parameter";
            else
                return(new Vector(v1.x+v2.x,v1.y+v2.y,v1.z+v2.z));
        } // end try
        
        catch(e) {
            console.log(e);
            return(new Vector(NaN,NaN,NaN));
        }
    } // end add static method

    static add3(v1, v2, v3) {
        try {
            if (!(v1 instanceof Vector) || !(v2 instanceof Vector)
                    || !(v3 instanceof Vector)) {
                throw "Vector.add: non-vector parameter";
            }
            else {
                return(new Vector(
                            v1.x + v2.x + v3.x,
                            v1.y + v2.y + v3.y,
                            v1.z + v2.z + v3.z));
            }
        } // end try
        catch(e) {
            console.log(e);
            return(new Vector(NaN, NaN, NaN));
        }
    } // end add3

    // static subtract method, v1-v2
    static subtract(v1,v2) {
        try {
            if (!(v1 instanceof Vector) || !(v2 instanceof Vector))
                throw "Vector.subtract: non-vector parameter";
            else {
                var v = new Vector(v1.x-v2.x,v1.y-v2.y,v1.z-v2.z);
                //v.toConsole("Vector.subtract: ");
                return(v);
            }
        } // end try
        
        catch(e) {
            console.log(e);
            return(new Vector(NaN,NaN,NaN));
        }
    } // end subtract static method

    // static scale method
    static scale(c,v) {
        try {
            if (!(typeof(c) === "number") || !(v instanceof Vector))
                throw "Vector.scale: malformed parameter";
            else
                return(new Vector(c*v.x,c*v.y,c*v.z));
        } // end try
        
        catch(e) {
            console.log(e);
            return(new Vector(NaN,NaN,NaN));
        }
    } // end scale static method
    
    // static normalize method
    static normalize(v) {
        try {
            if (!(v instanceof Vector))
                throw "Vector.normalize: parameter not a vector";
            else {
                var lenDenom = 1/Math.sqrt(Vector.dot(v,v));
                return(Vector.scale(lenDenom,v));
            }
        } // end try
        
        catch(e) {
            console.log(e);
            return(new Vector(NaN,NaN,NaN));
        }
    } // end scale static method

    // static cross method
    static cross(u, v) {
        try {
            if (!(u instanceof Vector) || !(v instanceof Vector))
                throw "Vector.cross: parameter not a vector.";
            else {
                var temp = new Vector(
                        u.y*v.z - u.z*v.y,
                        u.z*v.x - u.x*v.z,
                        u.x*v.y - u.y*v.x);
                return(Vector.normalize(temp));
            }
        } // end try
        catch(e) {
            console.log(e);
            return(new Vector(NaN, NaN, NaN));
        }
    } // end cross static method
    
} // end Vector class

class Triangle {
    constructor(json) {
        //console.log(typeof(json));
        this.set(json);
    }
    set(json) {

        // TODO:Maybe valid input here
        
        this.v0 = new Vector(json.v0[0], json.v0[1], json.v0[2]);
        this.v1 = new Vector(json.v1[0], json.v1[1], json.v1[2]);
        this.v2 = new Vector(json.v2[0], json.v2[1], json.v2[2]);

        this.u = new Vector(0,0,0);
        this.v = new Vector(0,0,0);

        this.u = Vector.subtract(this.v1, this.v0);
        this.v = Vector.subtract(this.v2, this.v0);

        this.n = Vector.normalize(Vector.cross(this.u, this.v));

        this.uu = Vector.dot(this.u, this.u);
        this.vv = Vector.dot(this.v, this.v);
        this.uv = Vector.dot(this.u, this.v);

        this.uv_square = this.uv * this.uv;

        this.ambient = new Array(3);
        this.diffuse = new Array(3);
        this.specular = new Array(3);
        for (var i=0; i<3; ++i) {
            this.ambient[i] = json.ambient[i];
            this.diffuse[i] = json.diffuse[i];
            this.specular[i] = json.specular[i];
        }

        this.n_spec = json.n;
    }
} // end Triangle class

/* utility functions */

// generate n integers in random order
// uses Fisher-Yates shuffle
function randPermutation(n) {
    var array = new Array(n);
    var bagSize = n, temp, randChoice;
    
    // fill the array 
    for (var i=0; i<n; i++)
        array[i] = i; 

    // while the bag isn't empty, pick from it
    while (bagSize !== 0) {
        randChoice = Math.floor(Math.random() * bagSize); // pick from bag
        bagSize--; // bag is less one
        temp = array[bagSize]; // remember what was at new bag slot
        array[bagSize] = array[randChoice]; // move old pick to new slot
        array[randChoice] = temp; // copy old element to old slot
    } // end while

    // for (i=0; i<n; i++)
    //    console.log(array[i]);

    return(array);
}

// get the JSON file from the passed URL
function getJSONFile(url,descr) {
    try {
        if ((typeof(url) !== "string") || (typeof(descr) !== "string"))
            throw "getJSONFile: parameter not a string";
        else {
            var httpReq = new XMLHttpRequest(); // a new http request
            httpReq.open("GET",url,false); // init the request
            httpReq.send(null); // send the request
            var startTime = Date.now();
            while ((httpReq.status !== 200) && (httpReq.readyState !== XMLHttpRequest.DONE)) {
                if ((Date.now()-startTime) > 3000)
                    break;
            } // until its loaded or we time out after three seconds
            if ((httpReq.status !== 200) || (httpReq.readyState !== XMLHttpRequest.DONE))
                throw "Unable to open "+descr+" file!";
            else
                return JSON.parse(httpReq.response); 
        } // end if good params
    } // end try    
    
    catch(e) {
        console.log(e);
        return(String.null);
    }
} // end get input spheres

// Solve quadratic. Return empty array if no solutions, 
// one t value if one solution, two if two solutions.
function solveQuad(a,b,c) {
    var discr = b*b - 4*a*c; 
    // console.log("a:"+a+" b:"+b+" c:"+c);

    if (discr < 0) { // no solutions
        // console.log("no roots!");
        return([]); 
    } else if (discr == 0) { // one solution
        // console.log("root: "+(-b/(2*a)));
        return([-b/(2*a)]);
    } else { // two solutions
        var denom = 0.5/a;
        var term1 = -b;
        var term2 = Math.sqrt(discr)
        var tp = denom * (term1 + term2);
        var tm = denom * (term1 - term2);
        // console.log("root1:"+tp+" root2:"+tm);
        if (tm < tp)
            return([tm,tp]);
        else
            return([tp,tm]);
    } 
} // end solveQuad
    
// ray sphere intersection
// if no intersect, return NaN
// if intersect, return xyz vector and t value
// intersects in front of clipVal don't count
//
// Return example: {"exist": true, "xyz": isect, "t": qsolve[0]}
function raySphereIntersect(ray,sphere,clipVal) {
    try {
        if (!(ray instanceof Array) || !(sphere instanceof Object))
            throw "RaySphereIntersect: ray or sphere are not formatted well";
        else if (ray.length != 2)
            throw "RaySphereIntersect: badly formatted ray";
        else { // valid params
            // ray = [eye, Dir]
            var a = Vector.dot(ray[1],ray[1]); // dot(D,D)
            var origMctr = Vector.subtract(ray[0],new Vector(sphere.x,sphere.y,sphere.z)); // E-C
            var b = 2 * Vector.dot(ray[1],origMctr); // 2 * dot(D,E-C)
            var c = Vector.dot(origMctr,origMctr) - sphere.r*sphere.r; // dot(E-C,E-C) - r^2
            // if (clipVal == 0) {
            //     ray[0].toConsole("ray.orig: ");
            //     ray[1].toConsole("ray.dir: ");
            //     console.log("a:"+a+" b:"+b+" c:"+c);
            // } // end debug case
        
            var qsolve = solveQuad(a,b,c);
            if (qsolve.length == 0) 
                throw "no intersection";
            else if (qsolve.length == 1) { 
                if (qsolve[0] < clipVal)
                    throw "intersection too close";
                else {
                    // Eye + D*t
                    var isect = Vector.add(ray[0],Vector.scale(qsolve[0],ray[1]));
                    //console.log("t: "+qsolve[0]);
                    //isect.toConsole("intersection: ");
                    return({"exists": true, "xyz": isect,"t": qsolve[0]});  
                } // one unclipped intersection
            } else if (qsolve[0] < clipVal) { // Two intersection
                if (qsolve[1] < clipVal)
                    throw "intersections too close";
                else { 
                    var isect = Vector.add(ray[0],Vector.scale(qsolve[1],ray[1]));
                    //console.log("t2: "+qsolve[1]);
                    //isect.toConsole("intersection: ");
                    return({"exists": true, "xyz": isect,"t": qsolve[1]});  
                } // one intersect too close, one okay
            } else {
                var isect = Vector.add(ray[0],Vector.scale(qsolve[0],ray[1]));
                //console.log("t1: "+qsolve[0]);
                //isect.toConsole("intersection: ");
                return({"exists": true, "xyz": isect,"t": qsolve[0]});  
            } // both not too close
        } // end if valid params
    } // end try

    catch(e) {
        //console.log(e);
        return({"exists": false, "xyz": NaN, "t": NaN});
    }
} // end raySphereIntersect

// Ray Triangle intersection
function rayTriangleIntersect(ray, triangle, clipVal) {
    try {
        if (!(ray instanceof Array) || !(triangle instanceof Object))
            throw "RayTriangleIntersect: ray or triangle are not formatted well.";
        else if (ray.length != 2)
            throw "RayTriangleIntersect: badly formatted ray.";
        else { // valid params
            // ray = [eye, Dir]
            //console.log("triangle norm: " + triangle.n[0] + " " + triangle.n[1] + " " + triangle[2]);
            var w = Vector.subtract(triangle.v0, ray[0]); // v0 - Eye
            var t = Vector.dot(triangle.n, w) / Vector.dot(triangle.n, ray[1]);
            if (t > 0) {
                var hit = Vector.add(ray[0], Vector.scale(t, ray[1]));
                var toHit = Vector.subtract(hit, triangle.v0);
                var dot00 = triangle.uu; // ac is u
                var dot01 = triangle.uv; 
                var dot02 = Vector.dot(triangle.u, toHit);
                var dot11 = triangle.vv; // ab is v
                var dot12 = Vector.dot(triangle.v, toHit);
                var divide = dot00 * dot11 - dot01 * dot01;
                var u = (dot11 * dot02 - dot01 * dot12) /divide;
                var v = (dot00 * dot12 - dot01 * dot02) /divide;
                if (u >= 0 && v >= 0 && u+v <= 1) 
                    return({"exists": true, "xyz": hit, "t": t});
            }
            return({"exists": false, "xyz":NaN, "t":NaN});
        }
    }
    catch(e) {
        console.log(e);
        return({"exists":false, "xyz": NaN, "t":NaN});
    }
} // end ray triangle intersection

// draw a pixel at x,y using color
function drawPixel(imagedata,x,y,color) {
    try {
        if ((typeof(x) !== "number") || (typeof(y) !== "number"))
            throw "drawpixel location not a number";
        else if ((x<0) || (y<0) || (x>=imagedata.width) || (y>=imagedata.height))
            throw "drawpixel location outside of image";
        else if (color instanceof Color) {
            var pixelindex = (y*imagedata.width + x) * 4;
            imagedata.data[pixelindex] = color[0];
            imagedata.data[pixelindex+1] = color[1];
            imagedata.data[pixelindex+2] = color[2];
            imagedata.data[pixelindex+3] = color[3];
        } else 
            throw "drawpixel color is not a Color";
    } // end try
    
    catch(e) {
        console.log(e);
    }
} // end drawPixel
    
// draw random pixels
function drawRandPixels(context) {
    var c = new Color(0,0,0,0); // the color at the pixel: black
    var w = context.canvas.width;
    var h = context.canvas.height;
    var imagedata = context.createImageData(w,h);
    const PIXEL_DENSITY = 0.01;
    var numPixels = (w*h)*PIXEL_DENSITY; 
    
    // Loop over 1% of the pixels in the image
    for (var x=0; x<numPixels; x++) {
        c.change(Math.random()*255,Math.random()*255,
            Math.random()*255,255); // rand color
        drawPixel(imagedata,
            Math.floor(Math.random()*w),
            Math.floor(Math.random()*h),
                c);
    } // end for x
    context.putImageData(imagedata, 0, 0);
} // end draw random pixels

// put random points in the spheres from the class github
function drawRandPixelsInInputSpheres(context) {
    var inputSpheres = getJSONFile(INPUT_SPHERES_URL,"spheres");
    var w = context.canvas.width; //512
    var h = context.canvas.height;
    var imagedata = context.createImageData(w,h);
    const PIXEL_DENSITY = 0.1;
    var numCanvasPixels = (w*h)*PIXEL_DENSITY; 
    
    if (inputSpheres != String.null) { 
        var x = 0; var y = 0; // pixel coord init
        var cx = 0; var cy = 0; // init center x and y coord
        var sphereRadius = 0; // init sphere radius
        var numSpherePixels = 0; // init num pixels in sphere
        var c = new Color(0,0,0,0); // init the sphere color
        var n = inputSpheres.length;
        //console.log("number of spheres: " + n);

        // Loop over the spheres, draw rand pixels in each
        for (var s=0; s<n; s++) {
            //console.log("s: "+s);
            cx = w*inputSpheres[s].x; // sphere center x
            cy = h*inputSpheres[s].y; // sphere center y
            sphereRadius = Math.round(w*inputSpheres[s].r); // radius
            numSpherePixels = sphereRadius*4*Math.PI; // sphere area
            numSpherePixels *= PIXEL_DENSITY; // percentage of sphere on
            numSpherePixels = Math.round(numSpherePixels);
            //console.log("sphere radius: "+sphereRadius);
            //console.log("num sphere pixels: "+numSpherePixels);
            c.change(
                inputSpheres[s].diffuse[0]*255,
                inputSpheres[s].diffuse[1]*255,
                inputSpheres[s].diffuse[2]*255,
                255); // rand color
            for (var p=0; p<numSpherePixels; p++) {
                do {
                    x = Math.random()*2 - 1; // in unit square 
                    y = Math.random()*2 - 1; // in unit square
                } while (Math.sqrt(x*x + y*y) > 1)
                var tempx = cx+Math.round(x*sphereRadius);
                //console.log("x: "+x+", cx: "+tempx);
                drawPixel(imagedata,
                    cx+Math.floor(x*sphereRadius),
                    cy+Math.floor(y*sphereRadius),c);
                //console.log("color: ("+c[0]+","+c[1]+","+c[2]+")");
                //console.log("x:" + tempx);
                //console.log("x: "+Math.round(w*inputSpheres[s].x));
                //console.log("y: "+Math.round(h*inputSpheres[s].y));
            } // end for pixels in sphere
        } // end for spheres
        context.putImagedata(imagedata, 0, 0);
    } // end if spheres found
} // end draw rand pixels in input spheres

// draw 2d projections read from the JSON file at class github
function drawInputSpheresUsingArcs(context) {
    var inputSpheres = getJSONFile(INPUT_SPHERES_URL,"spheres");
    
    if (inputSpheres != String.null) { 
        var c = new Color(0,0,0,0); // the color at the pixel: black
        var w = context.canvas.width;
        var h = context.canvas.height;
        var n = inputSpheres.length; 
        //console.log("number of spheres: " + n);

        // Loop over the spheres, draw each in 2d
        for (var s=0; s<n; s++) {
            context.fillStyle = 
                "rgb(" + Math.floor(inputSpheres[s].diffuse[0]*255)
                +","+ Math.floor(inputSpheres[s].diffuse[1]*255)
                +","+ Math.floor(inputSpheres[s].diffuse[2]*255) +")"; // rand color
            context.beginPath();
            context.arc(
                Math.round(w*inputSpheres[s].x),
                Math.round(h*inputSpheres[s].y),
                Math.round(w*inputSpheres[s].r),
                0,2*Math.PI);
            context.fill();
            //console.log(context.fillStyle);
            //console.log("x: "+Math.round(w*inputSpheres[s].x));
            //console.log("y: "+Math.round(h*inputSpheres[s].y));
            //console.log("r: "+Math.round(w*inputSpheres[s].r));
        } // end for spheres
    } // end if spheres found
} // end draw input spheres

function isLightOccluded(L,isectPos,isectSphere,spheres) {
    var s=0; // which sphere
    var lightOccluded = false; // if light is occluded
    var occluderIsect = {}; // occluder intersect details
    
    var numSpheres = spheres.length;
    while ((!lightOccluded) && (s<numSpheres)) {
        if (s == isectSphere) {
            // Skip the sphere itself.
            s++;
            continue;
        }
        occluderIsect = raySphereIntersect([isectPos, L], spheres[s], 0);
        if (!occluderIsect.exists) { // no intersection
            s++;
        } else if (occluderIsect.t > 1) { // light in front of intersection
            s++;
        } else {
            lightOccluded = true;
        }
    }
     
    return(lightOccluded);
} // end is light occluded

// 
function shadeIsectTriangle(isect, triangle, lights, spheres, addAmbient) {
    try {
        if (!(isect instanceof Object) || !(triangle instanceof Triangle) 
            || !(lights instanceof Array)) {
            throw "shadeIsectTriangle: bad parameter passed.";
        }   
        
        var c = new Color(0,0,0,255);

        var lightOccluded = false;
        var Lloc = new Vector(0,0,0);
        for (var l=0; l<lights.length; l++) {
            // Add ambient colors
            for (var i=0; i<3; i++) {
                c[i] += lights[l].ambient[i] * triangle.ambient[i];
            }
            Lloc.set(lights[l].x, lights[l].y, lights[l].z);
            var L = Vector.subtract(Lloc, isect.xyz);
            //L.toConsole("Light path: ");

            // if light isn't occluded
            if (!isLightOccluded(L, isect.xyz, -1, spheres)) { // -1 means test all spheres

                // normal is plane normal
                var N = triangle.n;
                var diffFactor = Math.max(0, Math.abs(Vector.dot(N, Vector.normalize(L))));
                if (diffFactor > 0) {
                    for (var i=0; i<3; i++) {
                        c[i] += lights[l].diffuse[i] * triangle.diffuse[i] * diffFactor;
                    }
                } // end nonzero diffuse factor

                // add in the specular light
                var V = Vector.normalize(Vector.subtract(Eye, isect.xyz)); // View vector
                var H = Vector.normalize(Vector.add(L, V)); // half vector

                var specFactor = Math.max(0, Math.abs(Vector.dot(N, H)));
                if (specFactor > 0) {
                    var newSpecFactor = specFactor;
                    for (var s=1; s<triangle.n_spec; s++) 
                        newSpecFactor *= specFactor;
                    for (var i=0; i<3; i++) 
                        c[i] += lights[l].specular[i] * triangle.specular[i] * newSpecFactor;
                }
            } // end if light is not occluded
        } // done all light
        for (var i=0; i<3; i++) 
            // 1 is black
            c[i] = 255*Math.min(1, c[i]);
        return(c);
    }
    catch(e) {
        console.log("shadeIsectTriangle error: "+e);
        return(Object.null);
    }
}

// color the passed intersection and sphere
function shadeIsect(isect,isectSphere,lights,spheres, addAmbient) {
    try {
        if (   !(isect instanceof Object) || !(typeof(isectSphere) === "number") 
            || !(lights instanceof Array) || !(spheres instanceof Array))
            throw "shadeIsect: bad parameter passed";
        else {
            var c = new Color(0,0,0,255); // init the sphere color to black
            var sphere = spheres[isectSphere]; // sphere intersected by eye
            // console.log("shading pixel");

            // add light for each source
            var lightOccluded = false; // if an occluder is found
            var Lloc = new Vector(0,0,0);
            for (var l=0; l<lights.length; l++) {

                // add in the ambient light
                c[0] += lights[l].ambient[0] * sphere.ambient[0]; // ambient term r
                c[1] += lights[l].ambient[1] * sphere.ambient[1]; // ambient term g
                c[2] += lights[l].ambient[2] * sphere.ambient[2]; // ambient term b
                
                // check each other sphere to see if it occludes light
                Lloc.set(lights[l].x,lights[l].y,lights[l].z);
                var L = Vector.subtract(Lloc,isect.xyz); // light vector unnorm'd
                // L.toConsole("L: ");
                // console.log("isect: "+isect.xyz.x+", "+isect.xyz.y+", "+isect.xyz.z);
                
                // if light isn't occluded
                if (!isLightOccluded(L,isect.xyz,isectSphere,spheres)) {
                    // console.log("no occlusion found");
                    
                    // add in the diffuse light
                    var sphereCenter = new Vector(sphere.x,sphere.y,sphere.z);
                    var N = Vector.normalize(Vector.subtract(isect.xyz,sphereCenter)); // surface normal
                    var diffFactor = Math.max(0,Vector.dot(N,Vector.normalize(L)));
                    if (diffFactor > 0) {
                        c[0] += lights[l].diffuse[0] * sphere.diffuse[0] * diffFactor;
                        c[1] += lights[l].diffuse[1] * sphere.diffuse[1] * diffFactor;
                        c[2] += lights[l].diffuse[2] * sphere.diffuse[2] * diffFactor;
                    } // end nonzero diffuse factor

                    // add in the specular light
                    var V = Vector.normalize(Vector.subtract(Eye,isect.xyz)); // view vector
                    // TODO:
                    // L need to be normalize?
                    var H = Vector.normalize(Vector.add(L,V)); // half vector
                    var specFactor = Math.max(0,Vector.dot(N,H)); 
                    if (specFactor > 0) {
                        var newSpecFactor = specFactor;
                        for (var s=1; s<spheres[isectSphere].n; s++) // mult by itself if needed
                            newSpecFactor *= specFactor;
                        c[0] += lights[l].specular[0] * sphere.specular[0] * newSpecFactor; // specular term
                        c[1] += lights[l].specular[1] * sphere.specular[1] * newSpecFactor; // specular term
                        c[2] += lights[l].specular[2] * sphere.specular[2] * newSpecFactor; // specular term
                    } // end nonzero specular factor
                    
                } // end if light not occluded
            } // end for lights
            
            c[0] = 255 * Math.min(1,c[0]); // clamp max value to 1
            c[1] = 255 * Math.min(1,c[1]); // clamp max value to 1
            c[2] = 255 * Math.min(1,c[2]); // clamp max value to 1

            return(c);
        } // if have good params
    } // end throw
    
    catch(e) {
        console.log(e);
        return(Object.null);
    }
}

function hitCheck(viewpoint, Dir, inputSpheres, triangles, inputLights, frontWall) {
    try {
        var closestT = Number.MAX_VALUE;
        var c = new Color(0,0,0,255);
        var numSpheres = inputSpheres.length;
        var numTriangles = triangles.length;
        // Ignore front wall
        if (!(frontWall)) 
            numTriangles -= 2;
        //Dir.copy(Vector.subtract(new Vector(wx, wy, WIN_Z), Eye));
        // Keep the intersection information, compute color afterward
        var hitSphere = false;
        var tempIsect = {"exists":false};
        var sphereIdx = -1;
        var tempTriangle;
        // Check spheres
        for (var s=0; s<numSpheres; s++) {
            isect = raySphereIntersect([viewpoint,Dir],inputSpheres[s],1); 
            if (isect.exists) {// there is an intersect
                if (isect.t < closestT) { // it is the closest yet
                    hitSphere = true;
                    closestT = isect.t; // record closest t yet
                    //c = shadeIsect(isect,s,inputLights,inputSpheres); 
                    tempIsect = isect;
                    sphereIdx = s;
                } // end if closest yet
            }
        } // done all spheres
        // Check trangles
        for (var t=0; t<numTriangles; t++) { // loop triangles
            isect = rayTriangleIntersect([viewpoint, Dir], triangles[t], 1);
            if (isect.exists) {
                if (isect.t < closestT) {
                    hitSphere = false;
                    closestT = isect.t;
                    tempIsect = isect;
                    tempTriangle = triangles[t];
                }
                // Triangle can break since it cannot hit two triangles in this case.
                break;
            }
        } // done all triangles

        // Compute color 
        if (!(tempIsect.exists)) {
            throw "tracer not hit anything";
        }
        return ({"isect": tempIsect, "hitSphere": hitSphere, "sphereIdx":sphereIdx, "triangle":tempTriangle});

    }
    catch(e) {
        //console.log(e);
        return ({"isect": {"exists":false}, "hitSphere": false, "sphereIdx":-1, "triangle":null});
    }
}

function radiance(Dir, tempIsect, hitSphere, sphereIdx, tempTriangle, triangles, inputLights, inputSpheres, russian, numSample, addAmbient) {
    try{
        if (russian >= 1.001) 
            throw "Russian roulette should be less than 1.";
        var c = new Color(0,0,0,0);
        if (hitSphere) {
            // Direct
            c = shadeIsect(tempIsect, sphereIdx, inputLights, inputSpheres, addAmbient);
            // Russian roulette
            if (Math.random() < russian) {
                var tempColor = new Color(0,0,0,0);
                // Sampling
                for (var s=0; s<numSample; ++s) {
                    // Find the intersection surface
                    var sphereCenter = new Vector(
                            inputSpheres[sphereIdx].x,
                            inputSpheres[sphereIdx].y,
                            inputSpheres[sphereIdx].z);
                    // Finding orthogonal vector on intersection surface
                    var surface_n = Vector.normalize(
                            Vector.subtract(tempIsect.xyz, sphereCenter)); 
                    var surface_x = Vector.cross(Vector.normalize(Dir), surface_n);
                    var surface_y = Vector.cross(surface_x, surface_n);

                    // Random sampling over hemisphere
                    var sample_dir = Vector.normalize(Vector.add3(
                            Vector.scale(Math.random(), surface_n),
                            Vector.scale(Math.random()*2 - 1, surface_x),
                            Vector.scale(Math.random()*2 - 1, surface_y)));
                    var sample_factor = Vector.dot(sample_dir, surface_n);
                    if (sample_factor < 0 || sample_factor > 1) 
                        throw "Random sample direction fails.";

                    // Skip small light
                    if (sample_factor < 0.2) continue; // Skip small light 

                    // Sampled color
                    var sample_c = newTracer(tempIsect.xyz, sample_dir, inputSpheres,
                            triangles, inputLights, 0.2, true, 1, true);
                    
                    for (var i=0; i<3; ++i) {
                        tempColor[i] += sample_factor*sample_c[i];
                    }
                }
                for (var i=0; i<3; ++i) {
                    tempColor[i] /= numSample;
                    c[i] = 0.8*c[i] + 0.7*tempColor[i];
                }
            }
            return c;
        } // End intersect with sphere
        else {
            c = shadeIsectTriangle(tempIsect, tempTriangle, inputLights, inputSpheres, addAmbient);
            // Russian Roulette
            if (Math.random() < russian) {
                var tempColor = new Color(0,0,0,0);
                for (var s=0; s<numSample; ++s) {
                    var DirNormDot = Vector.dot(tempTriangle.n, Dir);
                    // Reverse the norm if point to the same direction as Dir
                    if (DirNormDot > 0) 
                        tempTriangle.n = Vector.scale(-1.0, tempTriangle.n);

                    // Finding orthogonal vector on surface
                    var surface_x = Vector.cross(tempTriangle.n, Vector.normalize(Dir));
                    var surface_y = Vector.cross(tempTriangle.n, surface_x);

                    // Random sampling over hemisphere
                    var sample_dir = Vector.normalize(Vector.add3(
                                Vector.scale(Math.random(), tempTriangle.n),
                                Vector.scale(Math.random()*2 - 1, surface_x),
                                Vector.scale(Math.random()*2 - 1, surface_y)));

                    var sample_factor = Vector.dot(sample_dir, tempTriangle.n);
                    if (sample_factor < 0 || sample_factor > 1) 
                        throw "Random sampling direction fail triangle";

                    var sample_c = newTracer(tempIsect.xyz, sample_dir, inputSpheres,
                            triangles, inputLights, 0.2, true, 1, true);

                    for (var i=0; i<3; ++i) {
                        tempColor[i] += sample_c[i];
                    }
                }
                for (var i=0; i<3; ++i) {
                    tempColor[i] /= numSample;
                    c[i] = 0.8*c[i] + 0.7*tempColor[i];
                }
            }
        } // End intersect with triangle
        return c;
    }
    catch(e) {
        console.log(e);
        return new Color(25,25,25,255);
    }
}


function newTracer(viewpoint, Dir, inputSpheres, triangles, inputLights, russian, frontWall, numSample, addAmbient) {
    try {
        var hit = hitCheck(viewpoint, Dir, inputSpheres, triangles, inputLights, frontWall);
        //return ({"isect": tempIsect, "hitSphere": hitSphere, "sphereIdx":sphereIdx, "triangle":tempTriangle});
        if (!(hit.isect.exists)) // Not hit anything
            throw "new tracer not hit.";
        // Hit something
        var color = radiance(Dir, hit.isect, hit.hitSphere, hit.sphereIdx, 
                hit.triangle, triangles, inputLights, inputSpheres, russian, 
                numSample, addAmbient);
        return color;
    }
    catch(e) {
        //console.log(e);
        //return (new Color(255,255,255,255));
        return (new Color(90,90,90,255));
    }
}



function tracer(viewpoint, Dir, inputSpheres, triangles, inputLights, recursive, frontWall, numSample) {
    try {
        var closestT = Number.MAX_VALUE;
        var c = new Color(0,0,0,255);
        var numSpheres = inputSpheres.length;
        var numTriangles = triangles.length;
        // Ignore front wall
        if (!(frontWall)) 
            numTriangles -= 2;
        //Dir.copy(Vector.subtract(new Vector(wx, wy, WIN_Z), Eye));
        // Keep the intersection information, compute color afterward
        var hitSphere = false;
        var tempIsect = {"exists":false};
        var sphereIdx = -1;
        var tempTriangle;
        // Check spheres
        for (var s=0; s<numSpheres; s++) {
            isect = raySphereIntersect([viewpoint,Dir],inputSpheres[s],1); 
            if (isect.exists) {// there is an intersect
                if (isect.t < closestT) { // it is the closest yet
                    hitSphere = true;
                    closestT = isect.t; // record closest t yet
                    tempIsect = isect;
                    sphereIdx = s;
                } // end if closest yet
            }
        } // done all spheres
        // Check trangles
        for (var t=0; t<numTriangles; t++) { // loop triangles
            isect = rayTriangleIntersect([viewpoint, Dir], triangles[t], 1);
            if (isect.exists) {
                if (isect.t < closestT) {
                    hitSphere = false;
                    closestT = isect.t;
                    tempIsect = isect;
                    tempTriangle = triangles[t];
                }
                // Triangle can break since cannot hit two triangles in this case.
                break;
            }
        } // done all triangles

        var black = new Color(0, 0, 0, 255);

        // Compute color 
        if (!(tempIsect.exists)) {
            throw "tracer not hit anything";
        }
        else if (hitSphere) {
            // Find the color on the intersection
            c = shadeIsect(tempIsect, sphereIdx, inputLights, inputSpheres, true);

            if (recursive) {
                var tempColor = new Color(0,0,0,0);
                for (var s=0; s<numSample; ++s) {
                    // Find the intersection surface
                    var sphereCenter = new Vector(
                            inputSpheres[sphereIdx].x,
                            inputSpheres[sphereIdx].y,
                            inputSpheres[sphereIdx].z);
                    var surface_n = Vector.normalize(
                            Vector.subtract(tempIsect.xyz, sphereCenter)); 
                    var surface_x = Vector.cross(Vector.normalize(Dir), surface_n);
                    var surface_y = Vector.cross(surface_x, surface_n);

                    // Random sampling
                    var sample_dir = Vector.normalize(Vector.add3(
                            Vector.scale(Math.random(), surface_n),
                            Vector.scale(Math.random()*2 - 1, surface_x),
                            Vector.scale(Math.random()*2 - 1, surface_y)));
                    var assertDot = Vector.dot(sample_dir, surface_n);
                    var test_norm = Vector.cross(surface_x, surface_y);
                    if (assertDot < 0 || assertDot > 1) 
                        throw "assertDot in sphere fails";
                    //sphereCenter.toConsole();
                    //surface_n.toConsole();
                    //sample_dir.toConsole();

                    // Sampled color
                    var sample_c = tracer(tempIsect.xyz, sample_dir, inputSpheres, 
                            triangles, inputLights, false, true, 1);
                    if (sample_c[1] == 255) 
                        // It shouldn't happen.
                        sample_c = black;
                    
                    for (var i=0; i<3; ++i) {
                        tempColor[i] += sample_c[i];
                    }
                }
                for (var i=0; i<3; ++i) {
                    tempColor[i] /= numSample;
                    c[i] += 0.7*tempColor[i];
                }
            }
        } // End intersect with sphere
        else {
            c = shadeIsectTriangle(tempIsect, tempTriangle, inputLights, inputSpheres, addAmbient);
            if (recursive) {
                var tempColor = new Color(0,0,0,0);
                for (var s=0; s<numSample; ++s) {
                    var DirNormDot = Vector.dot(tempTriangle.n, Dir);
                    // Reverse the norm if point to the same direction as Dir
                    if (DirNormDot > 0) 
                        tempTriangle.n = Vector.scale(-1.0, tempTriangle.n);

                    // Random sampling
                    var surface_x = Vector.cross(tempTriangle.n, Vector.normalize(Dir));
                    var surface_y = Vector.cross(tempTriangle.n, surface_x);

                    // Random sampling
                    var sample_dir = Vector.normalize(Vector.add3(
                                Vector.scale(Math.random(), tempTriangle.n),
                                Vector.scale(Math.random()*2 - 1, surface_x),
                                Vector.scale(Math.random()*2 - 1, surface_y)));
                    var sample_c = tracer(tempIsect.xyz, sample_dir, inputSpheres, 
                            triangles, inputLights, false, true, 1);
                    if (sample_c[1] == 1) 
                        sample_c = black;

                    for (var i=0; i<3; ++i) {
                        tempColor[i] += sample_c[i];
                    }
                }
                for (var i=0; i<3; ++i) {
                    tempColor[i] /= numSample;
                    c[i] += 0.7*tempColor[i];
                }
            }
        } // End intersect with triangle
        return c;
    }
    catch(e) {
        console.log(e);
        return new Color(25,25,25,255);
    }
}

// Cast everything
function rayCastAll(context, numberSample, russian) {
    var inputTriangles = getJSONFile(INPUT_TRIANGLES_URL, "triangles");
    var inputLights = getJSONFile(INPUT_LIGHTS_URL, "lights");
    var inputSpheres = getJSONFile(INPUT_SPHERES_URL, "spheres");

    var w = context.canvas.width;
    var h = context.canvas.height;
    var imagedata = context.createImageData(w, h);

    var numTriangles = inputTriangles.length;
    var numSpheres   = inputSpheres.length;

    var triangles = new Array(numTriangles);
    for (var i=0; i<numTriangles; i++) {
        triangles[i] = new Triangle(inputTriangles[i]);
    }

    // Loop over pixels 
    var x = 0; var y = 0;
    var Dir = new Vector(0,0,0);
    var c = new Color(0,0,0,0);
    // Actual position in image
    var wx = WIN_LEFT;
    var wxd = (WIN_RIGHT - WIN_LEFT) * 1/(w-1);
    var wy = WIN_TOP;
    var wyd = (WIN_BOTTOM - WIN_TOP) * 1/(h-1);
    for (y=0; y<h; y++) { // Loop over pixel
        wx = WIN_LEFT;
        for (x=0; x<h; x++) {
            Dir.copy(Vector.subtract(new Vector(wx, wy, WIN_Z), Eye));
            // Trace the light
            //c = tracer(Eye, Dir, inputSpheres, triangles, inputLights, true, false, numberSample);
            c = newTracer(Eye, Dir, inputSpheres, triangles, inputLights, russian, false, numberSample);
            //c.toConsole("returned color: ");
            // Draw color in that pixel
            drawPixel(imagedata,x,y,c);
            wx += wxd;
        }
        wy += wyd;
    }
    context.putImageData(imagedata, 0, 0);
} // End of rayCastAll

// use ray casting with triangles to get pixel colors
function rayCastTriangles(context) {
    var inputTriangles = getJSONFile(INPUT_TRIANGLES_URL, "triangles");
    var inputLights = getJSONFile(INPUT_LIGHTS_URL, "lights");
    var w = context.canvas.width;
    var h = context.canvas.height;
    var imagedata = context.createImageData(w, h);
    if (inputTriangles != String.null) {
        var x = 0; var y = 0;
        var n = inputTriangles.length;
        var Dir = new Vector(0, 0, 0);
        var closestT = Number.MAX_VALUE;
        var c = new Color(0,0,0,0);
        var isect = {};
        //console.log("Number of triangles: "+n);

        // Build Triangle array
        var triangles = new Array(n);
        for (var i=0; i<n; ++i) {
            triangles[i] = new Triangle(inputTriangles[i]);
        }
        triangles[0].n.toConsole();

        // Loop over the pixels and triangles
        var wx = WIN_LEFT;
        var wxd = (WIN_RIGHT - WIN_LEFT) * 1/(w-1);
        var wy = WIN_TOP;
        var wyd = (WIN_BOTTOM - WIN_TOP) * 1/(h-1);
        for (y=0; y<h; y++) {
            wx = WIN_LEFT;
            for (x=0; x<h; x++) {
                closestT = Number.MAX_VALUE;
                c.change(0,0,0,255);
                Dir.copy(Vector.subtract(new Vector(wx, wy, WIN_Z), Eye));
                //Dir.toConsole("Dir: ");
                for (var t=0; t<n; t++) { // loop triangles
                    isect = rayTriangleIntersect([Eye, Dir], triangles[t], 1);
                    //console.log(isect.exists);
                    if (isect.exists) {
                        //console.log("intersect");
                        if (isect.t < closestT) {
                            closestT = isect.t;
                            c = shadeIsectTriangle(isect, triangles[t], inputLights, true);
                        }
                        break;
                    }
                } // done all triangles
                drawPixel(imagedata,x,y,c);
                wx += wxd;
            } // x increase
            wy += wyd;
        } // y increase
        context.putImageData(imagedata, 0, 0);
    }
} // end ray cast triangles

// use ray casting with spheres to get pixel colors
function rayCastSpheres(context) {
    var inputSpheres = getJSONFile(INPUT_SPHERES_URL,"spheres");
    var inputLights = getJSONFile(INPUT_LIGHTS_URL,"lights");
    var w = context.canvas.width;
    var h = context.canvas.height;
    var imagedata = context.createImageData(w,h);
    // console.log("casting rays");

    if (inputSpheres != String.null) { 
        var x = 0; var y = 0; // pixel coord init
        var n = inputSpheres.length; // the number of spheres
        var Dir = new Vector(0,0,0); // init the ray direction
        var closestT = Number.MAX_VALUE; // init the closest t value
        var c = new Color(0,0,0,0); // init the pixel color
        var isect = {}; // init the intersection
        //console.log("number of spheres: " + n);

        // Loop over the pixels and spheres, intersecting them
        var wx = WIN_LEFT; // init world pixel xcoord
        var wxd = (WIN_RIGHT-WIN_LEFT) * 1/(w-1); // world pixel x differential
        var wy = WIN_TOP; // init world pixel ycoord
        var wyd = (WIN_BOTTOM-WIN_TOP) * 1/(h-1); // world pixel y differential
        for (y=0; y<h; y++) {
            wx = WIN_LEFT; // init w
            for (x=0; x<h; x++) {
                closestT = Number.MAX_VALUE; // no closest t for this pixel
                c.change(0,0,0,255); // set pixel to background color
                Dir.copy(Vector.subtract(new Vector(wx,wy,WIN_Z),Eye)); // set ray direction
                //Dir.toConsole("Dir: ");
                for (var s=0; s<n; s++) {
                // for (var s=0; s<1; s++) {
                    isect = raySphereIntersect([Eye,Dir],inputSpheres[s],1); 
                    if (isect.exists) // there is an intersect
                        if (isect.t < closestT) { // it is the closest yet
                            closestT = isect.t; // record closest t yet
                            c = shadeIsect(isect,s,inputLights,inputSpheres, true); 
                        } // end if closest yet
                } // end for spheres
                drawPixel(imagedata,x,y,c);
                wx += wxd; 
                //console.log(""); // blank per pixel
            } // end for x
            wy += wyd; 
        } // end for y
        context.putImageData(imagedata, 0, 0);
    } // end if spheres found
} // end ray cast spheres

// given a pixel position, calculate x and y pixel and world coords
function getPixelLocat(pixelNum, w, h) {
    var y = Math.floor(pixelNum/w);
    var x = pixelNum - y*w; 
    var wx = WIN_LEFT + x/w * (WIN_RIGHT - WIN_LEFT);
    var wy = WIN_TOP + y/h * (WIN_BOTTOM - WIN_TOP);
    
    //console.log("pixelNum: "+pixelNum+", x:"+x+", y:"+y+", wx:"+wx+", wy:"+wy);
    
    return ({"x": x, "y": y, "wx": wx, "wy": wy});
}

// use frameless ray casting with spheres to get pixel colors
function framelessRayCastSpheres(context) {
    var inputSpheres = getJSONFile(INPUT_SPHERES_URL,"spheres");
    var inputLights = getJSONFile(INPUT_LIGHTS_URL,"lights");

    if ((inputSpheres != String.null) && (inputLights != String.null)) { 
        var n = inputSpheres.length; // the number of spheres
        var w = context.canvas.width;
        var h = context.canvas.height;
        var numPixels = w*h; 
        var pixelOrder = randPermutation(numPixels); // rand order for visiting pixels
        var imagedata = context.createImageData(1,1); //  just one pixel at a time
        imagedata.data[3] = 255; // pixels are always opaque
        var intervalID = 0; // the ID returned by the last setInterval call
        var p = 0; // start at first rand pixel 
        //console.log("num pixels: "+numPixels);
        //console.log("number of spheres: " + n);
        
        // update a frame with the next set of random rays
        function framelessUpdate() { 
            var endTime = performance.now() + 0.9;
            var pixelLocat; // where the pixel is located on the screen
            var Dir = new Vector(0,0,0); // init the ray direction
            var closestT = Number.MAX_VALUE; // init the closest t value
            var isect = {}; // init the intersection
            var c = new Color(0,0,0,255); // declare the pixel color

            // Loop over the pixels and spheres, intersecting them
            while (performance.now() < endTime) {
                closestT = Number.MAX_VALUE; // no closest t for this pixel
                c.change(0,0,0,255); // set pixel to background color
                pixelLocat = getPixelLocat(pixelOrder[p],w,h); // get pixel location
                Dir.copy(Vector.subtract(new Vector(pixelLocat.wx,pixelLocat.wy,WIN_Z),Eye)); // set ray direction
                //Dir.toConsole("Dir: ");
                for (var s=0; s<n; s++) { // for each sphere
                    // for (var s=0; s<1; s++) {
                    isect = raySphereIntersect([Eye,Dir],inputSpheres[s],1); 
                    if (isect.exists) // there is an intersect
                        if (isect.t < closestT) { // it is the closest yet
                            closestT = isect.t; // record closest t yet
                            c = shadeIsect(isect,s,inputLights,inputSpheres, true); 
                        } // end if closest yet
                } // end for spheres
                imagedata.data[0] = c[0];
                imagedata.data[1] = c[1];
                imagedata.data[2] = c[2];
                context.putImageData(imagedata,pixelLocat.x,pixelLocat.y);
                p++; // next pixel
                if (p >= numPixels) { // back to first pixel if finished
                    p = 0;
                    //console.log("restart rand pixel order: p=0");
                } // end if reached max pixels
            } // end while in frame
        } // end frameless update
        
        // update the current frame using frameless updates
        function frameUpdate() {
            
                // if frameless update is in progress, interrupt it.
            if (intervalID != 0) // an update is in progress 
                window.clearInterval(intervalID); 
            
                // now the end of frame is over, do 
            window.setInterval(framelessUpdate,1);
            window.requestAnimationFrame(frameUpdate);
        } // end frameUpdate
    
        window.requestAnimationFrame(frameUpdate);
    } // end if spheres found 
} // end frameless ray cast spheres


/* constants and globals */

const WIN_Z = 0; const WIN_Z_BACK = -1;
const WIN_LEFT = 0; const WIN_RIGHT = 1;
const WIN_BOTTOM = 0; const WIN_TOP = 1; 
const INPUT_SPHERES_URL = 
    "https://junioryi.github.io/csc562-programIO/spheres.json";
    //"https://junioryi.github.io/csc562-programIO/test_spheres.json";
    //"https://ncsucgclass.github.io/prog1/spheres.json";
    //"https://pages.github.ncsu.edu/bwatson/introcg-prog1/spheres.json";
const INPUT_LIGHTS_URL = 
    "https://junioryi.github.io/csc562-programIO/lights.json";
    //"https://ncsucgclass.github.io/prog1/lights.json";
    //"https://pages.github.ncsu.edu/bwatson/introcg-prog1/lights.json";
const INPUT_TRIANGLES_URL =
    "https://junioryi.github.io/csc562-programIO/triangles.json";
        
var Eye = new Vector(0.5,0.5,-0.5); // set the eye position


/* main -- here is where execution begins after window load */

function main() {

    // Get the canvas and context
    var canvas = document.getElementById("viewport"); 
    var context = canvas.getContext("2d");
 
    // Create the image
    //drawRandPixels(context);
      // shows how to draw pixels
    
    //drawRandPixelsInInputSpheres(context);
      // shows how to draw pixels and read input file
      
    //drawInputSpheresUsingArcs(context);
      // shows how to read input file, but not how to draw pixels
      
    //rayCastSpheres(context); 
    //rayCastTriangles(context);
    var numberSample = 20;
    rayCastAll(context, numberSample, 1.0);
    
    //framelessRayCastSpheres(context);
}
