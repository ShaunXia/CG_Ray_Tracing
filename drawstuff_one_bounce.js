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




    //A = [a1, a2, a3] and B = [b1, b2, b3]
    //cross(A, B) = [ a2 * b3 - a3 * b2, a3 * b1 - a1 * b3, a1 * b2 - a2 * b1 ]
    static cross(v1,v2) {
    return (new Vector(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x));
    }

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
    
} // end Vector class


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
function raySphereIntersect(ray,sphere,clipVal) {
    try {
        if (!(ray instanceof Array) || !(sphere instanceof Object))
            throw "RaySphereIntersect: ray or sphere are not formatted well";
        else if (ray.length != 2)
            throw "RaySphereIntersect: badly formatted ray";
        else { // valid params
            var a = Vector.dot(ray[1],ray[1]); // dot(D,D)

            var origMctr = Vector.subtract(ray[0],new Vector(sphere.x,sphere.y,sphere.z)); // E-C

            var b = 2 * Vector.dot(ray[1],origMctr); // 2 * dot(D,E-C)
            var c = Vector.dot(origMctr,origMctr) - sphere.r*sphere.r; // dot(E-C,E-C) - r^2

            // Ray-Sphere intersection equation and solution.
            // t = -(d*v)+-sqrt((d*v)^2-v^2+r^2)

            // if (clipVal == 0) {
            //     ray[0].toConsole("ray.orig: ");
            //     ray[1].toConsole("ray.dir: ");
            //     console.log("a:"+a+" b:"+b+" c:"+c);
            // } // end debug case
        
            var qsolve = solveQuad(a,b,c);

            var sphereCenter = new Vector(sphere.x,sphere.y,sphere.z);

            if (qsolve.length == 0) 
                throw "no intersection";
            else if (qsolve.length == 1) { 
                if (qsolve[0] < clipVal)
                    throw "intersection too close";
                else {
                    var isect = Vector.add(ray[0],Vector.scale(qsolve[0],ray[1]));
                    var N = Vector.normalize(Vector.subtract(isect,sphereCenter)); // surface normal
                    //console.log("t: "+qsolve[0]);
                    //isect.toConsole("intersection: ");
                    return({"exists": true, "tri_sphere":"sphere","xyz": isect,"t": qsolve[0],"normal":N});  
                } // one unclipped intersection
            } else if (qsolve[0] < clipVal) {
                if (qsolve[1] < clipVal)
                    throw "intersections too close";
                else { 
                    var isect = Vector.add(ray[0],Vector.scale(qsolve[1],ray[1]));
                    var N = Vector.normalize(Vector.subtract(isect,sphereCenter)); 
                    //console.log("t2: "+qsolve[1]);
                    //isect.toConsole("intersection: ");
                    return({"exists": true, "tri_sphere":"sphere","xyz": isect,"t": qsolve[1],"normal" : N});  
                } // one intersect too close, one okay
            } else {
                var isect = Vector.add(ray[0],Vector.scale(qsolve[0],ray[1]));
                var N = Vector.normalize(Vector.subtract(isect,sphereCenter)); 
                //console.log("t1: "+qsolve[0]);
                //isect.toConsole("intersection: ");
                return({"exists": true,"tri_sphere":"sphere", "xyz": isect,"t": qsolve[0],"normal" : N});  
            } // both not too close
        } // end if valid params
    } // end try
                    

    catch(e) {
        //console.log(e);
        return({"exists": false,"tri_sphere":NaN, "xyz": NaN, "t": NaN,"normal" : NaN});
    }
} // end raySphereIntersect

function rayTriangleIntersect(ray,triangle,clipVal)
{
	a = new Vector(triangle.vertices[0][0],triangle.vertices[0][1],triangle.vertices[0][2]);
	b= new  Vector(triangle.vertices[1][0],triangle.vertices[1][1],triangle.vertices[1][2]);
	c = new  Vector(triangle.vertices[2][0],triangle.vertices[2][1],triangle.vertices[2][2]);

	var ab = Vector.subtract(b,a);
	var ac = Vector.subtract(c,a);
	//var normal = Vector.normalize(Vector.cross(ab,ac)); // get the normal of triangle plane
    var normal = new Vector(triangle.normals[0],triangle.normals[1],triangle.normals[2]); // surface normal

    //origin = ray[0]
	var t = Vector.dot(normal,Vector.subtract(a,ray[0]))/Vector.dot(normal,ray[1]);	//normal.dot(a.subtract(origin)) / normal.dot(ray);

	if (t > 0)
	{
		var hit = Vector.add(ray[0],Vector.scale(t,ray[1]));//origin.add(ray.multiply(t));
		var toHit = Vector.subtract(hit,a);//hit.subtract(a);
		var dot00 = Vector.dot(ac,ac);//ac.dot(ac);
		var dot01 = Vector.dot(ac,ab);//ac.dot(ab);
		var dot02 = Vector.dot(ac,toHit);//ac.dot(toHit);
		var dot11 = Vector.dot(ab,ab);//ab.dot(ab);
		var dot12 = Vector.dot(ab,toHit);//ab.dot(toHit);
		var divide = dot00 * dot11 - dot01 * dot01;
		var u = (dot11 * dot02 - dot01 * dot12) / divide;
		var v = (dot00 * dot12 - dot01 * dot02) / divide;
		if (u >= 0 && v >= 0 && u + v <= 1) return ({"exists": true,"tri_sphere":"triangle", "xyz": hit,"t": t,"normal" : normal});
  	}

  return({"exists": false,"tri_sphere":NaN, "xyz": NaN, "t": NaN, "normal":NaN});

}

function getIntersectInScene(ray,triangles,spheres,clipVal)
{
            //Dir.toConsole("Dir: ");
        var isect;
        var closestT = Number.MAX_VALUE; 
        var hitID=0;
        var finalIsect={};
        for (var s=0; s<spheres.length; s++) {
        // for (var s=0; s<1; s++) {
            isect = raySphereIntersect(ray,inputSpheres[s],clipVal); 
            if (isect.exists) // there is an intersect
                if (isect.t < closestT) { // it is the closest yet
                    closestT = isect.t; // record closest t yet
                    finalIsect= isect;
                    hitID=s;
                } // end if closest yet
        } // end for spheres



        for (var s=0; s<triangles.length; s++) 
        {
            	//console.log("Triangle:"+s);
                isect_tri = rayTriangleIntersect(ray, inputTriangles[s],clipVal); 
                if (isect_tri.exists) // there is an intersect
                    if (isect_tri.t < closestT) 
                    { // it is the closest yet
                        closestT = isect_tri.t; // record closest t yet
                        isect=isect_tri;
                        hitID=s;
                        finalIsect=isect_tri; 

                    } // end if closest yet
        } // end for spheres

        // if(isect.exists&&isect["tri_sphere"]=="sphere")
        // {
        // 	console.log(hitID);
        // }

        finalIsect["hitID"]=(hitID);


        return finalIsect;

}

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
    var w = context.canvas.width;
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
                drawPixel(imagedata,
                    cx+Math.round(x*sphereRadius),
                    cy+Math.round(y*sphereRadius),c);
                //console.log("color: ("+c.r+","+c.g+","+c.b+")");
                //console.log("x: "+Math.round(w*inputSpheres[s].x));
                //console.log("y: "+Math.round(h*inputSpheres[s].y));
            } // end for pixels in sphere
        } // end for spheres
        context.putImageData(imagedata, 0, 0);
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
    // console.log("testing for occlusions");
    

    while ((!lightOccluded) && (s<spheres.length)) {
    	if(s==isectSphere) // do not check current sphere itself.
    	{
    		s++;
    		continue;
    	}
        occluderIsect = raySphereIntersect([isectPos,L],spheres[s],0);
        // console.log("oisect: "+occluderIsect);
        if (!occluderIsect.exists) { // no intersection
            s++; // on to next sphere
        } else if (occluderIsect.t > 1) { // light in front of intersection
            s++; // on to next sphere
        } else {
            lightOccluded = true;
            // console.log("occlusion found from sphere "+isectSphere+" to "+s);
        } // end if occlusion found
    } // while all lights after one intersected by eye
    
    return(lightOccluded);
} // end is light occluded

function reflectedRay(isect,hitID,lights,spheres,lev,ray)
{
	try {
        if (   !(isect instanceof Object) || !(typeof(hitID) === "number") 
            || !(lights instanceof Array) || !(spheres instanceof Array))
            throw "shadeIsect: bad parameter passed";
        else {
            var c = new Color(0,0,0,255); // init the sphere color to black
            var sphere = spheres[hitID]; // sphere intersected by eye
            // console.log("shading pixel");

            // add light for each source
            var lightOccluded = false; // if an occluder is found
            var Lloc = new Vector(0,0,0);
            for (var l=0; l<lights.length; l++) {

                // add in the ambient light
                c[0] += lights[l].ambient[0] * sphere.ambient[0]; // ambient term r
                c[1] += lights[l].ambient[1] * sphere.ambient[1]; // ambient term g
                c[2] += lights[l].ambient[2] * sphere.ambient[2]; // ambient term b

                // var L = Vector.subtract(ray,isect.xyz); // light vector unnorm'd
                var L = ray;
                // L.toConsole("L: ");
                // console.log("isect: "+isect.xyz.x+", "+isect.xyz.y+", "+isect.xyz.z);
                
                // if light isn't occluded
                if (!isLightOccluded(L,isect.xyz,hitID,spheres)) {
                    // console.log("no occlusion found");
                    
                    // add in the diffuse light
                    var sphereCenter = new Vector(sphere.x,sphere.y,sphere.z);
                    var N = Vector.normalize(Vector.subtract(isect.xyz,sphereCenter)); // surface normal

                    //var R = Vector.add(Vector.scale(-2*Vector.dot(N,L),N),L);
                    //var r_ray = Vector.subtract(R,isect.xyz);
                    // console.log("ori:"+L.x);
                    // console.log("r_ray:"+R.x);
                    // c=Vector.add(c,shadeIsect(isect,hitID,R,inputSpheres,2)); 

                    var diffFactor = Math.max(0,Vector.dot(N,Vector.normalize(L)));

                    if (diffFactor > 0) {
                        c[0] += lights[l].diffuse[0] * sphere.diffuse[0] * diffFactor;
                        c[1] += lights[l].diffuse[1] * sphere.diffuse[1] * diffFactor;
                        c[2] += lights[l].diffuse[2] * sphere.diffuse[2] * diffFactor;
                    } // end nonzero diffuse factor



                    // add in the specular light
                    // var V = Vector.normalize(Vector.subtract(Eye,isect.xyz)); // view vector
                    // var H = Vector.normalize(Vector.add(L,V)); // half vector
                    // var specFactor = Math.max(0,Vector.dot(N,H)); 
                    // if (specFactor > 0) {
                    //     var newSpecFactor = specFactor;
                    //     for (var s=1; s<spheres[hitID].n; s++) // mult by itself if needed
                    //         newSpecFactor *= specFactor;
                    //     c[0] += lights[l].specular[0] * sphere.specular[0] * newSpecFactor; // specular term
                    //     c[1] += lights[l].specular[1] * sphere.specular[1] * newSpecFactor; // specular term
                    //     c[2] += lights[l].specular[2] * sphere.specular[2] * newSpecFactor; // specular term
                    // } // end nonzero specular factor
                    


                } // end if light not occluded
            } // end for lights
            
          

            return(c);
        } // if have good params
    } // end throw
    
    catch(e) {
        console.log(e);
        return(Object.null);
    }
}

function radiance(ray,lights,spheres,triangles,lev,brdf,russian_flag=0)
{
	var isect={};
    isect=getIntersectInScene(ray,inputTriangles,inputSpheres,1);
	var c = new Color(0,0,0,255);
    if (isect.exists) {
		var dir_color=dirIllum(isect,lights,brdf);
    	var indir_color = indirIllum(isect,lev);
		c[0]=dir_color[0]+indir_color[0];
		c[1]=dir_color[1]+indir_color[1];
		c[2]=dir_color[2]+indir_color[2];
    	return (c);
    	
    }

    // no intersection
    return(c);
	
	
}

function dirIllum(isect,lights,brdf)
{
	var c = new Color(0,0,0,255); 
	var hitID=isect["hitID"];
	if (isect["tri_sphere"]=="sphere") {
		c=shadeIsect(isect,hitID,lights,inputSpheres,1,brdf); 
	}
	else
		c=shadeTriIsect(isect,hitID,lights,inputTriangles[hitID],inputSpheres,1,brdf);

	return(c);
}

function indirIllum(isect,lev)
{
	

	var sampleNum=100;
	var indirC = new Color(0,0,0,255); 
	var estRadiance=new Color(0,0,0,255); 
	// var p_bounce=0.5
	// if (Math.random()>p_bounce) {return(estRadiance);}

	if (lev>1) {return(estRadiance);}

    var smaplePoint=[];
    var sampleTheta=[];
    var scount=0;
    var totalcount=0;

    while(1)
    {
        if (scount==sampleNum) {
            break;
        }
        var randomX = Math.random();
        var randomY = Math.random();
        var randomZ = Math.random();
        var randomLoc = new Vector(randomX, randomY, randomZ);
        var Dir = Vector.subtract(randomLoc,isect.xyz);

        var R = Dir.x*Dir.x+Dir.y*Dir.y+Dir.z*Dir.z;
        // console.log(isect.tri_sphere+"--isect:"+isect.xyz.x+" "+isect.xyz.y+" "+isect.xyz.z);
        var cosTheta = Vector.dot(Vector.normalize(isect.normal),Vector.normalize(Dir));
        var acos = Math.acos(cosTheta);

        if(R<0.04&&acos<(Math.PI/2))//
            {
                ++scount;
                //console.log(scount);
                smaplePoint.push(Dir);
                sampleTheta.push(cosTheta);
            }
            // ++totalcount;

            // console.log("acos:"+acos+" R:"+R);
    }
    var brdf=1;

	for (var i = 0; i < smaplePoint.length; i++)
	{

    indirC = radiance([isect.xyz,smaplePoint[i]],inputLights,inputSpheres,inputTriangles,lev+1,brdf,0);

    //var diffFactor = Math.max(0,Vector.dot(N,Vector.normalize(L)));

    estRadiance[0] +=  Math.min(1,indirC[0]*sampleTheta[i]); // clamp max value to 1
    estRadiance[1] +=  Math.min(1,indirC[1]*sampleTheta[i]); // clamp max value to 1
    estRadiance[2] +=  Math.min(1,indirC[2]*sampleTheta[i]); // clamp max value to 1

	}
	//sampleNum*=p_bounce;
	estRadiance[0]/=sampleNum;
	estRadiance[1]/=sampleNum;
	estRadiance[2]/=sampleNum;
    //console.log(estRadiance);
	// console.log(estRadiance);
	return(estRadiance);
}


function indirIllumWithRussian(isect,lev)
{
	

	var sampleNum=100;
	var indirC = new Color(0,0,0,255); 
	var estRadiance=new Color(0,0,0,255); 
	var p_bounce=0.5
	if (Math.random()<p_bounce) {return(estRadiance);}

	//if (lev>1) {return(estRadiance);}

    var smaplePoint=[];
    var sampleTheta=[];
    var scount=0;
    var totalcount=0;

    while(1)
    {
        if (scount==sampleNum) {
            break;
        }
        var randomX = Math.random();
        var randomY = Math.random();
        var randomZ = Math.random();
        var randomLoc = new Vector(randomX, randomY, randomZ);
        var Dir = Vector.subtract(randomLoc,isect.xyz);

        var R = Dir.x*Dir.x+Dir.y*Dir.y+Dir.z*Dir.z;
        // console.log(isect.tri_sphere+"--isect:"+isect.xyz.x+" "+isect.xyz.y+" "+isect.xyz.z);
        var cosTheta = Vector.dot(Vector.normalize(isect.normal),Vector.normalize(Dir));
        var acos = Math.acos(cosTheta);

        if(R<0.04&&acos<(Math.PI/2))//
            {
                ++scount;
                //console.log(scount);
                smaplePoint.push(Dir);
                sampleTheta.push(cosTheta);
            }
            // ++totalcount;

            // console.log("acos:"+acos+" R:"+R);
    }
    var brdf=1;

	for (var i = 0; i < smaplePoint.length; i++)
	{

		// var r1=Math.random();
		// var r2=Math.random();
		// var sinTheta = Math.sqrt(1 - r1*r1);
		// var phi = 2 * 3.14 * r2;
		// var x = sinTheta * Math.cos(phi);
		// var z = sinTheta * Math.sin(phi);
		// var sample = new Vector(x,r1,z);

    indirC = radiance([isect.xyz,smaplePoint[i]],inputLights,inputSpheres,inputTriangles,lev+1,brdf,1);

    //var diffFactor = Math.max(0,Vector.dot(N,Vector.normalize(L)));

    estRadiance[0] +=  Math.min(1,indirC[0]*0.3); // clamp max value to 1
    estRadiance[1] +=  Math.min(1,indirC[1]*0.3); // clamp max value to 1
    estRadiance[2] +=  Math.min(1,indirC[2]*0.3); // clamp max value to 1

	}
	sampleNum*=p_bounce;
	estRadiance[0]/=sampleNum;
	estRadiance[1]/=sampleNum;
	estRadiance[2]/=sampleNum;
    //console.log(estRadiance);
	// console.log(estRadiance);
	return(estRadiance);
}

// color the passed intersection and sphere
function shadeIsect(isect,hitID,lights,spheres,lev,brdf) {
    try {
        if (   !(isect instanceof Object) || !(typeof(hitID) === "number") 
            || !(lights instanceof Array) || !(spheres instanceof Array))
            throw "shadeIsect: bad parameter passed";
        else {
            var c = new Color(0,0,0,255); // init the sphere color to black
            var sphere = spheres[hitID]; // sphere intersected by eye
            // console.log("shading pixel");

            // add light for each source
            var lightOccluded = false; // if an occluder is found
            var Lloc = new Vector(0,0,0);
            for (var l=0; l<lights.length; l++) {

                if(brdf==0)
                {
                // add in the ambient light
                c[0] += lights[l].ambient[0] * sphere.ambient[0]; // ambient term r
                c[1] += lights[l].ambient[1] * sphere.ambient[1]; // ambient term g
                c[2] += lights[l].ambient[2] * sphere.ambient[2]; // ambient term b
                }
                // check each other sphere to see if it occludes light
                Lloc.set(lights[l].x,lights[l].y,lights[l].z);
                var L = Vector.subtract(Lloc,isect.xyz); // light vector unnorm'd

                // L.toConsole("L: ");
                // console.log("isect: "+isect.xyz.x+", "+isect.xyz.y+", "+isect.xyz.z);
                
                var sphereCenter = new Vector(sphere.x,sphere.y,sphere.z);
                    var N = Vector.normalize(Vector.subtract(isect.xyz,sphereCenter)); // surface normal

                    //var R = Vector.add(Vector.scale(-2*Vector.dot(N,L),N),L);
                    //var r_ray = Vector.subtract(R,isect.xyz);
                    // console.log("ori:"+L.x);
                    // console.log("r_ray:"+R.x);

                   //var r_color=reflectedRay(isect,hitID,lights,inputSpheres,2,R); 
                   //console.log("r-color:",r_color);
                   // c[0]+=r_color[0];
                   // c[1]+=r_color[1];
                   // c[2]+=r_color[2];
                // if light isn't occluded
                if (!isLightOccluded(L,isect.xyz,hitID,spheres)) {
                    // console.log("no occlusion found");
                    
                    // add in the diffuse light
                    
                    var diffFactor = Math.max(0,Vector.dot(N,Vector.normalize(L)));

                    if (diffFactor > 0) {
                        c[0] += lights[l].diffuse[0] * sphere.diffuse[0] * diffFactor;
                        c[1] += lights[l].diffuse[1] * sphere.diffuse[1] * diffFactor;
                        c[2] += lights[l].diffuse[2] * sphere.diffuse[2] * diffFactor;
                    } // end nonzero diffuse factor

                    // add in the specular light
                    var V = Vector.normalize(Vector.subtract(Eye,isect.xyz)); // view vector
                    var H = Vector.normalize(Vector.add(L,V)); // half vector
                    var specFactor = Math.max(0,Vector.dot(N,H)); 
                    if (specFactor > 0&&brdf==0) {
                        var newSpecFactor = specFactor;
                        for (var s=1; s<spheres[hitID].n; s++) // mult by itself if needed
                            newSpecFactor *= specFactor;
                        c[0] += lights[l].specular[0] * sphere.specular[0] * newSpecFactor; // specular term
                        c[1] += lights[l].specular[1] * sphere.specular[1] * newSpecFactor; // specular term
                        c[2] += lights[l].specular[2] * sphere.specular[2] * newSpecFactor; // specular term
                    } // end nonzero specular factor
                    
                } // end if light not occluded
            } // end for lights
            
          

            return(c);
        } // if have good params
    } // end throw
    
    catch(e) {
        console.log(e);
        return(Object.null);
    }
}


// color the passed intersection and triangle
function shadeTriIsect(isect,isectSphere,lights,triangles,spheres,lev,brdf) {
    try {
        if (   !(isect instanceof Object) || !(typeof(isectSphere) === "number") 
            || !(lights instanceof Array) || !(triangles instanceof Object))
            throw "shadeIsect: bad parameter passed";
        else {
            var c = new Color(0,0,0,255); // init the sphere color to black
            // var triangles = triangles[isectSphere]; // sphere intersected by eye

            // add light for each source
            var lightOccluded = false; // if an occluder is found
            var Lloc = new Vector(0,0,0);
            for (var l=0; l<lights.length; l++) {

                if(brdf==0)
                
                {// add in the ambient light
                c[0] += lights[l].ambient[0] * triangles.material.ambient[0]; // ambient term r
                c[1] += lights[l].ambient[1] * triangles.material.ambient[1]; // ambient term g
                c[2] += lights[l].ambient[2] * triangles.material.ambient[2]; // ambient term b
                }
            	//console.log(isect.xyz);
                // check each other sphere to see if it occludes light
                //Lloc.set(0.5,1,-0.5);
                 Lloc.set(lights[l].x,lights[l].y,lights[l].z);
                var L = Vector.subtract(Lloc,isect.xyz); // light vector unnorm'd
                // L.toConsole("L: ");
                // console.log("isect: "+isect.xyz.x+", "+isect.xyz.y+", "+isect.xyz.z);
                
                // if light isn't occluded
                if (!isLightOccluded(L,isect.xyz,-1,spheres,false)) {
                    // console.log("no occlusion found");
                    
                    // add in the diffuse light

                    var N = new Vector(triangles.normals[0],triangles.normals[1],triangles.normals[2]); // surface normal

                    var diffFactor = Math.max(0,Vector.dot(N,Vector.normalize(L)));

                    if (diffFactor > 0) {
                        c[0] += lights[l].diffuse[0] * triangles.material.diffuse[0] * diffFactor;
                        c[1] += lights[l].diffuse[1] * triangles.material.diffuse[1] * diffFactor;
                        c[2] += lights[l].diffuse[2] * triangles.material.diffuse[2] * diffFactor;
                    } // end nonzero diffuse factor

                    // add in the specular light
                    var V = Vector.normalize(Vector.subtract(Eye,isect.xyz)); // view vector
                    var H = Vector.normalize(Vector.add(L,V)); // half vector
                    var specFactor = Math.max(0,Vector.dot(N,H)); 
                    if (specFactor > 0&&brdf==0) {
                        var newSpecFactor = specFactor;
                        for (var s=1; s<triangles.material.n; s++) // mult by itself if needed
                            newSpecFactor *= specFactor;
                        c[0] += lights[l].specular[0] * triangles.material.specular[0] * newSpecFactor; // specular term
                        c[1] += lights[l].specular[1] * triangles.material.specular[1] * newSpecFactor; // specular term
                        c[2] += lights[l].specular[2] * triangles.material.specular[2] * newSpecFactor; // specular term
                    } // end nonzero specular factor
                    
                } // end if light not occluded
            } // end for lights

            
            return(c);
        } // if have good params
    } // end throw
    
    catch(e) {
        console.log(e);
        return(Object.null);
    }
} // end shade triIsect

// use ray casting with spheres to get pixel colors
function rayCastSpheres(context) {

    // var inputSpheres = getJSONFile(INPUT_SPHERES_URL,"spheres");
    // var inputLights = getJSONFile(INPUT_LIGHTS_URL,"lights");
    // var inputTriangles = getJSONFile(INPUT_TRIANGLES_URL,"triangles");

    var w = context.canvas.width;
    var h = context.canvas.height;
    var imagedata = context.createImageData(w,h);
    // console.log("casting rays");

    if (inputSpheres != String.null) { 
        var x = 0; var y = 0; // pixel coord init
        var n = inputSpheres.length; // the number of spheres
        var Dir = new Vector(0,0,0); // init the ray direction
        // init the closest t value
        var c = new Color(0,0,0,0); // init the pixel color
        var isect = {}; // init the intersection
        //console.log("number of spheres: " + n);
        var closestT = Number.MAX_VALUE; 
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

                c=radiance([Eye,Dir],inputLights,inputSpheres,inputTriangles,1,0);
                
	  			c[0] = 255 * Math.min(1,c[0]); // clamp max value to 1
	            c[1] = 255 * Math.min(1,c[1]); // clamp max value to 1
	            c[2] = 255 * Math.min(1,c[2]); // clamp max value to 1
                drawPixel(imagedata,x,y,c);
                wx += wxd; 
                //console.log(""); // blank per pixel
            } // end for x
            console.log("y:"+y);
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
                            c = shadeIsect(isect,s,inputLights,inputSpheres); 
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

const WIN_Z = 0;
const WIN_LEFT = 0; const WIN_RIGHT = 1;
const WIN_BOTTOM = 0; const WIN_TOP = 1; 
const INPUT_SPHERES_URL = 
    "https://pages.github.ncsu.edu/sxia4/cg/spheres_new.json";
    //"https://pages.github.ncsu.edu/bwatson/introcg-prog1/spheres.json";
const INPUT_LIGHTS_URL = 
    "https://pages.github.ncsu.edu/sxia4/cg/lights.json";
    //"https://pages.github.ncsu.edu/bwatson/introcg-prog1/lights.json";

const INPUT_TRIANGLES_URL = 
    "https://pages.github.ncsu.edu/sxia4/cg/triangles.json";
    //"https://pages.github.ncsu.edu/bwatson/introcg-prog1/lights.json";

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
      
    rayCastSpheres(context); 
    
    //framelessRayCastSpheres(context);
}
