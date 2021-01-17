/*  Various functions

    This part gather various functions which will be usefull in the rest of the code
*/

function positivMod(n, mod) {
    // classic n%mod operator gives a negative number when n < 0. This function give a positive modulo in such cases .
    return (n%mod+mod)%mod;
}

/*  polygone classe 

    Come from polygon classes project (https://github.com/benpin31/polygon). Later, I with import it as js module, 
    But still have to be train to know how it works :o

    Module polygon contains some classes to manage polygon gemetry in code. Principal classes is 
    classes polygon which contains a method sat for (Separating Axes Theorem) which compute if two polygones are separated.

    Other classes are : 
    - point, vector, straightLine : they are requisite in polygon definitions and methods
    - square which extend polygon in the special case of square
*/

/* Various functions */

function keepNDec(x, nDec) {
    /*  "round" a number to its nDec decimal */

    return Math.round(x*10 ** nDec)/(10 ** nDec);
}


/* Polygones object */

class point {
    constructor(x,y) {
        /* Point are coded ny their two coordinates */
        this.x = x ;
        this.y = y ;
    }

    sameAbsciss(other) {
        /* return true if this and other have the same abscissa */
        return this.x === other.x ;
    }

    sameOrdinate(other) {
        /* return true if this and other have the same ordinate */
        return this.y === other.y ;
    }

    equal(other) {
        /* return trus if this and other are equal (same abscissa and odinate) */
        return this.sameAbsciss(other) && this.sameOrdinate(other) ;
    }

    distance(other) {
        /* return the distance between this and other */
        return Math.sqrt( (this.x - other.y) ** 2 + (this.y - other.y) ** 2);
    }

    addVector(v) {
        /* compute a new point which is the sum of this and a vector v. v is an elemnt of class vector */
        let res = new point(this.x + v.x, this.y + v.y) ;
        return res ;
    }

    translate(v) {
        /*  translate point this of vector v. contrary to addVector : the method change directly
            the attribute of the point and return nothing */
        this.x += v.x;
        this.y += v.y;
    }

    static segmentCenter(point1, point2) {
        /*  Compute the center of th segment defined by points point1 and point2 */
       let res = new point((point1.x+point2.x)/2, (point1.y+point2.y)/2)
       return res ;
    }
}

class vector {
    constructor(M,N) {
        /*  M and N can be : 
            - two points, in that case, vector coordinate are Difference of N and M coordinates
            - two coordinates
        */
        if (M instanceof point) {
            this.x = N.x - M.x ;
            this.y = N.y - M.y ;
        } else {
            /* vectors are coded ny their two coordinates */
            this.x = M ;
            this.y = N ;
        }
    }

    is0() {
        /* return true if the vector is null vector */
        return this.x === 0 && this.y === 0 ;
    }

    sum(other) {
        /* return the sum of this and another vector v */
        let res = new vector(this.x + other.x, this.y + other.y) ;
        return res ;
    }

    product(lambda) {
        /* return the external product of this per lambda*/
        let res = new vector(lambda * this.x, lambda * this.y) ;
        return res ;
    }

    scalarProduct(other) {
        /* compute the scalar product of this and other */
        return this.x*other.x + this.y*other.y ;
    } 

    norm() {
        /* compute the norm of this */
        return Math.sqrt(this.scalarProduct(this)) ;
    }

    orthogonalVector() {
        /*  Give the unique direct orthongonal vector of this with the same norm */
        let res = new vector(-this.y, this.x)
        return res;
    }

    othogonalProjection(v) {
        /* return the orthogonal projection of v on this */
        let res = this.product(this.scalarProduct(v) / (this.norm() ** 2));
        return res ;
    }

    polarCoordinate() {
        /*  Compute the polar coordintes of this */
        return [this.norm(), Math.atan(this.y/this.x)] ;
    }
}

class straightLine {
    constructor(point1, element) {
        /*  Straight lines can be definied two ways : 
            - by two points : in that case, element must be a point
            - by a point and a director vector, in that case, element must be a vector.
            In all the case, a straightLine has only two moints attributes. If defined by a vector, the constructor
            compute the other point*/

        if (element instanceof point) {
            if (!point1.equal(element)) {
                this.point1 = point1 ;
                this.point2 = element ;
            } else {
                throw "A straight line can't be definied by two identical points" ; 
            }
        } else {
            if (!element.is0()) {
                this.point1 = point1 ;
                this.point2 = point1.addVector(element) ;
            } else {
                throw "A straight line can't be definied by a null vector" ; 
            } 
        }

    }

    equation() {
        /*  return the equation of the straight line on the form of and array od three members. Equation ax+by+c = 0
            is coded by [a,b,c] */
        if (this.point1.sameAbsciss(this.point2)) {
            return([1,0,-this.point1.x])
        } else {
            let direction = ( this.point1.y - this.point2.y) / (this.point1.x - this.point2.x)
            let ordinateOrigin = this.point1.y - direction*this.point1.x ;
            return([-direction, 1, -ordinateOrigin])
        }
    }

    isBefore(point1, point2) {
        /*  Indicate in fot two points point1 and point2 on the line this, point1 < point2. To answer a question, we use the 
            order relation defined by M <= N if M.x <= N.x and line is not parallel to y-axe, and M.y < N.y otherwise*/
        if (this.point1.sameAbsciss(this.point2)) {
            return point1.y <= point2.y ;
        } else {
            return point1.x <= point2.x ;
        }
    }

    isBeforeStrict(point1, point2) {
        /*  strict vection of isBefore */
        if (this.point1.sameAbsciss(this.point2)) {
            return point1.y < point2.y ;
        } else {
            return point1.x < point2.x ;
        }
    }

    distanceToPoint(point) {
        /*  Compute the distance between line this and point */
        let lineEquation = this.equation() ;
        let numerator = Math.abs(lineEquation[0]*point.x + lineEquation[1]*point.y + lineEquation[2] );
        let denominator = Math.sqrt(lineEquation[0] ** 2 + lineEquation[1] ** 2);
        return numerator/denominator;
    }

    orthogonalProjection(point) {
        /*  return the orthogonal prohjection of point on the straight line this */
        let direction = new vector(this.point1, this.point2) ;
        let vectorToProject = new vector(this.point1, point) ;
        return this.point1.addVector(direction.othogonalProjection(vectorToProject)) ;
    }

    getOrthogonalLine() {
        /*  return the orthogonal line of this which pass through the point (0,0) */
        let lineEquation = this.equation() ;
        let originPoint = new point(0,0) ;
        let orthogonalDirection = new vector(lineEquation[0], lineEquation[1]) ;
        let res = new straightLine(originPoint, orthogonalDirection) ;

        return res ;
    }

}

class polygon {
    constructor(vertices, parallelEdge = null) {
        /*  vertices array is a list of point which determine the vertices of the polygon. The list must contain 
            only one exemplary of points. 
            One can precise the list of the edge that are parallel on the polygon, it will be usefull to optimise 
            sat algorithme. The structure of parralelEdge is an array of array constaining the number of the edge, 
            with the first edge defining by the first two points, and then following the order. parallelEdge must be null
            or contains all the edge, not only part of the edge. Two array in the same 
            subarray arre parallels. For example, let ABCDEF a polygon, the the array [[0,2],[1,4],[3],[5]] mean that AB//CD, 
            BC//EF , DE and FA are parallel to no other edges. */
        this.vertices = vertices ;
        this.parallelEdge = parallelEdge;
    }

    translate(translationVector) {
        /*  tranlaction of the polygon of vector v. Methods don't return new polygon, but modify directly the attributes*/
        this.vertices.forEach(point => {
            point.translate(translationVector) ;
        })
    }

    /*  SAT algorithm */

    getOrthogonalLines() {
        /*  return the set of orthogonal lines of the edges of this. If two edges of the polygone are parallels
            we return only one line */
        let verticesLine = [];
        let nbVertices = this.vertices.length;
        let line ;
        let nbLoop ;
        if (this.parallelEdge === null) {
            nbLoop = nbVertices ;
            for (let k = 0; k < nbLoop; k++) {
                line = new straightLine(this.vertices[k], this.vertices[(k+1)%nbVertices]) ;
                verticesLine.push(line) ;
            }
        } else {
            nbLoop = this.parallelEdge.length ;
            let firstPoint ;
            for (let k = 0; k < nbLoop; k++) {
                firstPoint = this.parallelEdge[k][0]
                line = new straightLine(this.vertices[firstPoint], this.vertices[(firstPoint+1)%nbVertices]) ;
                verticesLine.push(line) ;
            }
        }

        let orthogonalLines = []
        verticesLine.forEach(line => {
            orthogonalLines.push(line.getOrthogonalLine()) ;
        })

        return orthogonalLines ;
    }

    sat(other, areParallel = false) {
        /*  Separating Axes Theorem (S. Gottschalk. Separating axis theorem. Technical Report TR96-024,Department
            of Computer Science, UNC Chapel Hill, 1996) : 
                Two convex polytopes are disjoint iff there exists a separating axis orthogonal 
                to a face of either polytope or orthogonal to an edge from each polytope. 
            The argument areParallet indicate if all the edge of other are parallel to an edge of this. In that
            case, there is no need to compute orthogonal lines of edges of other */

        // 1. Get all the orthogonal axis : 

        let orthogonalLines = this.getOrthogonalLines() ;
        if (!areParallel) {
            let otherOrthogonalLine = other.getOrthogonalLines() ;
            otherOrthogonalLine.forEach(line => {
                orthogonalLines.push(line)
            })
        }

        // 2. point projection on othogonal line

        let isSeparate = false ;
        let cpt = 0;
        
        while (cpt < orthogonalLines.length && !isSeparate) {
            // we test sepration on each orthogonal lines. We continue testing separtion until we find
            // a line wich separte the two polygones
            let line = orthogonalLines[cpt] ;

            let maxThis = new point(-Infinity,-Infinity) ;
            let maxOther = new point(-Infinity,-Infinity) ; 
            let minThis = new point(Infinity,Infinity) ;
            let minOther = new point(Infinity,Infinity) ; 

            let proj
            this.vertices.forEach(vertice => {
                // projection of this on the line. We only keep min and max projection which will
                // be compared to min and max projection of other
                proj = line.orthogonalProjection(vertice) ;
                if(line.isBefore(maxThis, proj)) {
                    maxThis = proj ;
                }

                if(line.isBefore(proj, minThis)) {
                    minThis = proj ;
                }
            })
            
            other.vertices.forEach(vertice => {
                // projection of other on the line. We only keep min and max projection which will
                // be compared to min and max projection of other
                proj = line.orthogonalProjection(vertice) ;
                if(line.isBefore(maxOther, proj)) {
                    maxOther = proj ;
                }

                if(line.isBefore(proj, minOther)) {
                    minOther = proj ;
                }
            })

            if (line.isBeforeStrict(maxThis, minOther) || line.isBeforeStrict(maxOther, minThis)) {
                // the projections are separated on the line if maxThis < minOhter or maxOther < minThis
                isSeparate = true
            }

            cpt++ ;
        }

    return isSeparate ;

    }
}

class square extends polygon {
    /*  A square extend polygon class. Square attributes are
        - a set of 4 points
        - a description of parallelEdge (see class polygon) which is always [[0,2], [1,3]]
        - a center : the barycenter of the square 
        - a direction which is polar coordinates of the first edge of the square.
        Two last attribute are commod in order to rotate the square according to its center.*/
    constructor(element1, element2) {
        /*  there is tw way to constructs a square : 
            - Given the coordinates of the two limit point of one of its edge. In this case, the other point 
            are construct in direct order : edge 2 direction is edge 1 direction rotate from pi/2
            - given its center and the polar coordinates of one of its edge. In this case, the other point 
            are construct in direct order : edge 2 direction is edge 1 direction rotate from pi/2. polar coordinate
            is an array of a positive number (the length of the edge), and an angle in randiant*/
        let point1 ;
        let point2 ;
        let point3 ;
        let point4 ;

        let direction ;
        let polarDirection ;
        let center ;

        if (element2 instanceof point) {
            point1 = element1 ;
            point2 = element2 ;
            direction = new vector(point1, point2);
            polarDirection = direction.polarCoordinate() ;
            point3 = point2.addVector(direction.orthogonalVector()) ;
            point4 = point3.addVector(direction.orthogonalVector().orthogonalVector()) ;
            center = point.segmentCenter(point1, point3) ;
        } else {
            polarDirection = element2 ;
            direction = new vector(polarDirection[0]*Math.cos(polarDirection[1]), polarDirection[0]*Math.sin(polarDirection[1])) ;
            // create a first square with good direction, and first point = (0,0)
            point1 = new point(0,0) ;
            point2 = point1.addVector(direction) ;
            point3 = point2.addVector(direction.orthogonalVector()) ;
            point4 = point3.addVector(direction.orthogonalVector().orthogonalVector()) ;

            // translate the square to the good position
            let initialCenter = point.segmentCenter(point1, point3) ;
            center = element1 ;
            let translationVector = new vector(initialCenter, center) ;

            point1.translate(translationVector) ;
            point2.translate(translationVector) ;
            point3.translate(translationVector) ;
            point4.translate(translationVector) ;
        }
        super([point1, point2, point3, point4], [[0,2], [1,3]])
        this.center = center ;
        this.polarDirection = polarDirection ;
    }

    rotate(angle) {
        /* rotate the square according to its center. angle is in radiant */
        this.polarDirection[1] = angle ;
        let direction = new vector(this.polarDirection[0]*Math.cos(this.polarDirection[1]), 
            this.polarDirection[0]*Math.sin(this.polarDirection[1])) ;

        this.vertices[0] = new point(0,0) ;
        this.vertices[1] = this.vertices[0].addVector(direction) ;
        this.vertices[2] = this.vertices[1].addVector(direction.orthogonalVector()) ;
        this.vertices[3] = this.vertices[2].addVector(direction.orthogonalVector().orthogonalVector()) ;

        let initialCenter = point.segmentCenter(this.vertices[0], this.vertices[2]) ;
        let translationVector = new vector(initialCenter, this.center) ;

        this.vertices[0].translate(translationVector) ;
        this.vertices[1].translate(translationVector) ;
        this.vertices[2].translate(translationVector) ;
        this.vertices[3].translate(translationVector) ;
    }

    translate(transactionVector) {
        // extend class translate of polygon by translative square center plus the points or the polygon
        super.translate(transactionVector) ;
        this.center.translate(transactionVector) ;
    }

    getLowestPointIndex() {
        // indicate the index of the point(s) of minimal ordinates of the the square. If two points, we return
        // first the point with lowest abscissa.    
        let lowestPoint = new point(Infinity, Infinity) ;
        let lowestPointIndex = null ;

        for (let k = 0; k < this.vertices.length; k++) {
            if (keepNDec(this.vertices[k].y,6) < keepNDec(lowestPoint.y,6)) {
                // comparision are mad with 6 decimal to avoid precision error of java script
                lowestPoint = this.vertices[k] ;
                lowestPointIndex = k;
            }
        }

        let res = [] ;
        for (let k = 0; k < this.vertices.length; k++) {
            if (keepNDec(lowestPoint.y,6) === keepNDec(this.vertices[k].y,6)) {
                // comparision are mad with 6 decimal to avoid precision error of java script
                console.log(parseFloat(lowestPoint.y).toPrecision(2))
                res.push(k)
            }
        }

        if(res.length == 2) {
            // return lowest point by abscissa order
            if (this.vertices[res[0]].x > this.vertices[res[1]].x) {
                res = [res[1], res[0]] ;
            };
        } 

        return res ;
    }

}


/* Games elements classes */

class hero extends square {
    constructor() {
        /*  Body hitBox */
        let intialPosition = new point(3/2,1/2) ;
        super(intialPosition, [1,Math.PI/4])

        /*  foot hitBox */
        let footPoint = this.getLowestPointIndex()
        if (footPoint.length == 2) {
            let footPoint1 = this.vertices[footPoint[0]] ;
            let footPoint2 = this.vertices[footPoint[1]] ;
            let footPoint3 = point.segmentCenter(footPoint2, this.vertices[positivMod(footPoint[1]+1,4)]) ;
            let footPoint4 = point.segmentCenter(footPoint1, this.vertices[positivMod(footPoint[0]-1,4)]) ;

            this.foot = new polygon([footPoint1,footPoint2, footPoint3, footPoint4], [[0,2], [1,3]]);
        } else {
            let footPoint1 = this.vertices[footPoint[0]] ;
            let footPoint2 = point.segmentCenter(footPoint1, this.vertices[positivMod(footPoint[0]+1,4)]) ;
            this.foot = new square(footPoint1, footPoint2) ;
        }
    }
}


// main

const canvasGame = document.getElementById("canvas");
const ctx = canvasGame.getContext("2d");

let canvasWidth = document.getElementById("game-interface").offsetWidth

ctx.canvas.width  = canvasWidth/2;
ctx.canvas.height = canvasWidth/4;


