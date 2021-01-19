/*  Various functions

    This part gather various functions which will be usefull in the rest of the code
*/

function positivMod(n, mod) {
    // classic n%mod operator gives a negative number when n < 0. This function give a positive modulo in such cases .
    return (n%mod+mod)%mod;
}
/*  Module polygon contains some classes to manage polygon gemetry in code. Principal classes is 
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

}

class segment {
    constructor(point1, point2) {
        this.point1 = point1 ;
        this.point2 = point2 ;
    }

    center() {
        /* return the center of the segment */
        let res = new point((this.point1.x+this.point2.x)/2, (this.point1.y+this.point2.y)/2) ;
        return res ;
    }

    containPoint(point) {
        /*  for a point which lie on line (point1, point2) (be carreful, the methode doesn't verify it), 
            say if point belongs to the segment. */
        let vector1 = new vector(this.point1, this.point2) ;
        let vector2 = new vector(this.point1,point) ;

        let scalarProd1 = vector1.scalarProduct(vector2) ;
        let scalarProd2 = vector1.scalarProduct(vector1) ;

        return scalarProd1 >= 0 && scalarProd1 <= scalarProd2 ;
    }

    intersect(other) {
        /*  other is antoher segment of the line with direction this. Return true if intersections between 
            this and other is not null*/

        let segmentVector = new vector(this.point1, this.point2) ;
        let vector1 = new vector(this.point1, other.point1) ;
        let vector2 = new vector(this.point1, other.point2) ;
        let vector3 = new vector(this.point2, other.point1) ;
        let vector4 = new vector(this.point2, other.point2) ;

        let scalarProd1 = segmentVector.scalarProduct(vector1) ;
        let scalarProd2 = segmentVector.scalarProduct(vector2) ;
        let scalarProd3 = segmentVector.scalarProduct(vector3) ;
        let scalarProd4 = segmentVector.scalarProduct(vector4) ;

        return !( (scalarProd1 < 0 && scalarProd2 < 0) || (scalarProd3 > 0 && scalarProd4 > 0)) ;

    }

}

class polygon {
        /*  gather some methods relatives to polygones. Polygones are given by an array of en points (their vertices).
        Class polygon accept semgement as degenrate polygon. It can usefull to define segment as polygon to apply them
        sat algorithm. On the other hand, class polygon doesn't herit from segment, so segment method can't be applied
        directly to "segment polygon" */
    constructor(vertices) {
        /*  vertices array is a list of point which determine the vertices of the polygon. The list must contain 
            only one exemplary of points. 
         */
        this.vertices = vertices ;
    }

    translate(translationVector) {
        /*  tranlaction of the polygon of vector v. Methods don't return new polygon, but modify directly the attributes*/
        this.vertices.forEach(point => {
            point.translate(translationVector) ;
        })
    }

    edges() {
        /* return the list of edges of the polygon */
        let edges = [] 
        let nbVertices = this.vertices.length
        if (nbVertices > 2) {
            for (let k = 0; k < nbVertices; k++) {
                let edge = new segment(this.vertices[k], this.vertices[(k+1)%nbVertices]) ;
                edges.push(edge) ;
            }
        } else {
            let edge = new segment(this.vertices[0], this.vertices[1]) ;
            edges.push(edge) ;
        }
        return(edges)
    }

    isoBarycenter() {
        /* return ths isoBarycenter (with coefficients = 1) of the polygon */
        let barycenterAbscissa = 0 ;
        let barycenterOrdinate = 0 ;
        this.vertices.forEach(point => {
            barycenterAbscissa += point.x;
            barycenterOrdinate += point.y;
        }) 

        let res = new point(1/this.vertices.length * barycenterAbscissa, 1/this.vertices.length * barycenterOrdinate) ;
        return res ;
    }

    /*  SAT algorithm */

    static separation(other, edge, barycenter) {
        /*  Considers a polygone with barycenter "barycenter", and "other" another polygone. Suppose
            that edge is an edge of the first polygon. Separation methode return true if the straightLine with direction
            edge separte the two polygones. To do that : 
            - we consider the equation of the line generate by edge : ax+by+c = 0
            - we replace x and y by coordinate of the barycenter and get a number A
            - for all points of other, we replace x and y by point coordinates and get numbers B1,... Bn
            
            - If A = 0, it means that the first polygone is a segment. In that case, edge separate other if
            all B1,...,Bn have the same sign.
            - if A != 0, the edge seprate the polygones if B1,...,Bn have same sign, and that sign is different of A
            
            In some situations, we can have one or two of the B1,...,Bn which are = 0. Never more, because, it would say 
            that three points of other are align. In that case, if Bk beloongs to edge, or the segment [Bk, Bl] have
            non null intersection with edge, one have found common points for the two polygons, and they are not separate*/

        let otherNbVertices = other.vertices.length ;
        let segmentLine = new straightLine(edge.point1, edge.point2) ;
        let equation = segmentLine.equation()

        let thisSide = equation[0]*barycenter.x + equation[1]*barycenter.y + equation[2] ;

        let pointSideSet = [];
        let pointSide = [];

        let pointOnSepartor = [];

        for (let k = 0; k < otherNbVertices ; k++) {
            pointSide = equation[0]*other.vertices[k].x + equation[1]*other.vertices[k].y + equation[2] ;
            pointSideSet.push(pointSide) ;
            if( keepNDec(pointSide, 10) === 0) {
                pointOnSepartor.push(other.vertices[k])  
            }
        }

        let commonPoint = false ;
        if (pointOnSepartor.length == 1) {
            if (edge.containPoint(pointOnSepartor[0])) {
                commonPoint = true ;
            }
        } else if (pointOnSepartor.length == 2) {
            let alignSegment = new segment(pointOnSepartor[0], pointOnSepartor[1]);
            if (edge.intersect(alignSegment)) {
                commonPoint = true ;
            }
        }

        if (commonPoint) {
            return false
        } else  {
            let minPointSide = Math.min.apply(Math, pointSideSet) ;
            let maxPointSide = Math.max.apply(Math, pointSideSet)
            
            if (keepNDec(thisSide, 10) == 0) {
                return keepNDec(minPointSide, 10)* keepNDec(maxPointSide, 10) >= 0 ;
            } else {
                return keepNDec(thisSide, 10)* keepNDec(maxPointSide, 10) <= 0 && 
                            keepNDec(minPointSide, 10)* keepNDec(thisSide, 10) <= 0 ;
                ;
    
            }
        }
    }

    sat(other) {
    /*  Separating Axes Theorem (S. Gottschalk. Separating axis theorem. Technical Report TR96-024,Department
        of Computer Science, UNC Chapel Hill, 1996) : 
            Two convex polytopes are disjoint iff there exists a separating axis orthogonal 
            to a face of either polytope or orthogonal to an edge from each polytope.
            
        Our version of sat can also sepate segments which are degenerate polygons.
     */

        let thisEdges = this.edges() ;
        let otherEdges = other.edges() ;

        let thisBarycenter = this.isoBarycenter() ;
        let otherBarycenter = other.isoBarycenter() ;

        let isSeparated = false ;
        let cpt = 0 ;
        do {
            /* try to find a separator wicth edge of this */
            isSeparated = polygon.separation(other, thisEdges[cpt], thisBarycenter) ;
            cpt ++;
        } while (cpt < thisEdges.length & !isSeparated) 

        if (!isSeparated) {
            /* if no edges of this ae separting, one try with edges of other */
            cpt = 0;
            do {
                isSeparated = polygon.separation(this, otherEdges[cpt], otherBarycenter) ;
                cpt ++;
            } while (cpt < otherEdges.length & !isSeparated) 
        }

        return isSeparated ;
    }

}

class square extends polygon {
    /*  A square extend polygon class. Square attributes are
        - a set of 4 points
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
        let diagonal ;
        let center ;

        if (element2 instanceof point) {
            point1 = element1 ;
            point2 = element2 ;
            direction = new vector(point1, point2);
            polarDirection = direction.polarCoordinate() ;
            point3 = point2.addVector(direction.orthogonalVector()) ;
            point4 = point3.addVector(direction.orthogonalVector().orthogonalVector()) ;
            diagonal = new segment(point1, point3)
            center = diagonal.center() ;
        } else {
            polarDirection = element2 ;
            direction = new vector(polarDirection[0]*Math.cos(polarDirection[1]), polarDirection[0]*Math.sin(polarDirection[1])) ;
            // create a first square with good direction, and first point = (0,0)
            point1 = new point(0,0) ;
            point2 = point1.addVector(direction) ;
            point3 = point2.addVector(direction.orthogonalVector()) ;
            point4 = point3.addVector(direction.orthogonalVector().orthogonalVector()) ;

            // translate the square to the good position
            diagonal = new segment(point1, point3)
            let initialCenter = diagonal.center() ;
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

        let diagonal = new segment(this.vertices[0], this.vertices[2]) ;
        let initialCenter = diagonal.center() ;
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

class hero {
    constructor() {
        /*  a hero is the set of a body which is a square, and foot which will be usefull to test if the hearo land
            on a block (the foot touch the rough of the block) or not. foot is and array of segments declared a 
            polygon to apply them sat algorithm. If the square is horizontal, the array contains only one segment : the 
            lowest, else it contains two segments : those around the lowest point */
        /*  Body hitBox */
        let intialPosition = new point(3/2,3/2) ;
        this.body = new square(intialPosition,[1,0] );

        /*  foot hitBox */
        let footPoint = this.body.getLowestPointIndex()
            // the default body have angle 0, so footpoint return always 2 points
        let footPoint1 = this.body.vertices[footPoint[0]] ;
        let footPoint2 = this.body.vertices[footPoint[1]] ;
        let foot1 = new polygon([footPoint1,footPoint2]);
        this.foot = [foot1]
    }

    rotate(angle) {
        this.body.rotate(angle) ;

        let footPoint = this.body.getLowestPointIndex()
        if (footPoint.length == 2) {
            let footPoint1 = this.body.vertices[footPoint[0]] ;
            let footPoint2 = this.body.vertices[footPoint[1]] ;
            let foot1 = new polygon([footPoint1,footPoint2]);
            this.foot = [foot1]
        } else {
            let footPoint1 = this.body.vertices[footPoint[0]] ;
            let footPoint2 = this.body.vertices[positivMod(footPoint[0]+1,4)] ;
            let footPoint3 = this.body.vertices[positivMod(footPoint[0]-1,4)] ;
            let foot1 = new polygon([footPoint1, footPoint2]) ;
            let foot2 = new polygon([footPoint1, footPoint3]) ;
            this.foot = [foot1, foot2] ;
        }

    }

}

class drawingGrid {
    constructor() {
        let canvas = document.getElementById("canvas");
        this.ctx = canvas.getContext("2d") ;
        this.width = document.getElementById("game-interface").offsetWidth;
        this.height = this.width/2 ;
        this.unity = this.width/40;
    }

    drawGrid() {
        this.ctx.canvas.width  = this.width;
        this.ctx.canvas.height = this.height;
    }

    gridAbscissa(x) {
        /*  scalme abscissa to this.unity */
        return x*this.unity ;
    }

    gridOrdinate(y) {
        /*  in canvas, "ordinate 0" is top, and not bottom as usual in math. This methode transcript "math" 
            ordinate in "canvas" ordinate. Moreover, it scale to this.unity */
        return this.height-y*this.unity ;
    }

    plotHero(heroInstance) {
        let heroBody = new Path2D() ;
        heroBody.moveTo(this.gridAbscissa(heroInstance.body.vertices[0].x), this.gridOrdinate(heroInstance.body.vertices[0].y));
        heroBody.lineTo(this.gridAbscissa(heroInstance.body.vertices[1].x), this.gridOrdinate(heroInstance.body.vertices[1].y));
        heroBody.lineTo(this.gridAbscissa(heroInstance.body.vertices[2].x), this.gridOrdinate(heroInstance.body.vertices[2].y));
        heroBody.lineTo(this.gridAbscissa(heroInstance.body.vertices[3].x), this.gridOrdinate(heroInstance.body.vertices[3].y));
        heroBody.closePath();

        let heroFoot = new Path2D() ;
        heroFoot.moveTo(this.gridAbscissa(heroInstance.foot[0].vertices[0].x),this.gridOrdinate(heroInstance.foot[0].vertices[0].y) )
        heroFoot.lineTo(this.gridAbscissa(heroInstance.foot[0].vertices[1].x),this.gridOrdinate(heroInstance.foot[0].vertices[1].y) )
        if (heroInstance.foot.length === 2) {
            heroFoot.moveTo(this.gridAbscissa(heroInstance.foot[1].vertices[0].x),this.gridOrdinate(heroInstance.foot[1].vertices[0].y) )
            heroFoot.lineTo(this.gridAbscissa(heroInstance.foot[1].vertices[1].x),this.gridOrdinate(heroInstance.foot[1].vertices[1].y) )    
        } 
        console.log(this.gridAbscissa(heroInstance.foot[0].vertices[1].x)) ;
        console.log(this.gridOrdinate(heroInstance.foot[0].vertices[1].y)) ;

        this.ctx.fillStyle = "red" ;
        this.ctx.strokeStyle = "green" ;
        this.ctx.fill(heroBody) ;
        this.ctx.lineWidth = 3 ;
        this.ctx.stroke(heroFoot) ;
    }
}


// main

/*const canvasGame = document.getElementById("canvas");
const ctx = canvasGame.getContext("2d");

const canvasWidth = document.getElementById("game-interface").offsetWidth
const canvasDimension = {unity: canvasWidth/40, width: canvasWidth, height: canvasWidth/2 }

ctx.canvas.width  = canvasDimension.width;
ctx.canvas.height = canvasDimension.height;*/

const grid = new drawingGrid();
const heroInstance = new hero();
grid.drawGrid() ;
heroInstance.rotate(3*Math.PI/4)
grid.plotHero(heroInstance) ;


/*toto = new hero() ;

console.log(toto.body) ;
console.log(toto.foot[0]) ;
console.log(toto.foot[1]) ;

toto.rotate(Math.PI/4) ;
console.log(toto.body) ;
console.log(toto.foot[0]) ;
console.log(toto.foot[1]) ;*/







