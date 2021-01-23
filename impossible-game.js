/*  Various functions

    This part gather various functions which will be usefull in the rest of the code
*/

function positivMod(n, mod) {
    // classic n%mod operator gives a negative number when n < 0. This function give a positive modulo in such cases .
    return (n%mod+mod)%mod;
}


function keepNDec(x, nDec) {
    /*  "round" a number to its nDec decimal */

    return Math.round(x*10 ** nDec)/(10 ** nDec);
}

/*  Module polygon : come from a js module. I don't know actually how to use them, but in the future, this part will 
    be seprate from the main code

    Contains some classes to manage polygon gemetry in code. Principal classes is 
    classes polygon which contains a method sat for (Separating Axes Theorem) which compute if two polygones are separated.

    Other classes are : 
    - point, vector, straightLine : they are requisite in polygon definitions and methods
    - square which extend polygon in the special case of square
*/

class point {
    constructor(x,y) {
        /* Point are coded ny their two coordinates */
        this.x = x ;
        this.y = y ;
    }

    copy() {
        return new point(this.x, this.y) ;
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
        this.x = keepNDec(this.x, 10) ;
        this.y += v.y;
        this.y = keepNDec(this.y, 10) ;
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

    containPoint(point) {
        let equationLine = this.equation() ;
        return keepNDec(equationLine[0]*point.x + equationLine[1]*point.y + equationLine[2], 10) === 0 ;
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

    /*  SAT algorithm. Be carrefull : work only for convex polygons*/

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

        Be carrefull : work only for convex polygons.
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
        this.polarDirection[1] += angle ;
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

/*  Games elements classes 

    I have organize the game in three principal elements : 
     - the heros element : on class which contains all the method of hero
     - the grid element which contain all elements about the level. There is more than one class : some small class
       for element as peak or platform, and a class grid
     - the drawing element which allow to draw the game */


    /* hero */

class hero {
    constructor() {
        /*  a hero is the set of a body which is a square, and foot which will be usefull to test if the hero land
            on a block (the foot touch the roof of the block) or not. Soot is and array of segments declared a 
            polygon to apply them sat algorithm. If the square is horizontal, the array contains only one segment : the 
            lowest, else it contains two segments : those around the lowest point */

        /*  Body hitBox */
        let heroCenterCoordinate = [5/2,3/2]
        let intialPosition = new point(heroCenterCoordinate[0],heroCenterCoordinate[1]) ;
        this.body = new square(intialPosition,[1,0] );

        /*  foot hitBox */
        let footPoint = this.body.getLowestPointIndex()
            // the default body have angle 0, so footpoint return always 2 points
        let footPoint1 = this.body.vertices[footPoint[0]].copy() ;
        let footPoint2 = this.body.vertices[footPoint[1]].copy() ;
        let foot1 = new polygon([footPoint1,footPoint2]);
        this.foot = [foot1]

        /*  Physical attributes 
            Those parameters are computed in order the hero do a jump of 4 unity heigh and 3 unity long. 
            We use classical newtonian physic for the trajectories which says that he gravity center of the hero 
            follow the next trajectory (with t=0 as begning of a jump)
                x(t) = vx * t
                y(t) = -1/2*g*t^2 + vy0 * t + y0
            where 
                - vx is the horizontal speed or the hero : in the game it is a constant, so no need to take vx(0)
                - g id the gravitional constant
                - vy0 is the initial vertical speed when jump.
                - y0 is the initial y coordinate of the center of the hero
            This formula work equaly when the hero fall from a bloc (in that case vz0=0) and when the hero is on a 
            bloc (in that case, g and vy0 = 0).
            
            On the game, vx is fixed, and we choose g and vz0 so that a jump is xJump long and yJump height. 
            xJump and zJump are fixed too
        */

        this.vx = 10 ; // horizontal speed of the hero
        this.vy0 = 0 ; // vertical speed. when jump :  (2*this.zJump)/((this.xJump/(2*this.vx))
        this.xJump = 4 ; // length of a jump
        this.yJump = 3 ; // height of a jump
        this.g = 0 ; // gravitional constant. value when jump or fall from a platform : (2*this.zJump)/((this.xJump/(2*this.vx))**2) ;
        this.t = 0 ; // counter of time when the hero begin jumping or falling in order to use equation. must be 0 at the begning of a jump or fall
        this.isJumping = false ; 
            // use to avoid the gamer to do multiple jump, when jump, it become true and become false only when the hero land 
    
        /*  hero status */

        this.hasStarted = false ;
        this.isDead = false ; // for game over
        
    }

    rotate(angle) {
        /* rotate the hero : body + foot */
        this.body.rotate(angle) ;

        let footPoint = this.body.getLowestPointIndex()
        if (footPoint.length == 2) {
            let footPoint1 = this.body.vertices[footPoint[0]].copy() ;
            let footPoint2 = this.body.vertices[footPoint[1]].copy() ;
            let foot1 = new polygon([footPoint1,footPoint2]);
            this.foot = [foot1]
        } else {
            let footPoint1 = this.body.vertices[footPoint[0]].copy() ;
            let footPoint2 = this.body.vertices[positivMod(footPoint[0]+1,4)].copy() ;
            let footPoint3 = this.body.vertices[positivMod(footPoint[0]-1,4)].copy() ;
            let foot1 = new polygon([footPoint1, footPoint2]) ;
            let foot2 = new polygon([footPoint1.copy(), footPoint3]) ;
            this.foot = [foot1, foot2] ;
        } 

    }

    translate(transactionVector) {
        /* translate the hero : body + foot */
        this.body.translate(transactionVector) ;
        this.foot.forEach(footValue => {footValue.translate(transactionVector)}) ;
    }

    footContactWithRoof(previousFoot, platformInstance) {
        /*  Compute collision to know if the hero has land. The hero land on a platform if its foot have collision
            with platform roof (see class plateform). However, because the game is not continuous but dicrete, and 
            that roof and foot are lines, it could arive that at t, foot is above a roof and a t+1, it is below. In that
            case, the game will consider there is no collision. To avoid that. We construct a footPolygon which is the 
            polygon obtain by contatenate previous foot position and next foot position. Thank to that, if
            footPolygon have collision wich the roof, it means that bettween t and t+1, the hero have land on the roof. */
        let cpt = 0 ;

        let footPolygon ;
        for (let k = 0; k < this.foot.length; k++) {
            let lineTest = new straightLine(this.foot[k].vertices[1], previousFoot[k][0]) ;
            /*  Test if roof polygon is flat. In that case, footPolygon is a segment, and in order sat method of polygon work
                one need to keep only two points. sat work only if there is no more than two points align in polygon*/
            if (lineTest.containPoint(previousFoot[k][1])) {
                footPolygon = new polygon([this.foot[k].vertices[1],previousFoot[k][0]]) ;
                /*  only arrive when the hero have rectiling movment : when it is horizontal and placed on the ground.
                    In that case, foot have only one element and first point of the polygon is the leftest. So the segment
                    as defined is the larger possible segment*/
            } else {
                footPolygon = new polygon([this.foot[k].vertices[0], this.foot[k].vertices[1],
                    previousFoot[k][1], previousFoot[k][0]]) ;
            }
            /*  the order of point is important. When moving, each point is translate from the same vector, so we need 
                to choose order of the point to not have cross polygon */
            if (!footPolygon.sat(platformInstance.roof)) {
                /* if one collision : it means that the hero land, we conslude verifying cpt > 0 */
                cpt ++ ;
            }
        }

        return cpt > 0 ;
    }

    move(drawingInstance, gridInstance) {
        /*  The gravity center of the hero follow the next trajectory (with t=0 as begning of a jump)
                x(t) = vx * t
                y(t) = -1/2*g*t^2 + vy0 * t + y0
            where 
                - vx is the horizontal speed or the hero : in the game it is a constant, so no need to take vx(0)
                - g id the gravitaional constant
                - vy0 is the initial vertical speed when jump.
                - y0 is the initial y coordinate of the center of the hero
            This formula work equaly when the hero fall from a bloc (in that case vz0=0) and when the heri is on a bloc (in that case, g and vy0 = 0).

            On this methode, we consider finite diffrecne of this equation : y(t+dt) - y(t) and x(t+dt) - x(t). Because of the
            quadratiq nature of equation 2, the time t still appears in equation for y, and so we can't keep only dt value, we 
            need to use t.
            
            On the game, vx is fixed, and we choose g and vz0 so that a jump is xJump long and yJump height. xJump and zJump are fixed too
            */

        let previousFoot = [] ;
        this.foot.forEach(foot => {
            previousFoot.push([foot.vertices[0].copy(), foot.vertices[1].copy()]) ;
        })
            // save previous foot to verify if the hero land (see method footContactWithRoof)

        let dt = 1/drawingInstance.fps ;
            // time interval computed according to the game fps

        let translationVector = new vector(this.vx * dt, -1/2*this.g * dt * (2*this.t+dt) + this.vy0 * dt ) ;
        this.translate(translationVector) ;
            // translate the hero with finite diffrence equation

        this.t += dt ;
            /*  for next step, t become t+dt. If the movement is not a jump or a fall, t will be set to 0 next. Indeed
                t is not mandatory for recticlign horizontal movement, and that way, t will be 0 at the begning of a jump 
                or a fall as expected.*/

        /*  Grid interaction check : We test collision of the hero with its neihbourg envirenement. 
             - collision with a peak = dead
             - collision with a platform = dead if not collisions with the floor of the platform
        */

        let nbContact = 0 ;
            // at the end  of the collsion check, if not dead, use to know if there is at leat one contact (necessary)
            // a foot-roof contact, a modify physical constant in consequence
        let floorContactCoordinate = [] ;
            // if landing, use to replace the hero to good coordinate. Indeed, because of the non consinuity, the hero could
            // be below a roof, floorContactCoordinate contains the coOrdinates of the platform where the hero land,
            // and this way, we can corrEct y position of the hero

        let aroundGrid = gridInstance.grid.slice(Math.max(Math.floor(this.body.center.x-1),0), Math.floor(this.body.center.x+2)) ;
        // neighbour elements grid of the hero
        aroundGrid.forEach(col => {
            if (col != undefined) {
                col.forEach(element => {
                    if (element instanceof platform) {
                        if (!this.body.sat(element.platform)) {
                            nbContact ++ ;
                            if (this.footContactWithRoof(previousFoot,element)) {
                                floorContactCoordinate.push(element.platform.center) ;
                            } else {
                                this.isDead = true ;
                            }                        }
                    } else if (element instanceof peak) {
                        if (!this.body.sat(element.peak)) {
                            nbContact ++ ;
                            this.isDead = true ;
                        }
                    }
                })
            }
        })

        if (this.body.center.y < 0.5) {
            /*  The lower set of roof have coordinate 1, if the square fall under it, it means that it falls in a hole */
            this.isDead = true ;
        }

        if (!this.isDead) {
            if (nbContact > 0) {
                let newCenter = new point(this.body.center.x, floorContactCoordinate[0].y+1) ;
                let translateVector = new vector(this.body.center, newCenter) ;
                this.translate(translateVector) ;

                this.rotate(2*Math.PI - this.body.polarDirection[1]) ;

                this.g = 0 ; // g is compensated by newton 3rd law
                this.vy0 = 0 ; // the vertical movment stop
                this.t = 0 ; // t is not mandatory in equation anymore, and will be 0 at the begning of a jump or a fall
                this.isJumping = false ; // the gamer can jump
            } else {
                this.g = (2*this.yJump)/((this.xJump/(2*this.vx))**2) ; // no compensation by newton 3rd law
                this.rotate(-Math.PI/(drawingInstance.fps * this.xJump/(2*this.vx))) ; // to look pretty : the hero rotate when not on a roof
                this.isJumping = true ; // can't jump before the end of the jump /fall
            }
        }
    }

    jump() {
        /* Modify physical constant in order to the next move (methode move) is a jump */
        if(!this.isJumping) {
            this.isJumping = true ; // during a jump, the hero can't jump anymore
            // can't jump if already jumping
            this.g = (2*this.yJump)/((this.xJump/(2*this.vx))**2) ; // no compensation by newton 3rd law
            this.vy0 = (2*this.yJump)/(this.xJump/(2*this.vx)) ; // to look pretty : the hero rotate when not on a roof
            this.t = 0; // at the begning of a jump, t must be = 0 in order to the equation are ok
        }
    }

}

    /*  Grid elements : 
        
        Two kind of grid elements : 
        - platform on which one can land, or crash if foot hero don't touch platform roof
        - peak on which one can die*/

class platform {
    /*  A platform is 
         -  a square of edge length 1 and angle 0. 
         -  the the roof which is the upper edge of the square. It is used to verify if the hero land on the platforme : 
            the foot and the roof enter in collision, or not (in that case, if collision, it's game over)
            
        For each element one add col = floor(x-positon) attribute. Elements will be organized on a grid which is an array. In array
        cell n, we will place all element which col = n. That way, it is easy to get all element of the neibourhood
        of the hero to test collisions*/
    constructor(x, y) {
        /* x,y are the coordinates of the center of the platform */
        this.col = Math.floor(x) ;
        let platformCenter = new point(x, y) ;
        this.platform = new square(platformCenter,[1,0]) ;
        this.roof = new polygon([this.platform.vertices[2], this.platform.vertices[3]])
    }
}

class peak {
        /*  peak is a triangle. There is 4 kind of peak according to there orientation.
            
        For each element one add col = floor(x-positon) attribute. Elements will be organized on a grid which is an array. In array
        cell n, we will place all element which col = n. That way, it is easy to get all element of the neibourhood
        of the hero to test collisions*/
    constructor(x, y, orientation) {
        /* x, y are the coordinates of the lowest leftest point of the square in which the triangle is inscribed */
        let point1, point2, point3 ;
        switch(orientation) {
            case 'up' :
                point1 = new point(x, y) ;
                point2 = new point(x+1, y) ;
                point3 = new point(x+1/2, y+1) ;
                break ;
            case 'down' :
                point1 = new point(x, y+1) ;
                point2 = new point(x+1, y+1) ;
                point3 = new point(x+1/2, y) ;
                break ;
            case 'left' :
                point1 = new point(x+1, y+1) ;
                point2 = new point(x+1, y) ;
                point3 = new point(x, y+1/2) ;
                break ;
            case 'right' : 
                point1 = new point(x, y+1) ;
                point2 = new point(x, y) ;
                point3 = new point(x+1, y+1/2) ;
                break ;

        }
        this.col = Math.floor(x) ;
        this.peak = new polygon([point1, point2, point3])
    }
}

class grid {
    /*  grid is a discretisation of the level. All element have a col value which is the floor of there x-position.
         Those element will be organized on the grid which is an array. In array cell n, we will place all element
        which col = n. That way, it is easy to get all element of the neibourhood
        of the hero to test collisions*/
     
    constructor() {
        /* by default a grid is empty, we add element using methodes */
        this.grid = [] ;
    }

    addPlatform(platformInstance) {
        /*  Given a platform, place a platform in the grid */
        if (this.grid[platformInstance.col] != undefined) {
            this.grid[platformInstance.col].push(platformInstance) ;
        } else {
            this.grid[platformInstance.col] = [platformInstance] ;
        }
    } 

    addPeak(peakInstance) {
        /*  Given a peak, place a platform in the grid */
        if (this.grid[peakInstance.col] != undefined) {
            this.grid[peakInstance.col].push(peakInstance) ;
        } else {
            this.grid[peakInstance.col] = [peakInstance] ;
        }
    }


    removeCol(start, end) {
        /* remove the entire element of grid from col start to col end exclude (begning by 0) */
        for (let k = start; k < end; k++) {
            this.grid[k] = undefined ;
        }
    }

    defaultGrid(size) {
        /* create a default frid with a default ground of platform on first row. */
        this.grid = [] ;
        let platformInstance ;
        for (let k = 0 ; k < size ; k++) {
            platformInstance = new platform(k+1/2,1/2) ;
            this.addPlatform(platformInstance) ;
        }
    }
}

    /*  drawing class : gather all element used to draw each frame of the game */

class drawing {
    constructor() {
        let canvas = document.getElementById("canvas");
        this.ctx = canvas.getContext("2d") ;
        this.width = document.getElementById("game-interface").offsetWidth;
        this.height = this.width/3 ;
        this.unity = this.width/40;
        this.fps = 60 ;
    }

    drawGrid() {
        this.ctx.canvas.width  = this.width;
        this.ctx.canvas.height = this.height;
    }

    gridAbscissa(x) {
        /*  scale abscissa to this.unity */
        return x*this.unity ;
    }

    gridOrdinate(y) {
        /*  in canvas, "ordinate 0" is top, and not bottom as usual in math. This methode transcript "math" 
            ordinate in "canvas" ordinate. Moreover, it scale to this.unity */
        return this.height-y*this.unity ;
    }

    plotHero(heroInstance) {
        let heroCenterXPosition = heroInstance.body.center.x - 4 ;
        let heroBody = new Path2D() ;
        heroBody.moveTo(this.gridAbscissa(heroInstance.body.vertices[0].x - heroCenterXPosition), this.gridOrdinate(heroInstance.body.vertices[0].y));
        heroBody.lineTo(this.gridAbscissa(heroInstance.body.vertices[1].x - heroCenterXPosition), this.gridOrdinate(heroInstance.body.vertices[1].y));
        heroBody.lineTo(this.gridAbscissa(heroInstance.body.vertices[2].x - heroCenterXPosition), this.gridOrdinate(heroInstance.body.vertices[2].y));
        heroBody.lineTo(this.gridAbscissa(heroInstance.body.vertices[3].x - heroCenterXPosition), this.gridOrdinate(heroInstance.body.vertices[3].y));
        heroBody.closePath();

        /*let heroFoot = new Path2D() ;
        heroFoot.moveTo(this.gridAbscissa(heroInstance.foot[0].vertices[0].x),this.gridOrdinate(heroInstance.foot[0].vertices[0].y) )
        heroFoot.lineTo(this.gridAbscissa(heroInstance.foot[0].vertices[1].x),this.gridOrdinate(heroInstance.foot[0].vertices[1].y) )
        if (heroInstance.foot.length === 2) {
            heroFoot.moveTo(this.gridAbscissa(heroInstance.foot[1].vertices[0].x),this.gridOrdinate(heroInstance.foot[1].vertices[0].y) )
            heroFoot.lineTo(this.gridAbscissa(heroInstance.foot[1].vertices[1].x),this.gridOrdinate(heroInstance.foot[1].vertices[1].y) )    
        } */

        this.ctx.fillStyle = "red" ;
        /*this.ctx.strokeStyle = "green" ;*/
        this.ctx.fill(heroBody) ;
        /*this.ctx.stroke(heroFoot) ;*/
    }

    plotGrid(gridInstance, heroInstance) {
        let heroCenterXPosition = heroInstance.body.center.x - 4 ;
        let minColToKeep = Math.floor(Math.max(heroCenterXPosition, 0));
        let maxColToKeep = Math.floor(heroInstance.body.center.x+40)

        let platformDraw = new Path2D() ;
        let peakDraw = new Path2D() ;

        gridInstance.grid.slice(minColToKeep, maxColToKeep).forEach(gridRow => {
            if(gridRow != undefined) {
                gridRow.forEach(element => {
                    if (element instanceof platform) {
                        platformDraw.rect(this.gridAbscissa(element.platform.vertices[3].x - heroCenterXPosition), this.gridOrdinate(element.platform.vertices[3].y), this.unity, this.unity) ;
                        // By construction of square, the 4's point of a platform is the upper left
                        platformDraw.closePath() ;
                    } else if (element instanceof peak) {
                        peakDraw.moveTo(this.gridAbscissa(element.peak.vertices[0].x - heroCenterXPosition), this.gridOrdinate(element.peak.vertices[0].y)) ;
                        peakDraw.lineTo(this.gridAbscissa(element.peak.vertices[1].x - heroCenterXPosition), this.gridOrdinate(element.peak.vertices[1].y)) ;
                        peakDraw.lineTo(this.gridAbscissa(element.peak.vertices[2].x - heroCenterXPosition), this.gridOrdinate(element.peak.vertices[2].y)) ;
                        peakDraw.closePath() ;
                    }
                })
            }
        })

        this.ctx.fillStyle = "red" ;
        this.ctx.strokeStyle = "green" ;
        this.ctx.fill(platformDraw) ;
        this.ctx.fill(peakDraw) ;
        this.ctx.stroke(peakDraw) ;
        this.ctx.stroke(platformDraw) ;
    }
}


// main

const canvasGame = document.getElementById("canvas");
const ctx = canvasGame.getContext("2d");

const drawingInstance = new drawing();
const heroInstance = new hero();

const keys = {};

const gridInstance = new grid() ;
gridInstance.defaultGrid(1200) ;
let peak1 ;

for (let k = 4; k < 300; k = k+4) {
    peak1 = new peak(k,8,"down") ;
    gridInstance.addPeak(peak1) ;
}

drawingInstance.drawGrid() ;
drawingInstance.plotHero(heroInstance) ;
drawingInstance.plotGrid(gridInstance, heroInstance) ;

let cpt = 0
let stepTimeout = 0 ;
function moveAfterJump() {
    if (cpt < 10*drawingInstance.fps && !heroInstance.isDead) {
        if(keys.Space) {
            jump() ;
        }
        heroInstance.move(drawingInstance, gridInstance) ;
        ctx.clearRect(0,0, ctx.width, ctx.height)
        drawingInstance.drawGrid() ;
        drawingInstance.plotHero(heroInstance) ;
        drawingInstance.plotGrid(gridInstance, heroInstance) ;
        cpt ++;
        stepTimeout = setTimeout(moveAfterJump, 1/drawingInstance.fps * 1000);
    } else {
        titi = Date.now() ;
        console.log(titi - toto)
    }

}

function move() {
    moveAfterJump();
}

function print() {
    console.log(heroInstance)
}

function jump() {
    heroInstance.jump() ;
}

const angle = { angle: Math.PI/6} ;

function rotate() {
    heroInstance.rotate(angle.angle)
    ctx.clearRect(0,0, ctx.width, ctx.height)
    drawingInstance.drawGrid() ;
    drawingInstance.plotHero(heroInstance) ;
    drawingInstance.plotGrid(gridInstance, heroInstance) ;
}

function keyEventHandler(event){
    if (!heroInstance.hasStarted) {
        heroInstance.hasStarted = true ;
        move()
    }
    keys[event.code] = event.type === "keydown";
    event.preventDefault();
}
window.addEventListener("keydown",keyEventHandler);
window.addEventListener("keyup",keyEventHandler);



/*toto = new hero() ;

console.log(toto.body) ;
console.log(toto.foot[0]) ;
console.log(toto.foot[1]) ;

toto.rotate(Math.PI/4) ;
console.log(toto.body) ;
console.log(toto.foot[0]) ;
console.log(toto.foot[1]) ;

titi = new platform(3,3) ;
console.log(titi.platform.vertices)
console.log(titi.roof.vertices)
*/

/*toto = new grid() ;
toto.defaultGrid(10) ;
console.log(toto.grid[3][0])*/







