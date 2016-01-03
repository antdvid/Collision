#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/intersections.h>
#include <vector>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <FronTier.h>
#include "collid.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel   Kernel;
typedef Kernel::Point_3                                       Point_3;
typedef Kernel::Triangle_3                                    Triangle_3;

//define default parameters for collision detection
const double EPS = 0.000001;
const double DT  = 0.001;
double CollisionSolver::s_eps = EPS;
double CollisionSolver::s_thickness = 0.001;
double CollisionSolver::s_dt = DT;
double CollisionSolver::s_k = 1000;
double CollisionSolver::s_m = 0.01;
double CollisionSolver::s_lambda = 0.02;
bool   CollisionSolver3d::s_detImpZone = false;

double traitsForProximity::s_eps = EPS;
double traitsForCollision::s_eps = EPS;
double traitsForCollision::s_dt = DT;

//functions in the abstract base class
//set rounding tolerance
void CollisionSolver::setRoundingTolerance(double neweps){
	s_eps = neweps;
	traitsForProximity::s_eps = neweps;	
	traitsForCollision::s_eps = neweps;
}
double CollisionSolver::getRoundingTolerance(){return s_eps;}

//set fabric thickness
void CollisionSolver::setFabricThickness(double h){s_thickness = h;}
double CollisionSolver::getFabricThickness(){return s_thickness;}

//this function should be called at every time step
void CollisionSolver::setTimeStepSize(double new_dt){	
	s_dt = new_dt;
	traitsForCollision::s_dt = new_dt;
}
double CollisionSolver::getTimeStepSize(){return s_dt;}

//set spring constant
void   CollisionSolver::setSpringConstant(double new_k){s_k = new_k;}
double CollisionSolver::getSpringConstant(){return s_k;}

//set spring friction 
void   CollisionSolver::setFrictionConstant(double new_la){s_lambda = new_la;}
double CollisionSolver::getFrictionConstant(){return s_lambda;}

//set mass of fabric point
void   CollisionSolver::setPointMass(double new_m){s_m = new_m;}
double CollisionSolver::getPointMass(){return s_m;}

//functions in CollisionSolver3d
//assemble tris list from input intfc
void CollisionSolver3d::assembleFromInterface(
	const INTERFACE* intfc,const double dt)
{
	SURFACE** s;
	TRI *tri;
	setTimeStepSize(dt);
	trisList.clear();
	triPairList.clear();
	intfc_surface_loop(intfc,s)
	{
	    if (is_bdry(*s)) continue;
	    unsort_surf_point(*s);
	    surf_tri_loop(*s,tri)
	    {
	        trisList.push_back(tri);
	    }
	}
}

/*This function should be called before 
 *computing cloth internal dynamics*/
void CollisionSolver3d::recordOriginVelocity(){
	POINT* pt;
	STATE* sl, *sr;
	for (std::vector<TRI*>::iterator it = trisList.begin();
                it != trisList.end(); ++it){
	    for (int i = 0; i < 3; ++i){
		pt = Point_of_tri(*it)[i];
		sl = (STATE*)left_state(pt); 
		sr = (STATE*)right_state(pt);
		for (int j = 0; j < 3; ++j)
		    sl->x_old[j] = Coords(pt)[j];
	    }
	}
}

void CollisionSolver3d::computeAverageVelocity(){
	POINT* pt;
        STATE* sl; 
	double dt = getTimeStepSize();
        for (std::vector<TRI*>::iterator it = trisList.begin();
                it != trisList.end(); ++it){
            for (int i = 0; i < 3; ++i){
                pt = Point_of_tri(*it)[i];
                sl = (STATE*)left_state(pt); 
                for (int j = 0; j < 3; ++j)
		{
		    if (dt != 0)
                        sl->avgVel[j] = (Coords(pt)[j] - sl->x_old[j])/dt;
		    else
		        sl->avgVel[j] = sl->vel[j];
		}
            }
        }
	//restore coords of points to old coords !!!
	//x_old is the only valid coords for each point 
	//Coords(point) is for temporary judgement
	for (std::vector<TRI*>::iterator it = trisList.begin();
                it != trisList.end(); ++it){
            for (int i = 0; i < 3; ++i){
                pt = Point_of_tri(*it)[i];
                sl = (STATE*)left_state(pt);
                for (int j = 0; j < 3; ++j)
                    Coords(pt)[j] =  sl->x_old[j];
            }
        }
}

static void writeTriCoords(TRI* tri)
{
	for (int i = 0; i < 3; ++i)
	{
		std::cout << "[ ";
		for (int j = 0; j < 3; ++j)
		{
		    std::cout << Coords(Point_of_tri(tri)[i])[j] << " ";
		}
		std::cout << "] ";
	}
}

void CollisionSolver3d::printProximity()
{
	for (std::vector<TRI_PAIR>::iterator it = triPairList.begin(); 
		it != triPairList.end(); ++it)
	{
		std::cout<<"#" << it-triPairList.begin();
		writeTriCoords(it->first);
		std::cout << "intersects with " << std::endl;
		writeTriCoords(it->second);
		std::cout << std::endl;
	}
}

void CollisionSolver3d::printCollision()
{
	for (std::vector<TRI_PAIR>::iterator it = triPairList.begin(); 
		it != triPairList.end(); ++it)
	{
		std::cout<<"#" << it-triPairList.begin() << std::endl;
		writeTriCoords(it->first);
		std::cout << " collides with ";
		writeTriCoords(it->second);
	}
}

void CollisionSolver3d::getTriPairList(std::vector<std::pair<TRI*,TRI*> > &pairList)
{
	pairList =  triPairList;
}

void CollisionSolver3d::detectProximity()
{
	triPairList.clear();
	CGAL::box_self_intersection_d(trisList.begin(),trisList.end(),
                                     reportProximity<TRI*>(triPairList),traitsForProximity());
	updateAverageVelocity();
}

void CollisionSolver3d::detectCollision()
{
	bool is_collision = true; 
	const int MAX_ITER = 2;
	int niter = 0;

	std::cout<<"Starting collision handling: "<<std::endl;
	while(is_collision){
	    is_collision = false;
	    triPairList.clear();
	    CGAL::box_self_intersection_d(trisList.begin(),
		  trisList.end(),reportCollision<TRI*>(triPairList,is_collision),
		  traitsForCollision());
	    updateAverageVelocity();
	    std::cout<<"    #"<<niter << ": " << triPairList.size() 
		     << " pair of collision tris" << std::endl;
	    if (++niter > MAX_ITER) break;
	}
	if (is_collision) 
	    computeImpactZone();
}

void CollisionSolver3d::turnOffImpZone(){s_detImpZone = false;}
void CollisionSolver3d::turnOnImpZone(){s_detImpZone = true;}
bool CollisionSolver3d::getImpZoneStatus(){ return s_detImpZone;}
//this function is needed if collision still happens
//after several iterations;
void CollisionSolver3d::computeImpactZone()
{
	bool is_collision = true;
        int numZones = 0;
	int niter = 0;

        std::cout<<"Starting compute Impact Zone: "<<std::endl;
	turnOnImpZone();
	makeSet(trisList);             
        while(is_collision){
            is_collision = false;
            triPairList.clear();

	    //start UF alogrithm
	    //merge four pts if collision happens

	    CGAL::box_self_intersection_d(trisList.begin(),
                  trisList.end(),reportCollision<TRI*>(triPairList,is_collision),
                  traitsForCollision());

            updateAverageVelocity();
	    updateImpactZoneVelocity(numZones);
            std::cout <<"    #"<<niter++ << ": " << triPairList.size()
                      << " pair of collision tris" << std::endl;
	    std::cout <<"    # "<< numZones
		      <<" zones of impact" << std::endl;
        }
	turnOffImpZone();
	return;
}

void CollisionSolver3d::updateImpactZoneVelocity(int &numZones)
{
	POINT* pt;
	numZones = 0;
	unsortTriList(trisList);
	for (std::vector<TRI*>::iterator it = trisList.begin();
	     it != trisList.end(); ++it){
	    for (int i = 0; i < 3; ++i){
		pt = Point_of_tri(*it)[i];
		//skip traversed or isolated pts
		if (sorted(pt) ||
		    weight(findSet(pt)) == 1) continue;
		else{
		    updateImpactListVelocity(findSet(pt));
		    numZones++;
		}
	    }
	}	
}

inline double myDet3d(double a[3][3]){
    return  a[0][0]*(a[1][1]*a[2][2] - a[2][1]*a[1][2]) 
	  - a[0][1]*(a[1][0]*a[2][2] - a[2][0]*a[1][2]) 
	  + a[0][2]*(a[1][0]*a[2][1] - a[2][0]*a[1][1]);
}

void CollisionSolver3d::updateImpactListVelocity(POINT* head){
	STATE* sl = NULL;
	POINT* p = head;
	double m = getPointMass();
	double total_m = 0.0;
	double x_cm[3] = {0.0}, v_cm[3] = {0.0};
	double L[3] = {0.0}; //angular momentum
	double I[3][3] = {0.0}; //inertia tensor
	double tmp[3][3];
	int num_pts = 0;

	while(p){
		num_pts++; //debug
		sorted(p) = YES;
		sl = (STATE*)left_state(p);
		for (int i = 0; i < 3; ++i){
		    x_cm[i] += m*sl->x_old[i]; 
		    v_cm[i] += m*sl->avgVel[i];
		}
		total_m += m;
                p = next_pt(p);
        }
	//compute center and veclocity of impact Zone
	for (int i = 0; i < 3; ++i){
	    x_cm[i] /= total_m;
	    v_cm[i] /= total_m;
	}

	//compute angular momentum
	p = head;
	while(p){
	    double dx[3], dv[3], Li[3];
	    sl = (STATE*)left_state(p);
	    minusVec(sl->x_old,x_cm,dx);
	    minusVec(sl->avgVel,v_cm,dv); 	
	    Cross3d(dx,dv,Li);
	    scalarMult(m,Li,Li);
	    addVec(Li,L,L);    
	    p = next_pt(p);
	}
	//compute Inertia tensor
	p = head;
	while(p){
	    double dx[3], mag_dx = 0.0;
	    sl = (STATE*)left_state(p);
	    minusVec(sl->x_old,x_cm,dx);
	    mag_dx = Mag3d(dx);
	    for (int i = 0; i < 3; ++i)
	    for (int j = 0; j < 3; ++j){
		tmp[i][j] = -dx[i]*dx[j];
		if (i == j)
		    tmp[i][j] += mag_dx*mag_dx; 
	 	I[i][j] += tmp[i][j]*m;
	    } 
	    p = next_pt(p);
	}

	//compute angular velocity w: I*w = L;
	double w[3], mag_w = 0;
	for (int i = 0; i < 3; ++i){
	    memcpy(tmp,I,9*sizeof(double));
	    for (int j = 0; j < 3; j++)
		tmp[j][i] = L[j];
	    w[i] = myDet3d(tmp)/myDet3d(I);
	}
	mag_w = Mag3d(w);
	
	//compute average velocity for each point
	double dt = getTimeStepSize();
	p = head;
        while(p){
	    double x_new[3],dx[3];
	    double xF[3], xR[3];
	    sl = (STATE*)left_state(p);
	    minusVec(sl->x_old,x_cm,dx);
	    scalarMult(Dot3d(dx,w)/Dot3d(w,w),w,xF);
	    minusVec(dx,xF,xR);
	    for (int i = 0; i < 3; ++i)
	    {
		double wxR[3],tmpV[3];
		scalarMult(sin(dt*mag_w)/mag_w,w,tmpV);
		Cross3d(tmpV,xR,wxR);
	    	x_new[i] = x_cm[i] + dt*v_cm[i]
			 + xF[i] + cos(dt*mag_w)*xR[i]
			 + wxR[i];
		STATE* sl = (STATE*)left_state(p);
		sl->avgVel[i] = (x_new[i] - sl->x_old[i])/dt;
	    	if (isnan(sl->avgVel[i]))
		{ 
			printf("num_pts = %d, weight = %d\n",
			num_pts,weight(head));
			printf("nan vel, w = %f, mag_w = %f\n",
			w[i],mag_w);
			printf("L = [%f %f %f]\n",L[0],L[1],L[2]);
			printf("I = [%f %f %f;  %f %f %f; %f %f %f]\n",
			I[0][0],I[0][1],I[0][2],I[1][0],I[1][1],I[1][2],
			I[2][0],I[2][1],I[2][2]);
			printf("xF = %f %f %f, xR = %f %f %f\n",
			xF[0],xF[1],xF[2],xR[0],xR[1],xR[2]);
		}
	    }
	    p = next_pt(p);
	}
	//done!!!
}

//resolve collision in the input tris list
void CollisionSolver3d::resolveCollision()
{
	STATE* sl; 
	POINT* pt;
	pt = Point_of_tri(trisList[0])[0];
	sl = (STATE*)left_state(pt);
	//compute average velocity
	computeAverageVelocity();

	//test proximity for tris on surfaces
	detectProximity();

	//test collision for tri pairs
	detectCollision();

	//update position using average velocity
	updateFinalPosition();

	//update velocity using average velocity
	updateFinalVelocity();
}

extern bool isCollision(const TRI* a, const TRI* b){
	if (a->surf == b->surf && 
	    (wave_type(Hyper_surf(a->surf)) == NEUMANN_BOUNDARY ||
	     wave_type(Hyper_surf(a->surf)) == MOVABLE_BODY_BOUNDARY))
	     return false;
	for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            if (Point_of_tri(a)[i] == Point_of_tri(b)[j])
                return false;
        }
	if (MovingTriToTri(a,b,CollisionSolver3d::getRoundingTolerance())) 
	     return true;
	else
	     return false;			
}

extern bool isProximity(const TRI* a, const TRI* b){
	if (a->surf == b->surf &&
  	    (wave_type(Hyper_surf(a->surf)) == NEUMANN_BOUNDARY ||
  	    wave_type(Hyper_surf(a->surf)) == MOVABLE_BODY_BOUNDARY ||
	    wave_type(Hyper_surf(a->surf)) == FIRST_PHYSICS_WAVE_TYPE))
		return false;
	for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	{
	    if (Point_of_tri(a)[i] == Point_of_tri(b)[j])
		return false;
	}
	if (TriToTri(a,b,CollisionSolver3d::getFabricThickness()))
	    return true;
	else
	    return false;
}

//helper function to detect a collision between 
//a moving point and a moving triangle
//or between two moving edges 
static bool MovingTriToTri(const TRI* a,const TRI* b, const double h)
{
	bool is_detImpZone = CollisionSolver3d::getImpZoneStatus();
	POINT* pts[4];
	bool status = false;
	//detect point to tri collision
	for (int i = 0; i < 3; ++i){
	    for (int j = 0; j < 3; ++j)
	    	pts[j] = Point_of_tri(b)[j];
	    pts[3] = Point_of_tri(a)[i];
	    if(MovingPointToTri(pts,h)) status = true;
	    if (status && is_detImpZone)
		createImpZone(pts,4);

	    for (int j = 0; j < 3; ++j)
		pts[j] = Point_of_tri(a)[j];
	    pts[3] = Point_of_tri(b)[i];
	    if(MovingPointToTri(pts,h)) status = true;
	    if (status && is_detImpZone)
		createImpZone(pts,4);
	}

	//detect edge to edge collision
	for (int i = 0; i < 3; ++i)
        {
            pts[0] = Point_of_tri(a)[i];
            pts[1] = Point_of_tri(a)[(i+1)%3];
            for (int j = 0; j < 3; ++j){
                pts[2] = Point_of_tri(b)[j];
                pts[3] = Point_of_tri(b)[(j+1)%3];
                if(MovingEdgeToEdge(pts,h)) status = true;
	        if (status && is_detImpZone)
		    createImpZone(pts,4);
	    }
        }
	return status;
}

static void createImpZone(POINT* pts[], int num = 4){
	for (int i = 0; i < num-1; ++i)
	    mergePoint(pts[i],pts[i+1]); 
}

static bool MovingPointToTri(POINT* pts[],const double h){
	double dt = CollisionSolver3d::getTimeStepSize();
	double roots[4] = {-1,-1,-1,dt};
	STATE* sl;
	if (isCoplanar(pts,h,dt,roots)){
	    for (int i = 0; i < 4; ++i){
		if (roots[i] < 0) continue;
		for (int j = 0; j < 4; ++j){
		    sl = (STATE*)left_state(pts[j]);
		    for (int k = 0; k < 3; ++k)
		        Coords(pts[j])[k] = sl->x_old[k]+roots[i]*sl->avgVel[k];	
		}
		if (PointToTri(pts,h)) 
			return true;
	    }
	    return false;
	}
	else
	    return false;
}

static bool MovingEdgeToEdge(POINT* pts[],const double h){
	double dt = CollisionSolver3d::getTimeStepSize();
	double roots[4] = {-1,-1,-1,dt};
	STATE* sl;
	if (isCoplanar(pts,h,dt,roots)){
            for (int i = 0; i < 4; ++i){
                if (roots[i] < 0) continue;
                for (int j = 0; j < 4; ++j){
                    sl = (STATE*)left_state(pts[j]);
                    for (int k = 0; k < 3; ++k)
                        Coords(pts[j])[k] = sl->x_old[k]+roots[i]*sl->avgVel[k];
                }
                if (EdgeToEdge(pts,h))
                        return true;
            }
	    return false;
        }
	else
	    return false;
}

static bool isCoplanar(POINT* pts[], const double h, const double dt, double roots[])
{
	double v[4][3], x[4][3];
	for (int i = 0; i < 4; ++i){
	    STATE* sl = (STATE*)left_state(pts[i]);
	    for (int j = 0; j < 3; ++j){
	        v[i][j] = sl->avgVel[j];
		x[i][j] = sl->x_old[j];
	    }
	}
	for (int i = 1; i < 4; ++i){
	    for (int j = 0; j < 3; ++j){
		v[i][j] = v[0][j] - v[i][j];
		x[i][j] = x[0][j] - x[i][j];
	    }
	}
	//get roots "t" of a cubic equation
	//(x1+tv1)x(x2+tv2)*(x3+tv3) = 0
	//transform to at^3+bt^2+ct+d = 0
	double a, b, c, d;
	a = v[3][2]*(v[1][0]*v[2][1]-v[1][1]*v[2][0])-
            v[3][1]*(v[1][0]*v[2][2]-v[1][2]*v[2][0])+
            v[3][0]*(v[1][1]*v[2][2]-v[1][2]*v[2][1]);
        b = v[3][2]*(v[1][0]*x[2][1]-v[1][1]*x[2][0]-v[2][0]*x[1][1]+v[2][1]*x[1][0])-
            v[3][1]*(v[1][0]*x[2][2]-v[1][2]*x[2][0]-v[2][0]*x[1][2]+v[2][2]*x[1][0])+
            v[3][0]*(v[1][1]*x[2][2]-v[1][2]*x[2][1]-v[2][1]*x[1][2]+v[2][2]*x[1][1])+
            x[3][2]*(v[1][0]*v[2][1]-v[1][1]*v[2][0])-
            x[3][1]*(v[1][0]*v[2][2]-v[1][2]*v[2][0])+
            x[3][0]*(v[1][1]*v[2][2]-v[1][2]*v[2][1]);
        c = x[3][2]*(v[1][0]*x[2][1]-v[1][1]*x[2][0]-v[2][0]*x[1][1]+v[2][1]*x[1][0])-
            x[3][1]*(v[1][0]*x[2][2]-v[1][2]*x[2][0]-v[2][0]*x[1][2]+v[2][2]*x[1][0])+
            x[3][0]*(v[1][1]*x[2][2]-v[1][2]*x[2][1]-v[2][1]*x[1][2]+v[2][2]*x[1][1])+
            v[3][2]*(x[1][0]*x[2][1]-x[1][1]*x[2][0])-
            v[3][1]*(x[1][0]*x[2][2]-x[1][2]*x[2][0])+
            v[3][0]*(x[1][1]*x[2][2]-x[1][2]*x[2][1]);
        d = x[3][2]*(x[1][0]*x[2][1]-x[1][1]*x[2][0])-
            x[3][1]*(x[1][0]*x[2][2]-x[1][2]*x[2][0])+
            x[3][0]*(x[1][1]*x[2][2]-x[1][2]*x[2][1]);
	//solve equation using method from "Art of Scientific Computing"
	//transform equation to t^3+at^2+bt+c = 0
	if (a != 0){
	    b /= a; c /= a; d /= a;
	    a = b; b = c; c = d;
	    double Q, R, theta;
	    Q = (a*a-3*b)/9;
	    R = (2*a*a*a-9*a*b+27*c)/54;
	    if (R*R < Q*Q*Q){
		theta = acos(R/sqrt(Q*Q*Q));
		roots[0] = -2*sqrt(Q)*cos(theta/3)-a/3;
		roots[1] = -2*sqrt(Q)*cos((theta+2*M_PI)/3)-a/3;
		roots[2] = -2*sqrt(Q)*cos((theta-2*M_PI)/3)-a/3;	
	    }
	    else{
		double A, B;
		double sgn = (R > 0) ? 1.0 : -1.0;
		A = -sgn*pow(fabs(R)+sqrt(R*R-Q*Q*Q),1.0/3.0);
		B = (A == 0) ? 0.0 : Q/A;
		roots[0] = (A+B)-a/3.0;
		if (A == B)
		    roots[1] = roots[2] = -0.5*(A+B)-a/3.0; //multiple roots
		else
		    roots[1] = roots[2] = -1; //complex roots, discard
	    }
	}
	else{
		a = b; b = c; c = d;
	   	double delta = b*b-4.0*a*c;
	   	if (a != 0 && delta > 0){
		    roots[0] = (-b+sqrt(delta))/(2.0*a);
	    	    roots[1] = (-b-sqrt(delta))/(2.0*a);
		    roots[2] = -1;
	   	}
		else if (a == 0 && b != 0)
		{
		    roots[0] = -c/b;
		    roots[1] = roots[2] = -1;
	        }
	   	else
		    roots[0] = roots[1] = roots[2] = -1;//complex roots, discard
	}
	//select and sort roots;
	for (int i = 0; i < 3; ++i){
	    	if (roots[i] < 0 || roots[i] > dt) 
		    roots[i] = -1;
	}
	for (int i = 1; i < 3; ++i)
	for (int k = i; k > 0 && roots[k] < roots[k-1]; k--)
		std::swap(roots[k],roots[k-1]);

	for (int i = 1; i < 3; ++i)
	if (roots[i] != -1) return true;
	return false;
}

//helper function to detect proximity between elements 
static bool TriToTri(const TRI* tri1, const TRI* tri2, double h){
	POINT* pts[4];
	STATE* sl[2];
	bool status = false;
	//make sure the coords are old coords;
	for (int i = 0; i < 3; ++i){
	    pts[0] = Point_of_tri(tri1)[i];
	    sl[0] = (STATE*)left_state(pts[0]);
	    pts[1] = Point_of_tri(tri2)[i];
	    sl[1] = (STATE*)left_state(pts[1]);
	    for (int j = 0; j < 2; ++j)
	    for (int k = 0; k < 3; ++k)
	        Coords(pts[j])[k] = sl[j]->x_old[k];
	}

	for (int i = 0; i < 3; ++i)
	{
	    for (int j = 0; j < 3; ++j)
		pts[j] = Point_of_tri(tri2)[j];
	    pts[3] = Point_of_tri(tri1)[i];
	    if (PointToTri(pts,h))
		status = true;
	    for (int j = 0; j < 3; ++j)
                pts[j] = Point_of_tri(tri1)[j];
            pts[3] = Point_of_tri(tri2)[i];
	    if(PointToTri(pts,h))
                status = true;
	}
	for (int i = 0; i < 3; ++i)
	{
	    pts[0] = Point_of_tri(tri1)[i];
	    pts[1] = Point_of_tri(tri1)[(i+1)%3];
	    for (int j = 0; j < 3; ++j){
		pts[2] = Point_of_tri(tri2)[j];
		pts[3] = Point_of_tri(tri2)[(j+1)%3];
	  	if (EdgeToEdge(pts, h))
		    status = true;
	    }  
	}
	return status;
}

static bool EdgeToEdge(POINT** pts, double h)
{
/*	x1	x3
 *	/	 \
 *     /	  \
 * x2 /		   \ x4
 * solve equation
 * x21*x21*a - x21*x43*b = x21*x31
 * -x21*x43*a + x43*x43*b = -x43*x31
 */
	double x21[3], x43[3], x31[3];
	double a, b;
	double tmp[3];
	double v1[3],v2[3];
	double nor[3], nor_mag, dist;

	Pts2Vec(pts[0],pts[1],x21);	    
	Pts2Vec(pts[2],pts[3],x43);
	Pts2Vec(pts[0],pts[2],x31);
	Cross3d(x21,x43,tmp);
	if (Mag3d(tmp) < CollisionSolver3d::getRoundingTolerance())
	    return false;

	a = (Dot3d(x43,x43)*Dot3d(x21,x31)-Dot3d(x21,x43)*Dot3d(x43,x31))/
	    (Dot3d(x21,x21)*Dot3d(x43,x43)-Dot3d(x21,x43)*Dot3d(x21,x43)); 
	b = (Dot3d(x21,x43)*Dot3d(x21,x31)-Dot3d(x21,x21)*Dot3d(x43,x31))/
            (Dot3d(x21,x21)*Dot3d(x43,x43)-Dot3d(x21,x43)*Dot3d(x21,x43));
	a = std::max(std::min(a,1.0),0.0);	
	b = std::max(std::min(b,1.0),0.0);	
	scalarMult(a,x21,v1);
	scalarMult(b,x43,v2);
	addVec(x31,v2,v2);
	for (int i = 0; i < 3; ++i)
	    nor[i] = v1[i] - v2[i];
	nor_mag = Mag3d(nor);
	if (nor_mag == 0)
	{
	    //v1 == v2; 
	    //two edges intersect with each other
	    //normal direction is relative velocity
	    if (DEBUGGING)
	    printf("Edge intersects with edge, nor = %f\n",
		    nor_mag);
	    STATE* sl[4];
	    for (int j = 0; j < 4; ++j)
	    {
		sl[j] = (STATE*)left_state(pts[j]);
	    }
	    for (int j = 0; j < 3; ++j)
            {
                nor[j]  = (1.0-a) * sl[0]->avgVel[j] + a * sl[1]->avgVel[j];
                nor[j] -= (1.0-b) * sl[2]->avgVel[j] + b * sl[3]->avgVel[j];
            }
	    nor_mag = Mag3d(nor);
	    for (int j = 0; j < 3; ++j)
		nor[j] /= nor_mag;
	}
	else
	for (int i = 0; i < 3; ++i)
	     nor[i] /= nor_mag;
	dist = distBetweenCoords(v1,v2);
	if (dist > 0.1 * h)
	    return false;
	if (DEBUGGING){
	std::cout<<"edge to edge collision" << std::endl;
	printf("nor = %f %f %f, a = %f, b = %f, dist = %f\n",
		nor[0],nor[1],nor[2],a,b,dist);
	}
	EdgeToEdgeImpulse(pts, nor, a, b, dist);
	return true;
}

static bool PointToTri(POINT** pts, double h)
{
/*	x1
 *  	/\     x4 *
 *     /  \
 * x2 /____\ x3
 *
 * solve equation
 * x13*x13*w1 + x13*x23*w2 = x13*x43
 * x13*x23*w1 + x23*x23*w2 = x23*x43
 */
	double w[3];
	double x13[3], x23[3], x43[3];
	double v1[3], v2[3];
	double nor[3], nor_mag, dist;

	Pts2Vec(pts[0],pts[2],x13);
	Pts2Vec(pts[1],pts[2],x23);
	Pts2Vec(pts[3],pts[2],x43);
	Cross3d(x13, x23, nor);
	STATE* sl = (STATE*)left_state(pts[0]);
	nor_mag = Mag3d(nor);
	for (int i = 0; i < 3; ++i)
	    nor[i] /= nor_mag;
	if (fabs(Dot3d(x43, nor)) > 0.1 * h)
	    return false;

	w[0] = (Dot3d(x13,x43)*Dot3d(x23,x23)-Dot3d(x23,x43)*Dot3d(x13,x23))/
               (Dot3d(x13,x13)*Dot3d(x23,x23)-Dot3d(x13,x23)*Dot3d(x13,x23));
	w[1] = (Dot3d(x13,x13)*Dot3d(x23,x43)-Dot3d(x13,x23)*Dot3d(x13,x43))/
	       (Dot3d(x13,x13)*Dot3d(x23,x23)-Dot3d(x13,x23)*Dot3d(x13,x23));
	w[2] = 1 - w[0] - w[1];
	for (int i = 0; i < 3; ++i)
	{
	    v1[i] = Coords(pts[3])[i];
	    v2[i] = 0.0;
	    for (int j = 0; j < 3; ++j)
		v2[i] += w[j] * Coords(pts[j])[i];
	}
	dist = distBetweenCoords(v1,v2);
	
//	double clength = std::sqrt(tri_area(tri));
	for (int i = 0; i < 3; ++i)
	{
	    if (w[i] > 1+h || w[i] < -h) //h should be h/clength, clength is too small 
		return false;
	}
	if (DEBUGGING){
		std::cout<<"point to tri collision" << std::endl;
		printf("nor = %f %f %f, w = [%f %f %f], dist = %f\n",
		 	nor[0],nor[1],nor[2],w[0],w[1],w[2],dist);
	}
	PointToTriImpulse(pts, nor, w, dist);
	return true;
}

/* repulsion and friction functions, update velocity functions */
static void PointToTriImpulse(POINT** pts, double* nor, double* w, double dist)
{
	STATE *sl[4], *sr[4];
	for (int i = 0; i < 4; ++i)
	{
	    sl[i] = (STATE*)left_state(pts[i]);
	    sr[i] = (STATE*)right_state(pts[i]);
	}
	double v_rel[3] = {0.0}, vn = 0.0, vt = 0.0;
	double impulse = 0.0, m_impulse = 0.0;
	double k = 500, m = 0.01,lambda = 0.02, dt;
	k      = CollisionSolver::getSpringConstant();
	m      = CollisionSolver::getPointMass();
	dt     = CollisionSolver::getTimeStepSize();
	lambda = CollisionSolver::getFrictionConstant(); 

	/* it is supposed to use the average velocity*/
	for (int i = 0; i < 3; ++i)
	{
	    for (int j = 0; j < 3; ++j)
		v_rel[i] += w[j] * sl[j]->avgVel[i];
	    v_rel[i] -= sl[3]->avgVel[i];
	}
	vn = -Dot3d(v_rel, nor);
	vt = sqrt(Dot3d(v_rel, v_rel) - sqr(vn));
	impulse += vn * 0.5;
	if (vn * dt < 0.1 * dist)
	    impulse += - std::min(dt*k*dist/m, (0.1*dist/dt - vn));
	m_impulse = 2.0 * impulse / (1.0 + Dot3d(w, w));

	/* it is supposed to modify the average velocity*/
	for (int i = 0; i < 3; ++i)
	{
	    for(int j = 0; j < 3; ++j)
	    {
		sl[i]->collsnImpulse[j] += w[i] * m_impulse * nor[j];
		if (fabs(vt) > 1.0e-10)
		    sl[i]->friction[j] += std::max(-fabs(lambda * w[i] * 
			m_impulse/vt), -1.0) * (v_rel[j] - vn * nor[j]);
		sl[i]->collsn_num += 1;
		sr[i]->collsnImpulse[j] = sl[i]->collsnImpulse[j];
		sr[i]->friction[j] = sl[i]->friction[j];
		sr[i]->collsn_num = sl[i]->collsn_num;
	    }
	}
	for (int j = 0; j < 3; ++j)
	{
	    sl[3]->collsnImpulse[j] -= m_impulse * nor[j];
	    if (fabs(vt) > 1.0e-10)
	        sl[3]->friction[j] += std::max(-fabs(lambda * 
			m_impulse/vt), -1.0) * (v_rel[j] - vn * nor[j]);
	    sl[3]->collsn_num += 1;
	    sr[3]->collsnImpulse[j] = sl[3]->collsnImpulse[j];
	    sr[3]->friction[j] = sl[3]->friction[j];
	    sr[3]->collsn_num = sl[3]->collsn_num;
	}
}

static void EdgeToEdgeImpulse(POINT** pts, double* nor, double a, double b, double dist)
{
	STATE *sl[4], *sr[4];
	for (int i = 0; i < 4; ++i)
	{
	    sl[i] = (STATE*)left_state(pts[i]);
	    sr[i] = (STATE*)right_state(pts[i]);
	}
	double v_rel[3] = {0.0, 0.0, 0.0}, vn = 0.0, vt = 0.0;
	double impulse = 0.0, m_impulse = 0.0;
	double k = 500, m = 0.01, lambda = 0.02, dt;
	k      = CollisionSolver::getSpringConstant();
	m      = CollisionSolver::getPointMass();
	dt     = CollisionSolver::getTimeStepSize();
	lambda = CollisionSolver::getFrictionConstant(); 

	/* it is supposed to use the average velocity*/
	for (int j = 0; j < 3; ++j)
	{
	    v_rel[j] += (1.0-a) * sl[0]->avgVel[j] + a * sl[1]->avgVel[j];
	    v_rel[j] -= (1.0-b) * sl[2]->avgVel[j] + b * sl[3]->avgVel[j];
	}
	vn = -Dot3d(v_rel, nor);
	vt = sqrt(Dot3d(v_rel, v_rel) - sqr(vn));
	impulse += vn * 0.5;
	if (vn * dt < 0.1 * dist)
	    impulse += - std::min(dt*k*dist/m, (0.1*dist/dt - vn));
	m_impulse = 2.0 * impulse / (a*a + (1.0-a)*(1.0-a) + b*b + (1.0-b)*(1.0-b));

	/* it is supposed to modify the average velocity*/
	for (int j = 0; j < 3; ++j)
	{
	    //p[0]
	    sl[0]->collsnImpulse[j] += (1.0 - a) * m_impulse * nor[j];
	    if (fabs(vt) > 1.0e-10)
	        sl[0]->friction[j] += std::max(-fabs(lambda * (1.0-a) *
			m_impulse/vt), -1.0) * (v_rel[j] - vn * nor[j]);
	    sl[0]->collsn_num += 1;

	    //p[1]
	    sl[1]->collsnImpulse[j] += a * m_impulse * nor[j];
	    if (fabs(vt) > 1.0e-10)
	        sl[1]->friction[j] += std::max(-fabs(lambda * a * 
			m_impulse/vt), -1.0) * (v_rel[j] - vn * nor[j]);
	    sl[1]->collsn_num += 1;

	    //p[2]
	    sl[2]->collsnImpulse[j] -= (1.0 - b) * m_impulse * nor[j];
	    if (fabs(vt) > 1.0e-10)
	        sl[2]->friction[j] += std::max(-fabs(lambda * (1.0 - b) *
			m_impulse/vt), -1.0)*(v_rel[j] - vn * nor[j]);
	    sl[2]->collsn_num += 1;

	    //p[3]
	    sl[3]->collsnImpulse[j] -= b * m_impulse * nor[j];
	    if (fabs(vt) > 1.0e-10)
	        sl[3]->friction[j] += std::max(-fabs(lambda * b * 
			m_impulse/vt), -1.0) * (v_rel[j] - vn * nor[j]);
	    sl[3]->collsn_num += 1;
	}
	for (int i = 0; i < 4; ++i)
	{
	    for (int j = 0; j < 3; ++j)
	    {
		sr[i]->collsnImpulse[j] = sl[i]->collsnImpulse[j];
		sr[i]->friction[j] = sl[i]->friction[j];
	    }
	    sr[i]->collsn_num = sl[i]->collsn_num;
	}
}

void CollisionSolver3d::updateFinalPosition()
{
	POINT* pt;
	STATE* sl;
	double dt = getTimeStepSize();
	for (std::vector<TRI*>::iterator it = trisList.begin();
	     it != trisList.end(); ++it)
	{
	    if (isRigidBody(*it)) 
		continue;
	    for (int i = 0; i < 3; ++i){
		pt = Point_of_tri(*it)[i];
		sl = (STATE*)left_state(pt);
	    	for (int j = 0; j < 3; ++j)
		{
		    Coords(pt)[j] = sl->x_old[j]+sl->avgVel[j]*dt;
		    if (isnan(Coords(pt)[j]))
			printf("nan coords, x_old = %f, avgVel = %f\n",
				sl->x_old[j],sl->avgVel[j]);
		}
	    }
	}
}

void CollisionSolver3d::updateFinalVelocity()
{
	//TODO:avgVel is actually the velocity at t(n+1/2)
	//need to call spring solver to get velocity at t(n+1)
	//for simplicity now set v(n+1) = v(n+1/2)
	POINT* pt;
	STATE* sl;
	for (std::vector<TRI*>::iterator it = trisList.begin();
             it != trisList.end(); ++it)
        {
            for (int i = 0; i < 3; ++i){
                pt = Point_of_tri(*it)[i];
                sl = (STATE*)left_state(pt);
                for (int j = 0; j < 3; ++j)
                    sl->vel[j] = sl->avgVel[j];
            }
        }
}

void CollisionSolver3d::updateAverageVelocity()
{
	POINT *p;
	STATE *sl, *sr;
	int n = triPairList.size();
	double maxSpeed = 0;
	double* maxVel = NULL;

	/*debugging*/
	if (DEBUGGING)
	if (n == 1){
	    printf("tri %p and %p intersect\n",
	    (void*)triPairList[0].first,(void*)triPairList[0].second);
	    gviewplotTriPair("triPair",triPairList[0]);
	}

	unsortTriList(trisList);
	for (int i = 0; i < n; ++i)
	{
	    for (int j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(triPairList[i].first)[j];
		if (sorted(p)) continue;
		sl = (STATE*)left_state(p);
		sr = (STATE*)right_state(p);
		if (sl->collsn_num > 0)
		{
		    for (int k = 0; k < 3; ++k)
		    {
			sl->avgVel[k] += (sl->collsnImpulse[k] + sl->friction[k])
					/sl->collsn_num;
			sr->avgVel[k] = sl->avgVel[k];
			if (isinf(sl->avgVel[k]) || isnan(sl->avgVel[k])) 
			    printf("inf/nan vel[%d]: impulse = %f, friction = %f, collsn_num = %d\n",
				k,sl->collsnImpulse[k],sl->friction[k],sl->collsn_num);
			//reset impulse and fricition to 0
			//collision handling will recalculate the impulse
			sl->collsnImpulse[k] = sl->friction[k] = 0.0;
		    }
		    sl->collsn_num = 0;
		}
		//debugging: print largest speed
		double speed = Mag3d(sl->avgVel);
		if (speed > maxSpeed)
			maxVel = sl->avgVel;
		//
		sorted(p) = YES;
	    }
	    for (int j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(triPairList[i].second)[j];
		if (sorted(p)) continue;
		sl = (STATE*)left_state(p);
		sr = (STATE*)right_state(p);
		if (sl->collsn_num > 0)
		{
		    for (int k = 0; k < 3; ++k)
		    {
			sl->avgVel[k] += (sl->collsnImpulse[k] + sl->friction[k])
					 /sl->collsn_num;
			if (isinf(sl->avgVel[k]) || isnan(sl->avgVel[k])) 
			    printf("inf/nan vel[%d]: impulse = %f, friction = %f, collsn_num = %d\n",
				k,sl->collsnImpulse[k],sl->friction[k],sl->collsn_num);
			sr->avgVel[k] = sl->avgVel[k];
			//reset impulse and friction to 0
			//collision handling will recalculate the impulse
			sl->collsnImpulse[k] = sl->friction[k] = 0.0;
		    }
		    sl->collsn_num = 0;
		}
		//debugging: print largest speed
		double speed = Mag3d(sl->avgVel);
		if (speed > maxSpeed)
			maxVel = sl->avgVel;
		//
		sorted(p) = YES;
	    }
	}
	if (DEBUGGING)
	if (maxVel != NULL)
	    printf("    max velocity = [%f %f %f]\n",maxVel[0],maxVel[1],maxVel[2]);
}

//gview functions for debugging
void CollisionSolver3d::gviewplotPair(const char *dname)
{
	size_t num_tris = triPairList.size();
	if (num_tris == 0)
		return;
		
	static const char *indent = "    ";
	int        i,j;
	char       fname[256];
	FILE       *file;
	POINT      *p;

	if (create_directory(dname,YES) == FUNCTION_FAILED)
	{
	    screen("WARNING in gview_plot_tri_list(), "
	           "directory %s doesn't exist and can't be created\n",dname);
	    return;
	}

	for (i = 0; i < num_tris; ++i)
	{
	(void) sprintf(fname,"%s/tri_%03d.list",dname,i);
	if ((file = fopen(fname,"w")) == NULL)
	{
	    screen("WARNING in gview_plot_tri_list(), "
	           "can't open %s\n",fname);
	    return;
	}
	(void) fprintf(file,"{ LIST\n");
	double *BBL = topological_grid(triPairList[0].first->surf->interface).GL;
        double *BBU = topological_grid(triPairList[0].first->surf->interface).GU;
        gview_bounding_box(file,BBL,BBU,1,indent);
	(void) fprintf(file,"%s{\n%s%sOFF\n%s%s%6d %6d %6d\n",
			indent,indent,indent,
			indent,indent,3,1,0);
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(triPairList[i].first)[j];
		(void) fprintf(file, "%s%s%-9g %-9g %-9g\n",
			indent,indent,
			Coords(p)[0],Coords(p)[1],Coords(p)[2]);
	    }
	    (void) fprintf(file,"%s%s%-4d %-4d %-4d %-4d ",indent,indent,
		3,0,1,2);
	    (void) fprintf(file,"0.0 0.0 1.0 1.0\n");
	(void) fprintf(file,"%s}\n",indent);
	(void) fprintf(file,"%s{\n%s%sOFF\n%s%s%6d %6d %6d\n",
			indent,indent,indent,
			indent,indent,3,1,0);
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(triPairList[i].second)[j];
		(void) fprintf(file, "%s%s%-9g %-9g %-9g\n",
			indent,indent,
			Coords(p)[0],Coords(p)[1],Coords(p)[2]);
	    }
	    (void) fprintf(file,"%s%s%-4d %-4d %-4d %-4d ",indent,indent,
		3,0,1,2);
	    (void) fprintf(file,"0.0 1.0 0.0 1.0\n");
	(void) fprintf(file,"%s}\n",indent);

	(void) fprintf(file,"}\n");
	(void) fclose(file);
	}
}

void CollisionSolver3d::gviewplotPairList(const char *dname)
{
	size_t num_tris = triPairList.size();
	if (num_tris == 0)
		return;
	TRI** tris = new TRI*[num_tris*2];
	for (size_t i = 0; i < num_tris; ++i)
	{
	    tris[i*2] = triPairList[i].first;
	    tris[i*2+1] = triPairList[i].second;
	}
	
	static const char *indent = "    ";
	int        i,j;
	char       fname[256];
	FILE       *file;
	POINT      *p;

	if (create_directory(dname,YES) == FUNCTION_FAILED)
	{
	    screen("WARNING in gview_plot_tri_list(), "
	           "directory %s doesn't exist and can't be created\n",dname);
	    return;
	}
	(void) sprintf(fname,"%s/collid.list",dname);
	if ((file = fopen(fname,"w")) == NULL)
	{
	    screen("WARNING in gview_plot_tri_list(), "
	           "can't open %s\n",fname);
	    return;
	}
	(void) fprintf(file,"{ LIST\n");
	double *BBL = topological_grid(tris[0]->surf->interface).GL;
        double *BBU = topological_grid(tris[0]->surf->interface).GU;
	gview_bounding_box(file,BBL,BBU,1,indent);
	(void) fprintf(file,"%s{\n%s%sOFF\n%s%s%lu %lu %d\n",
			indent,indent,indent,
			indent,indent,3*num_tris,num_tris,0);
	for (i = 0; i < num_tris; ++i)
	{
	    for (j = 0; j < 3; ++j)
	    {
		p = Point_of_tri(tris[i])[j];
		(void) fprintf(file, "%s%s%-9g %-9g %-9g\n",
			indent,indent,
			Coords(p)[0],Coords(p)[1],Coords(p)[2]);
	    }
	}
	for (i = 0; i < num_tris; ++i)
	{
	    (void) fprintf(file,"%s%s%-4d %-4d %-4d %-4d ",indent,indent,
		3,3*i,3*i+1,3*i+2);
	    (void) fprintf(file,"0.0 0.0 1.0 1.0\n");
	}
	(void) fprintf(file,"%s}\n",indent);
	(void) fprintf(file,"}\n");
	(void) fclose(file);

	delete[] tris;
}

/* Function from CGAL libary, not used anymore */
static bool trisIntersect(const TRI* a, const TRI* b)
{
	Point_3 pts[3];
	Triangle_3 tri1,tri2;
	for (int i = 0; i < 3; ++i)
	   pts[i] = Point_3(Coords(Point_of_tri(a)[i])[0],
                            Coords(Point_of_tri(a)[i])[1],
                            Coords(Point_of_tri(a)[i])[2]);
	tri1 = Triangle_3(pts[0],pts[1],pts[2]);
	for (int i = 0; i < 3; ++i)
	   pts[i] = Point_3(Coords(Point_of_tri(b)[i])[0],
			    Coords(Point_of_tri(b)[i])[1],
			    Coords(Point_of_tri(b)[i])[2]);
	tri2 = Triangle_3(pts[0],pts[1],pts[2]);
	return CGAL::do_intersect(tri1,tri2);
}

/* The followings are helper functions for vector operations. */
static void Pts2Vec(const POINT* p1, const POINT* p2, double* v){
	for (int i = 0; i < 3; ++i)	
	    v[i] = Coords(p2)[i] - Coords(p1)[i];
}

static double distBetweenCoords(double* v1, double* v2)
{
	double dist = 0.0;
	for (int i = 0; i < 3; ++i){
		dist += sqr(v1[i]-v2[i]);
	}
	return std::sqrt(dist);
}

static void addVec(double* v1, double* v2, double* ans)
{
	for (int i = 0; i < 3; ++i)
	    ans[i] = v1[i]+v2[i];
}

static void minusVec(double* v1, double* v2, double* ans)
{
	for (int i = 0; i < 3; ++i)
	    ans[i] = v1[i]-v2[i];
}

static void scalarMult(double a,double* v, double* ans)
{
	for (int i = 0; i < 3; ++i)
            ans[i] = a*v[i];	
}

static void unsort_surf_point(SURFACE *surf)
{
        TRI *tri;
        POINT *p;
        int i;

        for (tri = first_tri(surf); !at_end_of_tri_list(tri,surf);
                        tri = tri->next)
        {
            for (i = 0; i < 3; ++i)
            {
                p = Point_of_tri(tri)[i];
                sorted(p) = NO;
            }
        }
}       /* end unsort_surf_point */

static void unsortTriList(std::vector<TRI*>& trisList){
	POINT* pt;
	for (std::vector<TRI*>::iterator it = trisList.begin();
		it != trisList.end(); ++it){
	    for (int i = 0; i < 3; ++i){
		sorted(Point_of_tri(*it)[i]) = NO;
	    }
	}
}

static bool isRigidBody(TRI* tri){
	if (wave_type(Hyper_surf(tri->surf)) == NEUMANN_BOUNDARY ||
	    wave_type(Hyper_surf(tri->surf)) == MOVABLE_BODY_BOUNDARY)
	    return true;
	else
	    return false;
}

static void gviewplotTriPair(const char dname[], const TRI_PAIR& tri_pair){
	TRI **tris = new TRI*[2];
	tris[0] = tri_pair.first;
	tris[1] = tri_pair.second;
	gview_plot_tri_list(dname,tris,2);
	delete[] tris;
}

//functions for UF alogrithm
inline int& weight(POINT* p){
	STATE* sl = (STATE*)left_state(p);
	return sl->impZone.num_pts;
}

inline POINT*& root(POINT* p){
	STATE* sl = (STATE*)left_state(p);
	return sl->impZone.root;
}

inline POINT*& next_pt(POINT* p){
	STATE* sl = (STATE*)left_state(p);
        return sl->impZone.next_pt;
}

inline POINT*& tail(POINT* p){
	STATE* sl = (STATE*)left_state(p);
	return sl->impZone.tail;
}

static void makeSet(std::vector<TRI*>& trisList){
	STATE* sl;
	POINT* pt;
        for (std::vector<TRI*>::iterator it = trisList.begin();
                it != trisList.end(); ++it){
            for (int i = 0; i < 3; ++i){
		pt = Point_of_tri(*it)[i];
                sorted(pt) = NO;
		sl = (STATE*)left_state(pt);
		sl->impZone.next_pt = NULL;
		sl->impZone.tail = pt;
		sl->impZone.root = pt;
		sl->impZone.num_pts = 1;
            }
        }
}

static POINT* findSet(POINT* p){
	if (root(p) != p){
		root(p) = findSet(root(p));
	}
	else
	    return p;
}

static void mergePoint(POINT* X, POINT* Y){
	POINT* PX = findSet(X);
	POINT* PY = findSet(Y);
	STATE* s1, *s2;
	if (PX == PY) return;
	if (weight(PX) > weight(PY)){
	    //update root after merge
	    weight(PX) += weight(PY);
	    root(PY) = PX;
	    //link two list, update tail
	    next_pt(tail(PX)) = PY;
	    tail(PX) = tail(PY); 
	}
	else{
	    //update root after merge
	    weight(PY) += weight(PX);
	    root(PX) = PY;
	    //link two list, update tail
	    next_pt(tail(PY)) = PX;
	    tail(PY) = tail(PX); 
	}
}
//end of UF functions
