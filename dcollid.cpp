#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/intersections.h>
#include <vector>
#include <fstream>
#include <FronTier.h>
#include "collid.h"

static bool trisIntersect(const TRI* a, const TRI* b);
static bool TriToTri(const TRI* tri1, const TRI* tri2, double h);
static bool PointToTri(POINT** pts, double h);
static bool EdgeToEdge(POINT** pts, double h);
static void PointToTriUpdateVel(POINT** pts, double* nor, double* w, double dist);
static void EdgeToEdgeUpdateVel(POINT** pts, double* nor, double a, double b, double dist);

static bool TriPairToTriPair(const TRI_PAIR& a,const TRI_PAIR &b);

static void Pts2Vec(const POINT* p1, const POINT* p2, double* v);
static void scalarMult(double a,double* v, double* ans);
static void addVec(double* v1, double* v2, double* ans);
static double distBetweenCoords(double* v1, double* v2);
static void unsort_surf_point(SURFACE *surf);

typedef CGAL::Exact_predicates_inexact_constructions_kernel   Kernel;
typedef Kernel::Point_3                                       Point_3;
typedef Kernel::Triangle_3                                    Triangle_3;

//define rounding tolerance for collision detection
const double EPS = 0.000001;
double CollisionSolver::eps = EPS;
double CollisionSolver::thickness = 0.001;
double Traits_of_FT<TRI*>::eps = EPS;
double Traits_of_FT<TRI_PAIR>::eps = EPS;

//functions in the abstract base class
//set rounding tolerance
void CollisionSolver::setRoundingTolerance(double neweps)
{
	eps = neweps;
	Traits_of_FT<TRI*>::eps = eps;	
	Traits_of_FT<TRI_PAIR>::eps = eps;
}

double CollisionSolver::getRoundingTolerance()
{
	return eps;
}

//set fabric thickness
void CollisionSolver::setFabricThickness(double h)
{
	thickness = h;
}

double CollisionSolver::getFabricThickness()
{
	return thickness;
}

//functions in CollisionSolver3d
//assemble tris list from input intfc
void CollisionSolver3d::assembleFromInterface(
	const INTERFACE* intfc)
{
	SURFACE** s;
	TRI *tri;
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

void CollisionSolver3d::assembleFromTwoInterface(
	const INTERFACE* oldIntfc,const INTERFACE* newIntfc)
{
	SURFACE **oldS, **newS;
	TRI *oldTri, *newTri;

	triPairList.clear();
	triTwoPairList.clear();
	for ((oldS) = (oldIntfc)->surfaces, 
	     (newS) = (newIntfc)->surfaces; 
	     (oldS) && *(oldS) &&
	     (newS) && *(newS); 
	     ++(oldS), ++(newS)) 
	{
	    if (is_bdry(*oldS) || is_bdry(*newS)) continue;
	    for ((oldTri) = first_tri((*oldS)),
		 (newTri) = first_tri((*newS)); 
		!at_end_of_tri_list((oldTri),(*oldS)) &&
		!at_end_of_tri_list((newTri),(*newS));
		(oldTri) = (oldTri)->next,
		(newTri) = (newTri)->next)
	    {
	        triPairList.push_back(std::make_pair<TRI*,TRI*>(oldTri,newTri));
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
		    std::cout << Coords(tri->__pts[i])[j] << " ";
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
		std::cout << "intersects with" << std::endl;
		writeTriCoords(it->second);
		std::cout << std::endl;
	}
}

void CollisionSolver3d::printCollision()
{
	int count = 0;
	for (std::vector<std::pair<TRI_PAIR,TRI_PAIR> >::iterator it = triTwoPairList.begin(); 
		it != triTwoPairList.end(); ++it)
	{
		std::cout<<"#" << ++count << std::endl;
		std::cout<<(std::ptrdiff_t)(it->first.first);
		std::cout << " intersects with ";
		std::cout<<(std::ptrdiff_t)(it->second.first) << std::endl;
	}
}

void CollisionSolver3d::getProximityTrisPairList(std::vector<std::pair<TRI*,TRI*> > &pairList)
{
	pairList =  triPairList;	
}

void CollisionSolver3d::detectProximity()
{
	CGAL::box_self_intersection_d(trisList.begin(),trisList.end(),
                                     Report<TRI*>(triPairList),Traits_of_FT<TRI*>());
}

void CollisionSolver3d::detectCollision()
{
	CGAL::box_self_intersection_d(triPairList.begin(),triPairList.end(),
                                     Report<TRI_PAIR>(triTwoPairList),Traits_of_FT<TRI_PAIR>());
}

//resolve collision in the input tris list
void CollisionSolver3d::resolveCollision()
{
	//test proximity for tris on surfaces
	detectProximity();
	//test collision for tri pairs
	detectCollision();
	//print proximity tri pairs
	printProximity();
}

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

extern bool doIntersect(const TRI_PAIR& a, const TRI_PAIR& b){
	if (a.first->surf == b.first->surf && 
	    (wave_type(Hyper_surf(a.first->surf)) == NEUMANN_BOUNDARY ||
	     wave_type(Hyper_surf(a.first->surf)) == MOVABLE_BODY_BOUNDARY))
	     return false;
	else if (TriPairToTriPair(a,b)) 
	     return true;
	else
	     return false;			
}

extern bool doIntersect(const TRI* a, const TRI* b){
	if (a->surf == b->surf &&
  	    (wave_type(Hyper_surf(a->surf)) == NEUMANN_BOUNDARY ||
  	    wave_type(Hyper_surf(a->surf)) == MOVABLE_BODY_BOUNDARY ||
	    wave_type(Hyper_surf(a->surf)) == FIRST_PHYSICS_WAVE_TYPE))
		return false;
	for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	{
	    if (a->__pts[i] == b->__pts[j])
		return false;
	}
	if (TriToTri(a,b,CollisionSolver3d::getFabricThickness()))
	    return true;
	else
	    return false;
}

static bool TriPairToTriPair(const TRI_PAIR& a,const TRI_PAIR &b)
{
	return true;
}

//helper function to detect intersection or proximity between geometries
static bool TriToTri(const TRI* tri1, const TRI* tri2, double h){
	POINT* pts[4];
	for (int i = 0; i < 3; ++i)
	{
	    for (int j = 0; j < 3; ++j)
		pts[j] = tri2->__pts[j];
	    pts[3] = tri1->__pts[i];
	    if (PointToTri(pts,h))
		return true;
	    for (int j = 0; j < 3; ++j)
                pts[j] = tri1->__pts[j];
            pts[3] = tri2->__pts[i];
	    if(PointToTri(pts,h))
                return true;
	}
	for (int i = 0; i < 3; ++i)
	{
	    pts[0] = tri1->__pts[i];
	    pts[1] = tri1->__pts[(i+1)%3];
	    for (int j = 0; j < 3; ++j){
		pts[2] = tri2->__pts[j];
		pts[3] = tri2->__pts[(j+1)%3];
	  	if (EdgeToEdge(pts, h))
		    return true;
	    }  
	}
	return false;
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
	for (int i = 0; i < 3; ++i)
	    nor[i] /= nor_mag;
	dist = distBetweenCoords(v1,v2);
	if (dist > 0.1 * h)
	    return false;
	EdgeToEdgeUpdateVel(pts, nor, a, b, dist);
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
	PointToTriUpdateVel(pts, nor, w, dist);
	return true;
}

/* repulsion and friction functions, update velocity functions */
static void PointToTriUpdateVel(POINT** pts, double* nor, double* w, double dist)
{
	std::cout << "In PointToTirUpdateVel()" << std::endl;
	std::cout << "pts[3] before vel " << pts[3]->vel[0] << " " 
		  << pts[3]->vel[1] << " " 
		  << pts[3]->vel[2] << std::endl;
	std::cout << "sorted(pts[3]) = " << sorted(pts[3]) << std::endl;

	double v_rel[3] = {0.0}, vn = 0.0;
	double impulse = 0.0, m_impulse = 0.0;
	double k = 50000, m = 0.01, dt = 0.01;

	/* it is supposed to use the average velocity*/
	for (int i = 0; i < 3; ++i)
	{
	    for (int j = 0; j < 3; ++j)
		v_rel[i] += w[j] * pts[j]->vel[i];
	    v_rel[i] -= pts[3]->vel[i];
	}
	std::cout << "v_rel = " << v_rel[0] << " " 
		  << v_rel[1] << " " << v_rel[2] << std::endl;
	vn = Dot3d(v_rel, nor);
	if (vn < 0.0)
	    impulse += vn * 0.5;
	if (vn * dt < 0.1 * dist)
	    impulse += - std::min(dt*k*dist/m, (0.1*dist/dt - vn));
	m_impulse = 2.0 * impulse / (1.0 + Dot3d(w, w));

	/* it is supposed to modify the average velocity*/
	for (int i = 0; i < 3; ++i)
	{
	    if (sorted(pts[i])) continue;
	    for (int j = 0; j < 3; ++j)
		pts[i]->vel[j] += w[i] * m_impulse * nor[j];
	}
	if (!sorted(pts[3]))
	{
	    for (int j = 0; j < 3; ++j)
		pts[3]->vel[j] -= m_impulse * nor[j];
	}
	for (int i = 0; i < 4; ++i)
	    sorted(pts[i]) = YES;
	std::cout << "pts[3] after vel " << pts[3]->vel[0] << " " 
		  << pts[3]->vel[1] << " " 
		  << pts[3]->vel[2] << "\n" << std::endl;
}

static void EdgeToEdgeUpdateVel(POINT** pts, double* nor, double a, double b, double dist)
{
	std::cout << "In EdgeToEdgeUpdateVel()" << std::endl;
	std::cout << "pts[3] before vel " << pts[3]->vel[0] << " "
                  << pts[3]->vel[1] << " "
                  << pts[3]->vel[2] << std::endl;
        std::cout << "sorted(pts[3]) = " << sorted(pts[3]) << std::endl;

	double v_rel[3] = {0.0}, vn = 0.0;
	double impulse = 0.0, m_impulse = 0.0;
	double k = 50000, m = 0.01, dt = 0.01;

	/* it is supposed to use the average velocity*/
	for (int i = 0; i < 3; ++i)
	{
	    v_rel[i] += (1.0-a) * pts[0]->vel[i] + a * pts[1]->vel[i];
	    v_rel[i] -= (1.0-b) * pts[2]->vel[i] + b * pts[3]->vel[i];
	}
	vn = Dot3d(v_rel, nor);
	if (vn < 0.0)
	    impulse += vn * 0.5;
	if (vn * dt < 0.1 * dist)
	    impulse += - std::min(dt*k*dist/m, (0.1*dist/dt - vn));
	m_impulse = 2.0 * impulse / (a*a + (1.0-a)*(1.0-a) + b*b + (1.0-b)*(1.0-b));

	/* it is supposed to modify the average velocity*/
	for (int j = 0; j < 3; ++j)
	{
	    if (!sorted(pts[0]))
		pts[0]->vel[j] += (1.0 - a) * m_impulse * nor[j];
	    if (!sorted(pts[1]))
		pts[1]->vel[j] += a * m_impulse * nor[j];
	    if (!sorted(pts[2]))
		pts[2]->vel[j] -= (1.0 - b) * m_impulse * nor[j];
	    if (!sorted(pts[3]))
		pts[3]->vel[j] -= b * m_impulse * nor[j];
	}
	for (int i = 0; i < 4; ++i)
	    sorted(pts[i]) = YES;
	std::cout << "pts[3] after vel " << pts[3]->vel[0] << " "
                  << pts[3]->vel[1] << " "
                  << pts[3]->vel[2] << "\n" << std::endl;
}

/* Function from CGAL libary, not used anymore */
static bool trisIntersect(const TRI* a, const TRI* b)
{
	Point_3 pts[3];
	Triangle_3 tri1,tri2;
	for (int i = 0; i < 3; ++i)
	   pts[i] = Point_3(Coords(a->__pts[i])[0],
                            Coords(a->__pts[i])[1],
                            Coords(a->__pts[i])[2]);
	tri1 = Triangle_3(pts[0],pts[1],pts[2]);
	for (int i = 0; i < 3; ++i)
	   pts[i] = Point_3(Coords(b->__pts[i])[0],
			    Coords(b->__pts[i])[1],
			    Coords(b->__pts[i])[2]);
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
