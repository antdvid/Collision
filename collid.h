//collid3d.h
#include <CGAL/Simple_cartesian.h>
#include <CGAL/box_intersection_d.h>
#include "../iFluid/ifluid_state.h"
#include <functional>

/*
user-defined state should include the following
struct UF{
	POINT* next_pt;
	POINT* root;
	POINT* tail;
	int num_pts;
};

struct STATE{
	double vel[3];
	double collsnImpulse[3];
	double friction[3];
	double avgVel[3];
	double x_old[3];
	int    collsn_num;
	bool   has_collsn;
	UF     impZone;
};*/

//abstract base class for hypersurface element(HSE)
//can be a point or a bond or a triangle
class CD_HSE{
public:
	virtual double max_static_coord(int) = 0;
	virtual double min_static_coord(int) = 0;
	virtual double max_moving_coord(int,double) = 0;
	virtual double min_moving_coord(int,double) = 0;
	virtual POINT* Point_of_hse(int) const  = 0;
	virtual int num_pts() const= 0;
	virtual ~CD_HSE(){};
};

//wrap class for triangle
class CD_TRI: public CD_HSE{
public:
	TRI* m_tri;
	CD_TRI(TRI* tri):m_tri(tri){}
	double max_static_coord(int);
	double min_static_coord(int);
	double max_moving_coord(int,double);
	double min_moving_coord(int,double);
	POINT* Point_of_hse(int) const;
	int num_pts() const {return 3;}
};

//wrap class for bond
class CD_BOND: public CD_HSE{
public:
	BOND* m_bond;
	CD_BOND(BOND* bond):m_bond(bond){}
	double max_static_coord(int);
	double min_static_coord(int);
	double max_moving_coord(int,double);
	double min_moving_coord(int,double);
	POINT* Point_of_hse(int) const;
	int num_pts()const{return 2;}
};

//wrap class for point
class CD_POINT: public CD_HSE{
public:
	POINT* m_point;
	CD_POINT(POINT* point):m_point(point){}
	double max_static_coord(int);
	double min_static_coord(int);
	double max_moving_coord(int,double);
	double min_moving_coord(int,double);
	POINT* Point_of_hse(int) const;
	int num_pts()const {return 1;}
};

//box traits structure for proximity detection 
template <int dim>
struct traitsForProximity{
	typedef double 		NT;
	typedef CD_HSE*		Box_parameter;
	typedef std::ptrdiff_t ID;

	static double s_eps;
	static int dimension(){ return dim;}
	static double min_coord(Box_parameter b, int d)
		{return b->min_static_coord(d)-s_eps;}
	static double max_coord(Box_parameter b, int d)
		{return b->max_static_coord(d)+s_eps;}
	static std::ptrdiff_t id(Box_parameter b) { return (std::ptrdiff_t)(b);}
};

//box traits structure for collision detection
template <int dim>
struct traitsForCollision{
	typedef double          		NT;
        typedef CD_HSE*      Box_parameter;

	static double s_eps;
	static double s_dt;
        static int dimension(){ return dim;}
        static double min_coord(Box_parameter b, int d)
		{return b->min_moving_coord(d,s_dt)-s_eps;}
        static double max_coord(Box_parameter b, int d)
		{return b->max_moving_coord(d,s_dt)+s_eps;}
	static std::ptrdiff_t id(Box_parameter b) { return (std::ptrdiff_t)(b);}
};

extern bool isProximity(const CD_HSE*,const CD_HSE*);
extern bool isCollision(const CD_HSE*,const CD_HSE*);

//callback functor to identify real collision
template <typename T>
struct reportProximity{
    int& num_pairs;
    reportProximity(int &npair): num_pairs(npair = 0){}
    // We write the elements with respect to 'boxes' to the output
    void operator()( const T &a, const T &b) {
	if(isProximity(a,b)){
	    num_pairs++;
	}
    }
};

template <typename T>
struct reportCollision{
    bool& is_collision;
    int&  num_pairs;
    reportCollision(bool &status, int &npairs) 
		   : is_collision(status), num_pairs(npairs = 0) {}
    // We write the elements with respect to 'boxes' to the output
    void operator()( const T &a, const T &b) {
	if (isCollision(a,b)){
	    num_pairs ++;
	    is_collision = true;
	}
    }
};

//abstract base class for collision detection and handling
class CollisionSolver {
private:
	//global parameters
	static double s_eps;
	static double s_thickness;
	static double s_dt;
	static double s_m;
	static double s_k;
	static double s_lambda;
protected:
	std::vector<CD_HSE*> hseList;
	void clearHseList();
public:
	static void setRoundingTolerance(double);
	static double getRoundingTolerance();
	static void setFabricThickness(double);
	static double getFabricThickness();
	static void setTimeStepSize(double);
	static double getTimeStepSize();
	static void setSpringConstant(double);
	static double getSpringConstant();
	static void setFrictionConstant(double);
	static double getFrictionConstant();
	static void setPointMass(double);
	static double getPointMass();
	virtual ~CollisionSolver() {}; //virtual destructor
	//pure virtual functions
	virtual void assembleFromInterface(const INTERFACE*,double dt) = 0;
	virtual void detectProximity() = 0;
	virtual void detectCollision() = 0;
	virtual void resolveCollision() = 0;
	virtual void recordOriginPosition() = 0;
};

//derived 2D-class for collision detection and handling
class CollisionSolver2d : public CollisionSolver {
public:
	void assembleFromInterface(const INTERFACE*,double dt);
	void detectProximity();
	void detectCollision();
	void resolveCollision();
	void getBondPairList(std::vector<std::pair<BOND*,BOND*> >&);
};

//derived 3D-class for collision detection and handling
class CollisionSolver3d : public CollisionSolver {
private:
	//input tris
	static bool s_detImpZone;
	void computeAverageVelocity();
	void updateAverageVelocity();
	void updateFinalVelocity();
	void updateFinalPosition();
	void computeImpactZone();
	void updateImpactZoneVelocity(int&);
	void updateImpactListVelocity(POINT*);
	void reduceSuperelast();
	bool reduceSuperelastOnce(int&);
public:
	//for debugging
	static int moving_edg_to_edg;
	static int moving_pt_to_tri;
	static int is_coplanar;
	static int edg_to_edg;
	static int pt_to_tri;
	static void printDebugVariable();
	void assembleFromInterface(const INTERFACE*,double dt);
	void recordOriginPosition();
	void detectProximity();
	void detectCollision();
	void resolveCollision();
	static void turnOffImpZone();
	static void turnOnImpZone();
	static bool getImpZoneStatus();
};

#if defined(isnan)
#undef isnan
#endif

#define DEBUGGING false
void initSurfaceState(SURFACE*,const double*);
void initTestModule(Front&, char*);
extern void testFunctions(POINT**,double);
inline void Pts2Vec(const POINT* p1, const POINT* p2, double* v); 
inline void scalarMult(double a,double* v, double* ans); 
inline void addVec(double* v1, double* v2, double* ans); 
inline void minusVec(double* v1, double* v2, double* ans); 
inline double distBetweenCoords(double* v1, double* v2);
