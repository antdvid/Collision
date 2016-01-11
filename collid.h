//collid3d.h
#include <CGAL/Simple_cartesian.h>
#include <CGAL/box_intersection_d.h>
#include "../iFluid/ifluid_state.h"

typedef std::pair<TRI*,TRI*> TRI_PAIR;
typedef std::pair<BOND*,BOND*> BOND_PAIR;
/*
user-defined state should include the following
#ifndef FT_STATE
#define FT_STATE
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
	UF     impZone;
};
#endif*/
//box traits structure for proximity detection 
struct traitsForProximity{
	typedef double 		NT;
	typedef TRI*		Box_parameter;
	typedef std::ptrdiff_t ID;

	static double s_eps;
	static int dimension(){ return 3;}
	static double min_coord(Box_parameter b, int d){
		POINT** pts = Point_of_tri(b);
		STATE* ls[3] = {(STATE*)left_state(pts[0]),
				(STATE*)left_state(pts[1]),
				(STATE*)left_state(pts[2])};
		return std::min(std::min(
				ls[0]->x_old[d],
				ls[1]->x_old[d]),
				ls[2]->x_old[d])-s_eps;
	}
	static double max_coord(Box_parameter b, int d){
		POINT** pts = Point_of_tri(b);
		STATE* ls[3] = {(STATE*)left_state(pts[0]),
				(STATE*)left_state(pts[1]),
				(STATE*)left_state(pts[2])};
		return std::max(std::max(
				ls[0]->x_old[d],
				ls[1]->x_old[d]),
				ls[2]->x_old[d])+s_eps;
	}
	static std::ptrdiff_t id(Box_parameter b) { return (std::ptrdiff_t)(b);}
};

//box traits structure for collision detection
struct traitsForCollision{
	typedef double          		NT;
        typedef TRI*      Box_parameter;

	static double s_eps;
	static double s_dt;
        static int dimension(){ return 3;}
        static double min_coord(Box_parameter b, int d){
		double ans = INT_MAX;
		POINT** pts = Point_of_tri(b);
		for (int i = 0; i < 3; ++i){
			STATE* ls = (STATE*)left_state(pts[i]);
			ans = std::min(ans,ls->x_old[d]);
			ans = std::min(ans,ls->x_old[d]+s_dt*ls->avgVel[i]);
		}
		return ans-s_eps;
        }
        static double max_coord(Box_parameter b, int d){
		double ans = INT_MIN;
		POINT** pts = Point_of_tri(b);
		for (int i = 0; i < 3; ++i){
			STATE* ls = (STATE*)left_state(pts[i]);
			ans = std::max(ans,ls->x_old[d]);
			ans = std::max(ans,ls->x_old[d]+s_dt*ls->avgVel[i]);
		}
		return ans+s_eps;
        }
	static std::ptrdiff_t id(Box_parameter b) { return (std::ptrdiff_t)(b);}
};

extern bool isProximity(const TRI*,const TRI*);
extern bool isCollision(const TRI*,const TRI*);

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
	Front *front;
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
private:
	//input bonds
	std::vector<BOND*> bondsList;
	//bond pairs
	std::vector<BOND_PAIR> bondPairList;
public:
	void assembleFromInterface(const INTERFACE*,double dt);
	void detectProximity();
	void detectCollision();
	void resolveCollision();
	void gviewplotPairList(const char*);
	void gviewplotPair(const char *);
	void getBondPairList(std::vector<std::pair<BOND*,BOND*> >&);
};

//derived 3D-class for collision detection and handling
class CollisionSolver3d : public CollisionSolver {
private:
	//input tris
	std::vector<TRI*> trisList;
	static bool s_detImpZone;
	void computeAverageVelocity();
	void updateAverageVelocity();
	void updateFinalVelocity();
	void updateFinalPosition();
	void computeImpactZone();
	void updateImpactZoneVelocity(int&);
	void updateImpactListVelocity(POINT*);
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

static bool trisIntersect(const TRI* a, const TRI* b);
static bool TriToTri(const TRI* tri1, const TRI* tri2, double h);
static bool PointToTri(POINT** pts, double h);
static bool EdgeToEdge(POINT** pts, double h);
static void PointToTriImpulse(POINT** pts, double* nor, double* w, double dist);
static void EdgeToEdgeImpulse(POINT** pts, double* nor, double a, double b, double dist);

static bool isCoplanar(POINT*[], const double, double[]);
static bool MovingTriToTri(const TRI*,const TRI*,double);
static bool MovingPointToTri(POINT*[],const double);
static bool MovingEdgeToEdge(POINT*[],const double);

inline void Pts2Vec(const POINT* p1, const POINT* p2, double* v);
inline void scalarMult(double a,double* v, double* ans);
inline void addVec(double* v1, double* v2, double* ans);
inline void minusVec(double* v1, double* v2, double* ans);
inline double distBetweenCoords(double* v1, double* v2);
static void unsort_surface_point(SURFACE *surf);
static void unsortTriList(std::vector<TRI*>&);
static bool isRigidBody(const TRI*);
static bool isRigidBody(const POINT*);
static void gviewplotTriPair(const char[], const TRI_PAIR&);

static void makeSet(std::vector<TRI*>&);
static POINT* findSet(POINT*);
static void mergePoint(POINT*,POINT*);
static void createImpZone(POINT*[],int);
inline int& weight(POINT*);
inline POINT*& root(POINT*);
inline POINT*& next_pt(POINT*);
inline POINT*& tail(POINT*);

void initSurfaceState(SURFACE*,const double*);
void initTestModule(Front&, char*);
static void initBalls(Front&);
static void initPlane(Front&);
extern void testFunctions(POINT**,double);
static void printPointList(POINT**, const int);
