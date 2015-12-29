//collid3d.h
#include <CGAL/Simple_cartesian.h>
#include <CGAL/box_intersection_d.h>

typedef std::pair<TRI*,TRI*> TRI_PAIR;
typedef std::pair<BOND*,BOND*> BOND_PAIR;

struct STATE{
	double vel[3];
	double impulse[3];
	double friction[3];
	double avgVel[3];
	double x_old[3];
	int collsn_num;
};

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
    std::vector<std::pair<T,T> > *output;
    reportProximity(std::vector<std::pair<T,T> > &pairList) : output(&pairList) {}
    // We write the elements with respect to 'boxes' to the output
    void operator()( const T &a, const T &b) {
	if (isProximity(a,b))
    	    output->push_back(std::make_pair<T,T>(a,b));
    }
};

template <typename T>
struct reportCollision{
    std::vector<std::pair<T,T> > *output;
    reportCollision(std::vector<std::pair<T,T> > &pairList) : output(&pairList) {}
    // We write the elements with respect to 'boxes' to the output
    void operator()( const T &a, const T &b) {
	if (isCollision(a,b))
    	    output->push_back(std::make_pair<T,T>(a,b));
    }
};

//abstract base class for collision detection and handling
class CollisionSolver {
private:
	//rounding tolerance
	static double s_eps;
	static double s_thickness;
	static double s_dt;
protected:
	Front *front;
public:
	static void setRoundingTolerance(double);
	static double getRoundingTolerance();
	static void setFabricThickness(double);
	static double getFabricThickness();
	static void setTimeStepSize(double);
	static double getTimeStepSize();
	virtual ~CollisionSolver() {}; //virtual destructor
	//pure virtual functions
	virtual void assembleFromInterface(const INTERFACE*) = 0;
	virtual void detectProximity() = 0;
	virtual void detectCollision() = 0;
	virtual void resolveCollision() = 0;
	virtual void printProximity() = 0;
	virtual void printCollision() = 0;
	virtual void gviewplotPairList(const char*) = 0;
	virtual void gviewplotPair(const char *) = 0;
};

//derived 2D-class for collision detection and handling
class CollisionSolver2d : public CollisionSolver {
private:
	//input bonds
	std::vector<BOND*> bondsList;
	//bond pairs
	std::vector<BOND_PAIR> bondPairList;
public:
	void assembleFromInterface(const INTERFACE*);
	void detectProximity();
	void detectCollision();
	void resolveCollision();
	void printProximity();
	void printCollision();
	void gviewplotPairList(const char*);
	void gviewplotPair(const char *);
	void getBondPairList(std::vector<std::pair<BOND*,BOND*> >&);
};

//derived 3D-class for collision detection and handling
class CollisionSolver3d : public CollisionSolver {
private:
	//input tris
	std::vector<TRI*> trisList;
	//tri pairs
	std::vector<TRI_PAIR> triPairList;
	void computeAverageVelocity();
	void updateAverageVelocity();
	void updateFinalVelocity();
	void updateFinalPosition();
	void recordOriginVelocity();
public:
	void assembleFromInterface(const INTERFACE*);
	void detectProximity();
	void detectCollision();
	void resolveCollision();
	void printProximity();
	void printCollision();
	void gviewplotPairList(const char*);
	void gviewplotPair(const char *);
	void getTriPairList(std::vector<TRI_PAIR>&);
};
