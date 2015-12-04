//collid3d.h
#include <CGAL/Simple_cartesian.h>
#include <CGAL/box_intersection_d.h>

typedef std::pair<TRI*,TRI*> TRI_PAIR;
typedef std::pair<BOND*,BOND*> BOND_PAIR;

struct STATE{
	double vel[3];
	double impulse[3];
	double avgVel[3];
};
//box traits structure for FronTier triangles
template<typename T> struct Traits_of_FT
{
	typedef double          NT;
        typedef T*            Box_parameter;
        typedef std::ptrdiff_t ID;
};

template<> struct Traits_of_FT<TRI*>{
	typedef double 		NT;
	typedef TRI*		Box_parameter;
	typedef std::ptrdiff_t ID;

	static double eps;
	static int dimension(){ return 3;}
	static double min_coord(Box_parameter b, int d){
		return std::min(std::min(Coords(b->__pts[0])[d],
			Coords(b->__pts[1])[d]),
			Coords(b->__pts[2])[d])-eps;	
	}
	static double max_coord(Box_parameter b, int d){
		return std::max(std::max(Coords(b->__pts[0])[d],
                        Coords(b->__pts[1])[d]),
                        Coords(b->__pts[2])[d])+eps;
	}
	static std::ptrdiff_t id(Box_parameter b) { return (std::ptrdiff_t)(b);}
};

//box traits structure for FronTier tri_pair
template<> struct Traits_of_FT<std::pair<TRI*,TRI*> >{
	typedef double          		NT;
        typedef std::pair<TRI*,TRI*>      Box_parameter;

	static double eps;
        static int dimension(){ return 3;}
        static double min_coord(Box_parameter b, int d){
		double ans = INT_MAX;
		for (int i = 0; i < 3; ++i){
			ans = std::min(ans,Coords((b.first)->__pts[i])[d]);
			ans = std::min(ans,Coords((b.second)->__pts[i])[d]);
		}
		return ans-eps;
        }
        static double max_coord(Box_parameter b, int d){
		double ans = INT_MIN;
		for (int i = 0; i < 3; ++i){
			ans = std::max(ans,Coords((b.first)->__pts[i])[d]);
			ans = std::max(ans,Coords((b.second)->__pts[i])[d]);
		}
		return ans+eps;
        }
	static std::ptrdiff_t id(Box_parameter b) { return (std::ptrdiff_t)(b.first);}
};

//box traits structure for FronTier bonds
template<> struct Traits_of_FT<BOND*>{
	typedef double          NT;
        typedef BOND*            Box_parameter;
        typedef std::ptrdiff_t ID;

	static double eps;
        static int dimension(){ return 2;}
        static double min_coord(Box_parameter b, int d){
                return std::min(Coords(b->start)[d],
                        Coords(b->end)[d])-eps;
        }
        static double max_coord(Box_parameter b, int d){
                return std::max(Coords(b->start)[d],
                        Coords(b->end)[d])+eps;
        }
	static std::ptrdiff_t id(Box_parameter b) { return (std::ptrdiff_t)(b);}
};

//box traits structure for FronTier bonds_pair
template<> struct Traits_of_FT<std::pair<BOND*,BOND*> >{
	typedef double          		NT;
        typedef std::pair<BOND*,BOND*>     Box_parameter;

	static double eps;
        static int dimension(){ return 2;}
        static double min_coord(Box_parameter b, int d){
                double ans =  std::min(Coords(b.first->start)[d],
                        Coords(b.first->end)[d]);
		return std::min(ans,std::min(Coords(b.second->start)[d],
                        Coords(b.second->end)[d]))-eps;
        }
        static double max_coord(Box_parameter b, int d){
                double ans =  std::max(Coords(b.first->start)[d],
                        Coords(b.first->end)[d]);
		return std::max(ans,std::max(Coords(b.second->start)[d],
                        Coords(b.second->end)[d]))+eps;
        }
	static std::ptrdiff_t id(Box_parameter b) { return (std::ptrdiff_t)(b.first);}
};

extern bool doIntersect(const TRI*,const TRI*);
extern bool doIntersect(const TRI_PAIR&,const TRI_PAIR&);

template <typename T>
struct Report {
    std::vector<std::pair<T,T> > *output;
    Report(std::vector<std::pair<T,T> > &pairList) : output(&pairList) {}
    // We write the elements with respect to 'boxes' to the output
    void operator()( const T &a, const T &b) {
	if (doIntersect(a,b))
    	    output->push_back(std::make_pair<T,T>(a,b));
    }
};

//abstract base class for collision detection and handling
class CollisionSolver {
private:
	//rounding tolerance
	static double eps;
	static double thickness;
protected:
	Front *front;
public:
	static void setRoundingTolerance(double);
	static double getRoundingTolerance();
	static void setFabricThickness(double);
	static double getFabricThickness();
	virtual ~CollisionSolver() {}; //virtual destructor
	//pure virtual functions
	virtual void assembleFromInterface(const INTERFACE*) = 0;
	virtual void assembleFromTwoInterface(const INTERFACE*,const INTERFACE*) = 0;
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
	//pair of bond pairs
	std::vector<std::pair<BOND_PAIR,BOND_PAIR> >bondTwoPairList;
public:
	void assembleFromInterface(const INTERFACE*);
	void assembleFromTwoInterface(const INTERFACE*,const INTERFACE*);
	void detectProximity();
	void detectCollision();
	void resolveCollision();
	void printProximity();
	void printCollision();
	void gviewplotPairList(const char*);
	void gviewplotPair(const char *);
	void getProximityBondsPairList(std::vector<std::pair<BOND*,BOND*> >&);
	void getCollisionBondsDoublePairList(std::vector<std::pair<BOND_PAIR,BOND_PAIR> >&);
};

//derived 3D-class for collision detection and handling
class CollisionSolver3d : public CollisionSolver {
private:
	//input tris
	std::vector<TRI*> trisList;
	//tri pairs
	std::vector<TRI_PAIR> triPairList;
	//pair of tri pairs
	std::vector<std::pair<TRI_PAIR,TRI_PAIR> > triTwoPairList;
public:
	void assembleFromInterface(const INTERFACE*);
	void assembleFromTwoInterface(const INTERFACE*,const INTERFACE*);
	void detectProximity();
	void detectCollision();
	void resolveCollision();
	void printProximity();
	void printCollision();
	void gviewplotPairList(const char*);
	void gviewplotPair(const char *);
	void getProximityTrisPairList(std::vector<std::pair<TRI*,TRI*> >&);
	void getCollisionTrisDoublePairList(std::vector<std::pair<TRI_PAIR,TRI_PAIR> >&);
};
