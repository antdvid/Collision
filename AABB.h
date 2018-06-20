#ifndef AABB_H_
#define AABB_H_

#include <vector>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include "collid.h"

// header file for AABB tree
// axis-aligned bounding box class
using CPoint = std::vector<double>;

namespace aabb {
    // type is used to separate proximity and collition check
    enum type{STATIC, MOVING};
}

class Node;
class AABBTree;

class AABB {
    friend class Node;
    friend class AABBTree;
    CPoint lowerbound;
    CPoint upperbound;
    // indices will store the index of points on the  
    // corresponding triangle or bond.
    std::vector<long> indices;
    void updateAABBInfo(double);
    //void updateAABBInfo(const std::unordered_map<long, POINT*>&);
    bool contain(const AABB*);
    CD_HSE* hse = nullptr;
    double dt;
    aabb::type abType;
public:
    // constructor
    AABB() {}
    AABB(CD_HSE*, aabb::type);
    AABB(CD_HSE*, aabb::type, double);
    AABB(const CPoint&, const CPoint&);
    // merge this with anther AABB to get a
    // merged AABB and construct the corresponding AABB tree
    AABB merge(const AABB&) const;
    // get the volume of the AABB
    double volume();
    bool isCollid(const AABB&);
};

// tree node corresponding to AABB
class AABBTree;

class Node {
    friend class AABBTree;
    // AABB stored in node. May store information for branch AABB
    // and may be adjusted for dynamic AABB 
    AABB box;
    // if leaf, point to the corresponding AABB
    // empty for branch
    AABB* data;
    // parent node
    Node* parent = nullptr;
    // left and right children node
    Node* left = nullptr;
    Node* right = nullptr;
    void updateBranch();
public:
    // make this node to be brance from two Node parameter
    void setBranch(Node*, Node*);
    // judge if this node is a leaf
    bool isLeaf();
    // set an AABB element to be a leaf
    void setLeaf(AABB*);
    bool isCollid(Node*);
    void updateAABB();
    Node* getSibling() const;
};

using CollidPairList = std::vector<std::pair<AABB*, AABB*>>;
using CollidPairSet = std::set<std::pair<Node*, Node*>>;

class AABBTree {
    Node* root = nullptr;
    // store the collied pair
    CollidPairList colldList;
    // used to check if certain pair is considered or not.
    CollidPairSet collidSet;
    // node needed to be removed and reinsert to the tree
    std::unordered_map<long, POINT*> ump;
    // map from object's indices (2 or 3 points' global indices)
    // to corresponding CD_HSE* in collision library 
    std::map<std::vector<long>, CD_HSE*> vhMap;
    std::unordered_set<Node*> nodeSet;
    std::vector<Node*> nodeArray;
    int count;
    int numLeaf = 0;
    double treeHeight(Node*); 
    double dt;
    bool isCollsn;
    void queryProximity(Node*, CollisionSolver*);
    bool queryCollision(Node*, CollisionSolver*);
public:
    // add an AABB element into a tree
    void addAABB(AABB*);
    // insert a node into the subtree with parent 
    // as the root
    void insertNode(Node*, Node*&);
    // query all collid pairs
    void query(CollisionSolver*, int);
    // check all collid elements for node parameter
    const CollidPairList& getCollidPair();
    int getCount() { return count; }
    double getVolume() { return root->box.volume(); } 
    void updateTreeStructure();
    void updatePointMap(const std::vector<CD_HSE*>&);
    void setTimeStep(double t) { dt = t; }
    bool getCollsnState() { return isCollsn; }
    void updateAABBTree(const std::vector<CD_HSE*>&);
};

#endif
