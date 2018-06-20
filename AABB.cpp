#include "AABB.h"
#include <stack>
#include <fstream>

//  for proximity detection
AABB::AABB(CD_HSE* h, aabb::type type) : hse(h), dt(0), abType(type), 
        lowerbound(3), upperbound(3) {
    if (type == aabb::type::STATIC)
        for (int i = 0; i < 3; i++) {
             lowerbound[i] = h->min_static_coord(i)-1e-6;
             upperbound[i] = h->max_static_coord(i)+1e-6;
        }
    else
        for (int i = 0; i < 3; i++) {
             lowerbound[i] = h->min_moving_coord(i, dt)-1e-6;
             upperbound[i] = h->max_moving_coord(i, dt)+1e-6;
        }
    for (int i = 0; i < h->num_pts(); i++)
         indices.push_back(h->Point_of_hse(i)->global_index);
}
// for collision detection
AABB::AABB(CD_HSE* h, aabb::type type, double t) : hse(h), dt(t), 
        abType(type), lowerbound(3), upperbound(3) {   
    if (type == aabb::type::STATIC) 
        for (int i = 0; i < 3; i++) {
             lowerbound[i] = h->min_static_coord(i)-1e-6;
             upperbound[i] = h->max_static_coord(i)+1e-6;
        }
    else 
        for (int i = 0; i < 3; i++) {
             lowerbound[i] = h->min_moving_coord(i, dt)-1e-6;
             upperbound[i] = h->max_moving_coord(i, dt)+1e-6;
        }
    for (int i = 0; i < h->num_pts(); i++) 
         indices.push_back(h->Point_of_hse(i)->global_index);
}

AABB::AABB(const CPoint& pl, const CPoint& pu) : lowerbound(pl), upperbound(pu) {
}

AABB AABB::merge(const AABB& ab) const {
    CPoint pl(3), pu(3);

    for (int i = 0; i < 3; i++) {
         pl[i] = std::min(lowerbound[i], ab.lowerbound[i]);
         pu[i] = std::max(upperbound[i], ab.upperbound[i]);
    }
    return AABB(pl, pu);
} 

double AABB::volume() {
    return (upperbound[0]-lowerbound[0])*(upperbound[1]-lowerbound[1])*
            (upperbound[2]-lowerbound[2]);
}

bool AABB::isCollid(const AABB& ab) {
    return (lowerbound[0] <= ab.upperbound[0] && upperbound[0] >= ab.lowerbound[0]) && 
           (lowerbound[1] <= ab.upperbound[1] && upperbound[1] >= ab.lowerbound[1]) && 
           (lowerbound[2] <= ab.upperbound[2] && upperbound[2] >= ab.lowerbound[2]); 
}

void AABB::updateAABBInfo(double dt) {
    if (abType == aabb::type::STATIC)
        for (int i = 0; i < 3; i++) {
             lowerbound[i] = hse->min_static_coord(i)-1e-6;
             upperbound[i] = hse->max_static_coord(i)+1e-6;
        }
    else 
        for (int i = 0; i < 3; i++) {
             lowerbound[i] = hse->min_moving_coord(i, dt)-1e-6;
             upperbound[i] = hse->max_moving_coord(i, dt)+1e-6;
        }
}

bool AABB::contain(const AABB* ab) {
    return lowerbound[0] <= ab->lowerbound[0] && lowerbound[1] <= ab->lowerbound[1] &&
        lowerbound[2] <= ab->lowerbound[2] && upperbound[0] >= ab->upperbound[0] &&
        upperbound[1] >= ab->upperbound[1] && upperbound[2] >= ab->upperbound[2];
}

void Node::setBranch(Node* n1, Node* n2) {
    n1->parent = this;
    n2->parent = this;
    left = n1;
    right = n2;
}

bool Node::isLeaf() {
    return left == nullptr && right == nullptr;
}

void Node::setLeaf(AABB* ab) {
    data = ab;
}

void Node::updateBranch() {
    if (isLeaf())
        return;
    for (int i = 0; i < 3; i++) {
         box.lowerbound[i] = std::min(left->box.lowerbound[i], right->box.lowerbound[i]);
         box.upperbound[i] = std::max(left->box.upperbound[i], right->box.upperbound[i]);
    }
}

void Node::updateAABB() {
    if (isLeaf()) {
        box.lowerbound = data->lowerbound;
        box.upperbound = data->upperbound;
    }
    else {
        // branch node has no AABB yet
        if (box.lowerbound.size() == 0)
            box = left->box.merge(right->box);
        else 
            updateBranch();
    }
}

bool Node::isCollid(Node* n) {
    return box.isCollid(n->box);
}

Node* Node::getSibling() const {
    if (!this->parent) 
        return nullptr;
    return this == this->parent->left ? this->parent->right : this->parent->left;
}

void AABBTree::addAABB(AABB* ab) {
    if (root) {
        Node* node = new Node;

        node->setLeaf(ab);
        node->updateAABB();
        insertNode(node, root);
        nodeArray.push_back(node);
        numLeaf++;
    }
    else {
        root = new Node;
        root->setLeaf(ab);
        root->updateAABB();
        nodeArray.push_back(root);
        numLeaf++;
    }
}

// reorganize the tree structure
void AABBTree::updateTreeStructure() {
    for (int i = numLeaf; i< nodeArray.size(); i++) {
         delete nodeArray[i];
         nodeArray[i] = nullptr;
    }
    root = nullptr;
    for (auto node : nodeArray) {
         if (root) 
             insertNode(node, root);
         else 
             root = node;
    }
}

void AABBTree::insertNode(Node* n, Node*& parent) {
    Node* p = parent;
    // if parent is a leaf node, then create a branch
    // with n and parent to be two children
    if (p->isLeaf()) {
        Node* newParent = new Node;

        
        newParent->parent = p->parent;
        if (p->parent)
            p->parent->left == p? p->parent->left = newParent : 
                p->parent->right = newParent;
        newParent->setBranch(n, p);
        parent = newParent;
    }
    // we have to decide which subtree to insert to
    // the rule is insert to the subtree with smaller volume 
    else {
        AABB& abl = p->left->box;
        AABB& abr = p->right->box;
        // get volume after inserting current node to 
        // left or right subtree
        double vdiff1 = abl.merge(n->box).volume()-abl.volume();
        double vdiff2 = abr.merge(n->box).volume()-abr.volume();
        // insert to left subtree
        if (vdiff1 < vdiff2) {
            insertNode(n, p->left);
        }
        else
            insertNode(n, p->right);
    }
    // this will guarantee all relavent ancestor will be 
    // updated
    parent->updateAABB();
}

void AABBTree::updatePointMap(const std::vector<CD_HSE*>& hseList) {
    vhMap.clear();
    nodeSet.clear();
    count = 0; 
    for (auto it : hseList) {
         std::vector<long> ids;
         for (int i = 0; i < it->num_pts(); i++) 
              ids.push_back(it->Point_of_hse(i)->global_index);
         vhMap.insert({ids, it});
    }
    
}

void AABBTree::updateAABBTree(const std::vector<CD_HSE*>& hseList) {
    updatePointMap(hseList);

    std::stack<Node*> sn;
    Node* cur = root;

    if (!root) return;
    // iterative postorder traverse
    do {
        while (cur) {
            if (cur->right)
                sn.push(cur->right);
            sn.push(cur);
            cur = cur->left;
        }
        cur = sn.top();
        sn.pop();
        if (cur->right && !sn.empty() && cur->right == sn.top()) {
            sn.pop();
            sn.push(cur);
            cur = cur->right;
        }
        else {
            if (cur->isLeaf()) {
                cur->data->hse = vhMap[cur->data->indices];
                cur->data->updateAABBInfo(dt);
                cur->updateAABB();    
            }
            if (!cur->isLeaf())
                cur->updateBranch();
            cur = nullptr;
        }
    } while (!sn.empty());
}

double AABBTree::treeHeight(Node* root) {
    if (!root)
        return 0;
    return std::max(treeHeight(root->left), treeHeight(root->right))+1;
}
// inorder traverse the tree and whenever come up with a leaf node, 
// find collided pairs correspond to it.
void AABBTree::query(CollisionSolver* collsn_solver, int type) {

    Node* cur = root;
    std::stack<Node*> sn;

    while (cur || !sn.empty()) {
        while (cur) {
            sn.push(cur);
            cur = cur->left;
        }
        cur = sn.top();
        sn.pop();
        
        if (cur->isLeaf()) {
            if (type == aabb::type::STATIC)
                queryProximity(cur, collsn_solver);
            else
                isCollsn = queryCollision(cur, collsn_solver);
            nodeSet.insert(cur);
        }
        
        cur = cur->right;
    }
}

// for AABB inside Node n, find all collied AABBs
// and add it to the list
// preorder traverse the tree and if find a collided node to be 
// (1) leaf, find a pair and add to the list
// (2) branch, push two children into the stack
void AABBTree::queryProximity(Node* n, CollisionSolver* collsn_solver) {
    std::stack<Node*> sn;
    Node* cur = root;

    while (cur || !sn.empty()) {
        while (cur) {
            if (cur->isCollid(n)) {
                if (cur->isLeaf() && n != cur) {
                    if (nodeSet.find(cur) == nodeSet.end()) {
                        CD_HSE* a = cur->data->hse;
                        CD_HSE* b = n->data->hse;

                        if (collsn_solver->isProximity(a, b))
                            count++; 
                    }
                }
                sn.push(cur);
                cur = cur->left;
            }   
            else 
                break;    
        }
        if (sn.empty())
            break;
        cur = sn.top();
        sn.pop();
        cur = cur->right;
    }
}

bool AABBTree::queryCollision(Node* n, CollisionSolver* collsn_solver) {
    std::stack<Node*> sn;
    Node* cur = root;

    while (cur || !sn.empty()) {
        while (cur) {
            if (cur->isCollid(n)) {
                if (cur->isLeaf() && n != cur) {
                    if (nodeSet.find(cur) == nodeSet.end()) {
                        CD_HSE* a = cur->data->hse;
                        CD_HSE* b = n->data->hse;

                        if (collsn_solver->isCollision(a, b)) 
                            count++;
                    }
                }
                sn.push(cur);
                cur = cur->left;
            }   
            else 
                break;    
        }
        if (sn.empty())
            break;
        cur = sn.top();
        sn.pop();
        cur = cur->right;
    }
    return count > 0;
}
const CollidPairList& AABBTree::getCollidPair() {
    return colldList;
}
