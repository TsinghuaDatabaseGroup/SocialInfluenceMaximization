#ifndef __QUADTREE__
#define __QUADTREE__

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>

using namespace std;

struct XY {
    XY(float _x, float _y) : x(_x), y(_y) {}  // x as longitude, y as latitude
    float x, y;
};

struct AABB {
    AABB(XY _center, XY _halfDimension) : center(_center), halfDimension(_halfDimension) {}

    bool containsPoint(float x, float y) {
        XY p(x, y);
        return ((p.x >= center.x - halfDimension.x && p.x <= center.x + halfDimension.x) &&
                (p.y >= center.y - halfDimension.y && p.y <= center.y + halfDimension.y));
    }

    bool intersectsAABB(AABB other) {
        // intersect or contain
        return fabs(other.center.x - center.x) < (halfDimension.x + other.halfDimension.x) && (fabs(other.center.y - center.y) < (halfDimension.y + other.halfDimension.y));
    }

    bool containsAABB(AABB other) {
        return (
        center.x + halfDimension.x > other.center.x + other.halfDimension.x && 
        center.x - halfDimension.x < other.center.x - other.halfDimension.x && 
        center.y + halfDimension.y > other.center.y + other.halfDimension.y && 
        center.y - halfDimension.y < other.center.y - other.halfDimension.y);
    }

    XY center;
    XY halfDimension;
};

class Quadtree {
    public:
        static int count;

        //  Quadtree(XY, XY);
        Quadtree(XY, XY, int);

        bool insert(int, float, float);
        void subdivide();
        void queryRange(std::vector<int> &, AABB);
        void queryRangeFO(std::vector<int> &, std::vector<int> &, vector<AABB>, AABB);
        void queryRangeLocations(std::vector<int> &, std::vector<float> &, std::vector<float> &, AABB);

        // queryRangeRegion should contain all regions that intersects with AABB
        // not necessarily contained by AABB
        void queryRangeRegion(std::vector<int> &, std::vector<AABB> &, AABB);
        void overall(std::vector<int> &);
        void print();

        void clear();

        AABB boundary;
        int nodeCapacity;
        int id;
        int population;

        // leaves
        Quadtree * NW;
        Quadtree * NE;
        Quadtree * SW;
        Quadtree * SE;

        // data
        std::vector<int> ids;
        std::vector<float> xs;
        std::vector<float> ys;
};

#endif // __QUADTREE__
