#include <cassert>
#include "Quadtree.h"

// Quadtree::Quadtree(XY _center, XY _halfDimension)
//   : boundary(_center, _halfDimension), nodeCapacity(4)
// {
//   NW = NE = SW = SE = NULL;
//   points.reserve(nodeCapacity);
// }

int Quadtree::count = 0;

Quadtree::Quadtree(XY _center, XY _halfDimension, int _nodeCapacity)
: boundary(_center, _halfDimension), nodeCapacity(_nodeCapacity)
{
    id = count++;
    NW = NE = SW = SE = NULL;
    population = 0;
    ids.reserve(nodeCapacity);
    xs.reserve(nodeCapacity);
    ys.reserve(nodeCapacity);
}

bool Quadtree::insert(int id, float x, float y) {
    if (!boundary.containsPoint(x, y))
        return false;

    if (ids.size() < nodeCapacity) {
        ids.push_back(id);
        xs.push_back(x);
        ys.push_back(y);
        population += 1;
        return true;
    }

    if (NW == NULL) subdivide();

    if (NW->insert(id, x, y) || 
        NE->insert(id, x, y) ||
        SW->insert(id, x, y) ||
        SE->insert(id, x, y)) {
        population += 1;
        return true;
    }

    assert(false);
    return false; // should never happen
}

void Quadtree::subdivide() {
    XY center = boundary.center;
    XY newDim(boundary.halfDimension.x / 2, boundary.halfDimension.y / 2);

    NW = new Quadtree(XY(center.x - newDim.x, center.y + newDim.y), newDim, nodeCapacity);
    NE = new Quadtree(XY(center.x + newDim.x, center.y + newDim.y), newDim, nodeCapacity);
    SW = new Quadtree(XY(center.x - newDim.x, center.y - newDim.y), newDim, nodeCapacity);
    SE = new Quadtree(XY(center.x + newDim.x, center.y - newDim.y), newDim, nodeCapacity);
}

void Quadtree::queryRange(std::vector<int> & list, AABB range) {
    if (!boundary.intersectsAABB(range)) return ; // list is empty

    for (int i = 0; i < ids.size(); ++i)
        if (range.containsPoint(xs[i], ys[i])) {
            list.push_back(ids[i]);
        }

    if (NW == NULL) return ;
    NW->queryRange(list, range);
    NE->queryRange(list, range);
    SW->queryRange(list, range);
    SE->queryRange(list, range);
}

void Quadtree::queryRangeLocations(std::vector<int> & list, vector<float>& nativesX, vector<float>& nativesY, AABB range) {
    if (!boundary.intersectsAABB(range)) return ; // list is empty

    for (int i = 0; i < ids.size(); ++i)
        if (range.containsPoint(xs[i], ys[i])) {
            list.push_back(ids[i]);
            nativesX.push_back(xs[i]);
            nativesY.push_back(ys[i]);
        }

    if (NW == NULL) return ;
    NW->queryRangeLocations(list, nativesX, nativesY, range);
    NE->queryRangeLocations(list, nativesX, nativesY, range);
    SW->queryRangeLocations(list, nativesX, nativesY, range);
    SE->queryRangeLocations(list, nativesX, nativesY, range);
}

void Quadtree::queryRangeRegion(std::vector<int> & list, std::vector<AABB> & boundaries, AABB range)
{
//    printf("XY:(%f %f), halfDimension:(%f %f)\n", boundary.center.x, boundary.center.y, boundary.halfDimension.x, boundary.halfDimension.y);
//    for (int i = 0; i < ids.size(); i++) {
//        printf("%d \n", ids[i]);
//    }
    if (!boundary.intersectsAABB(range)) return;
//    if (population > MAX_POPULATION) {
//        // too large, ask subregions directly.
//        if (NW != NULL) {
//            NW->queryRangeRegion(list, boundaries, populations, range);
//            NE->queryRangeRegion(list, boundaries, populations, range);
//            SW->queryRangeRegion(list, boundaries, populations, range);
//            SE->queryRangeRegion(list, boundaries, populations, range);
//        }
//    }
    else {  // moderate region, have intersection with query range
        if (range.containsAABB(boundary)) {  // totally covered by query range
            list.push_back(id);
            boundaries.push_back(boundary);
//            printf("XY:(%f %f), halfDimension:(%f %f)\n", boundary.center.x, boundary.center.y, boundary.halfDimension.x, boundary.halfDimension.y);
        }
        else if (NW != NULL) {
            // not covered by query range, but it is a intermediate region.
            NW->queryRangeRegion(list, boundaries, range);
            NE->queryRangeRegion(list, boundaries, range);
            SW->queryRangeRegion(list, boundaries, range);
            SE->queryRangeRegion(list, boundaries, range);
        }
    }
}

void Quadtree::queryRangeFO(std::vector<int> &list, std::vector<int> &others, vector<AABB> boundaries, AABB range)
{
    if (!boundary.intersectsAABB(range)) return ; // list is empty

    int numOfRegions = boundaries.size();
    for (int i = 0; i < ids.size(); ++i)
        if (range.containsPoint(xs[i], ys[i])) {
            list.push_back(ids[i]);
            int j = 0;
            for (; j < numOfRegions; j++)
                if (boundaries[j].containsPoint(xs[i], ys[i]))
                    break;
            if (j == numOfRegions)
                others.push_back(ids[i]);
        }

    if (NW == NULL) return ;
    NW->queryRange(list, range);
    NE->queryRange(list, range);
    SW->queryRange(list, range);
    SE->queryRange(list, range);
}

void Quadtree::overall(std::vector<int> & list)
{
    for (int i = 0; i < ids.size(); ++i) {
        list.push_back(ids[i]);
    }
    if (NW == NULL) return ;
    NW->overall(list);
    NE->overall(list);
    SW->overall(list);
    SE->overall(list);
}

void Quadtree::print()
{
    printf("XY:(%f %f), halfDimension:(%f %f)\n", boundary.center.x, boundary.center.y, boundary.halfDimension.x, boundary.halfDimension.y);
    for (int i = 0; i < ids.size(); ++i) {
        printf("%d ", ids[i]);
    }
    printf("\n");
    if (NW == NULL) return ;
    NW->print();
    NE->print();
    SW->print();
    SE->print();
}

// scan the tree and remove all node/Item*
void Quadtree::clear() {
    count = 0;
    if (this == NULL) return ;
    ids.clear();
    ys.clear();
    xs.clear();

    NW->clear();
    delete NW;
    NW = NULL;
    NE->clear();
    delete NE;
    NE = NULL;
    SW->clear();
    delete SW;
    SW = NULL;
    SE->clear();
    delete SE;
    SE = NULL;
}
