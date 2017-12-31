#ifndef MIP_H
#define MIP_H

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include "limit.h"
#include "graph.h"
#include "heap.h"

using namespace std;

class Mip
{
    private:
	    static int n;
        // maximum influence path propagation probability.
        // resize to n
        static vector<double *> mip;
        static vector<int *> inflees;
        static vector<int> numInflees;

        static vector<bool> changed;
        static vector<double *> bmip;
        static vector<int *> bInflees;
        static vector<int> bNumInflees;

        static vector<int *> outPath;
        static vector<double *> outAp;
        static vector<double *> outB;

        static vector<vector<int> > cdds;
        static vector<int> numCdds;

        static Heap *maxHeap;  // \Heap in paper
        static Heap *lMaxHeap;  // lower bound \Heap in paper

        // seeds currently selected
	    static vector<bool> used;  // selected as topk, global structure

        // activation probability
        // resize to n
        static vector<double> ap;
        static vector<double> sap;

        // paths from influencees to seeds.
        // resize to n
        static vector<vector<int> > path;
        static vector<vector<int> > children;
        static vector<vector<double> > bb;
        static vector<double> b;
        static vector<int> numInfrs;  // total number of influencers
        static vector<int> numchild; 

        // dijkstra
	    static vector<int> S; 
        static vector<int> parent;
        static vector<int> queue;
	    static vector<int> heap;
        static vector<double> distance;
        static vector<int> childlist;

        static vector<int> status;
        static double longest;

        // hints
//        static vector<double> hap;
//        static vector<double> uap;  // upper bound ap
//        static double lowerHint;

        // helper functions
        // spread from one node
        static void shortestPathsFrom(int fromNode);

        // really incremental influence
        static double incInfl(int s);
        // use \Heap to return next best seed in \Heap
        static int nextSeed(double *inc);
        // accept s as seed. It means to update influencees' influence path to
        // include s in them.
        static void accept(int s);
        static void gc();
        static double estimateIncInfl(int s);
        static double estimateLowerIncInfl(const vector<int>&);
        static double estimateUpperIncInfl(const vector<int>&);
        static double lowerIncInfl(int s);

        static double getHints(vector<int>&hint, const vector<vector<int> >& topTau, const vector<vector<double> >& dpTau, int num);
    
    public:
        // initialization process, calculate all MIPs
        static void init(double longest);

        // influence spread function
//        static double spread(const unordered_set<int>& kset, const unordered_set<int>& influencees);

        // select
        static vector<int> select(int k);
        static void cacheTau(int k, vector<int>& topTau, vector<double>& tauDp);
        static vector<int> assembly(int k, const vector<vector<int> >& candidates, const vector<vector<double> > &initDp);
        static vector<int> hint(int k, double acceptRatio, const vector<vector<int> >& candidates, const vector<vector<double> > &initDp, const vector<vector<int> >& topTau, const vector<vector<double> > &dpTau);
        static vector<int> dumbHint(int k, double acceptRatio, const vector<vector<int> >& candidates, const vector<vector<double> > &initDp);
        static void reset_native();
        static void exit();
        static void sortInit(const vector<int>& vertices, vector<int>& candidates,
                vector<double>& cddsDp);

        // helper variables
        static int count;
        static int estCount;
        static double incInflTimer;
        static vector<bool> is_native;  // if this node is a native user.
        static void reset();
};

#endif
