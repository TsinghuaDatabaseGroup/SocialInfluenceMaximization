#ifndef CACHE_H
#define CACHE_H

#include <vector>
#include <map>

using namespace std;

class Cache
{
    private:
        static vector<vector<int> >topTau;  // set of local topk users
        static vector<vector<int> >cdds;  // set of candidates in a Quadtree node
        // initial incremental influence of each candidate
        static vector<vector<double> >dp; 
        static vector<vector<double> >dpTau; 

        static vector<bool> is_init_cached;
        static vector<bool> is_tau_cached;

    public:
        static void init(int numOfRegions);
        static bool isInitCached(int rid);
        static bool isTauCached(int rid);
        static void setInitCached(int rid, bool isCached);
        static void setTauCached(int rid, bool isCached);

        static bool getTopTau(vector<int>& list, int rid, int num);
        static bool getCandidates(vector<int>& list, int rid, int num);
        static bool getDP(vector<double>& list, int rid, int num);
        static bool getDPTau(vector<double>& list, int rid, int num);

        static void setCandidates(int rid, vector<int> list);
        static void setTopTau(int rid, vector<int> list);
        static void setDP(int rid, vector<double> list);
        static void setDPTau(int rid, vector<double> list);
};

#endif
