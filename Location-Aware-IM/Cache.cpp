#include "limit.h"
#include "Cache.h"

#include <stdio.h>
#include <string.h>
#include <vector>
#include <map>
#include <cassert>

vector<vector<int> > Cache::topTau;
vector<vector<int> > Cache::cdds;
vector<vector<double> > Cache::dp;
vector<vector<double> > Cache::dpTau;
vector<bool> Cache::is_init_cached;
vector<bool> Cache::is_tau_cached;

void Cache::init(int numOfRegions)
{
    cdds.resize(numOfRegions);
    topTau.resize(numOfRegions);
    dp.resize(numOfRegions);
    dpTau.resize(numOfRegions);
    is_init_cached.resize(numOfRegions, false);
    is_tau_cached.resize(numOfRegions, false);
    for (int _ = 0; _ < numOfRegions; _++) {
        cdds[_].clear();
        topTau[_].clear();
        dpTau[_].clear();
        dp[_].clear();
    }
}

bool Cache::isInitCached(int rid)
{  // if rid region has been cached
    return is_init_cached[rid];
}

bool Cache::isTauCached(int rid)
{  // if rid region has been cached
    return is_tau_cached[rid];
}

void Cache::setInitCached(int rid, bool isCached)
{
    is_init_cached[rid] = isCached;
}

void Cache::setTauCached(int rid, bool isCached)
{
    is_tau_cached[rid] = isCached;
}

bool Cache::getTopTau(vector<int>& list, int rid, int num)
{   // return false if this region is not cached
    list.clear();
    assert(rid >= 0 && rid < topTau.size());

    if (num < 0 || num > topTau[rid].size()) {
        num = topTau[rid].size();
    }
    list.resize(num);
    for (int i = 0; i < num; i++) {
        list[i] = topTau[rid][i];
    }

    if (list.size() == 0)
        return false;
    else 
        return true;
}

bool Cache::getCandidates(vector<int>& list, int rid, int num)
{   // return false if this region is not cached
    list.clear();
    assert(rid >= 0 && rid < cdds.size());

    if (num < 0 || num > cdds[rid].size()) {
        num = cdds[rid].size();
    }
    list.resize(num);
    for (int i = 0; i < num; i++) {
        list[i] = cdds[rid][i];
    }

    if (list.size() == 0)
        return false;
    else 
        return true;
}

bool Cache::getDP(vector<double>& list, int rid, int num)
{   // return false if this region is not cached
    list.clear();
    assert(rid >= 0 && rid < dp.size());

    if (num < 0 || num > dp[rid].size()) {
        num = dp[rid].size();
    }
    list.resize(num);
    for (int i = 0; i < num; i++) {
        list[i] = dp[rid][i];
    }

    if (list.size() == 0)
        return false;
    else 
        return true;
}

bool Cache::getDPTau(vector<double>& list, int rid, int num)
{   // return false if this region is not cached
    list.clear();
    if (num < 0 || num > dpTau[rid].size()) {
        num = dpTau[rid].size();
    }
    list.resize(num);
    for (int i = 0; i < num; i++) {
        list[i] = dpTau[rid][i];
    }

    if (list.size() == 0)
        return false;
    else 
        return true;
}

void Cache::setCandidates(int rid, vector<int> list)
{
    // rid must be a valid region ID
    assert(rid >= 0 && rid < cdds.size());
    cdds[rid] = list;
    is_init_cached[rid] = true;
}

void Cache::setTopTau(int rid, vector<int> list)
{
    assert(rid >= 0 && rid < topTau.size());
    topTau[rid] = list;
    is_tau_cached[rid] = true;
}

void Cache::setDP(int rid, vector<double> list)
{
    assert(rid >= 0 && rid < dp.size());
    dp[rid] = list;
}

void Cache::setDPTau(int rid, vector<double> list)
{
    assert(rid >= 0 && rid < dpTau.size());
    dpTau[rid] = list;
}
