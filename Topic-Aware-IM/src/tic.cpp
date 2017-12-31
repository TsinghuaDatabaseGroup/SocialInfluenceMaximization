#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>     /* assert */
#include <algorithm>    // std::make_heap, std::pop_heap, std::push_heap, std::sort_heap
#include <set>
#include <map>
#include <time.h>

#include "tic.h"

using namespace std;

int TIC::n_ = 0;
int TIC::m_ = 0;
int TIC::Z_ = 0;
//vector<int> TIC::indegree_;
//vector<int> TIC::outdegree_;
vector<int> TIC::inindex_;
vector<int> TIC::outindex_;
vector<Edge> TIC::inedges_;
vector<Edge> TIC::outedges_;
vector<double> TIC::incache_;
vector<double> TIC::outcache_;
vector<bool> TIC::outcflag_;  // is cached ?
vector<bool> TIC::incflag_;
map<int, int> TIC::cvt_;

typedef std::pair<double, int> di;

int TIC::GetN() {
    return n_;
}

int TIC::GetM() {
    return m_;
}

int TIC::GetZ() {
    return Z_;
}

void TIC::resetCache() {
    for (int i = 0; i < m_; i++) {
        outcflag_[i] = false;
        incflag_[i] = false;
    }
}

void TIC::init() {
    n_ = m_ = Z_ = 0;
    cvt_.clear();
    inindex_.clear();
    outindex_.clear();
    inedges_.clear();
    outedges_.clear();
    outcache_.clear();
    outcflag_.clear();
    incache_.clear();
    incflag_.clear();
}

void TIC::Build2TICFromFile(char *filename) {
#ifdef VERBOSE
    printf("building TIC model...\n");
#endif
    init();

    FILE* f = fopen(filename, "r");
    fscanf(f, "%d %d", &m_, &Z_);
    inedges_.resize(m_);
    outedges_.resize(m_);

    int u, v; double r;
    for (int i = 0; i < m_; i++) {
        fscanf(f, "%d\t%d", &u, &v);
        if (cvt_.count(u) == 0) cvt_[u] = n_++;
        if (cvt_.count(v) == 0) cvt_[v] = n_++;

        outedges_[i].u = inedges_[i].u = cvt_.at(u);
        outedges_[i].v = inedges_[i].v = cvt_.at(v);

        inedges_[i].w.resize(Z_);
        outedges_[i].w.resize(Z_);
        for (int z = 0; z < Z_; z++) {
            fscanf(f, "%lg", &r);
            inedges_[i].w[z] = outedges_[i].w[z] = r;
        }
    }

    qsort_inedges(0, m_-1);  // sort out edges
    qsort_outedges(0, m_-1);  // sort in edges

    RemoveMultiPaths();

    IndexNodesOnEdges();

    incache_.resize(m_);
    outcache_.resize(m_);
    outcflag_.resize(m_, false);
    incflag_.resize(m_, false);

    fclose(f);

    Stats();
}

void TIC::qsort_outedges(int h, int t) {
	if (h<t) 
	{
		int i = h, j = t;
		Edge mid = outedges_[(i+j)/2];
		outedges_[(i+j)/2] = outedges_[i];

		while (i<j) 
		{
			while ((i<j) && ((outedges_[j].u>mid.u)||((outedges_[j].u==mid.u)&&(outedges_[j].v>mid.v))))
				j--;
			if (i<j) {
				outedges_[i] = outedges_[j];
				i++;
			}
			while ((i<j) && ((outedges_[i].u<mid.u)||((outedges_[i].u==mid.u)&&(outedges_[i].v<mid.v))))
				i++;
			if (i<j) {
				outedges_[j] = outedges_[i];
				j--;
			}
		}

		outedges_[i] = mid;
		qsort_outedges(h, i-1);
		qsort_outedges(i+1, t);
	}
}

void TIC::qsort_inedges(int h, int t) {
	if (h<t) 
	{
		int i = h, j = t;
		Edge mid = inedges_[(i+j)/2];
		inedges_[(i+j)/2] = inedges_[i];

		while (i<j) 
		{
			while ((i<j) && ((inedges_[j].v>mid.v)||((inedges_[j].v==mid.v)&&(inedges_[j].u>mid.u))))
				j--;
			if (i<j) {
				inedges_[i] = inedges_[j];
				i++;
			}
			while ((i<j) && ((inedges_[i].v<mid.v)||((inedges_[i].v==mid.v)&&(inedges_[i].u<mid.u))))
				i++;
			if (i<j) {
				inedges_[j] = inedges_[i];
				j--;
			}
		}

		inedges_[i] = mid;
		qsort_inedges(h, i-1);
		qsort_inedges(i+1, t);
	}
}

void TIC::RemoveMultiPaths() {
	int m1 = 0, m2 = 0;
	for (int i=1; i<m_; i++)
	{
		if ((outedges_[i].u != outedges_[m1].u) || (outedges_[i].v != outedges_[m1].v))
		{
			m1++;
			outedges_[m1] = outedges_[i];
		}
	}
	for (int i=1; i<m_; i++)
	{
		if ((inedges_[i].u != inedges_[m2].u) || (inedges_[i].v != inedges_[m2].v))
		{
			m2++;
			inedges_[m2] = inedges_[i];
		}
	}
    assert(m1 == m2);    
	if (m_!=0)
		m_ = m1+1;
}

void TIC::IndexNodesOnEdges() {
    outindex_.resize(n_);   
	for (int i=0; i<m_; i++)
		outindex_[outedges_[i].u] = i;
	for (int i=1; i<n_; i++)
		if (outindex_[i] < outindex_[i-1])
			outindex_[i] = outindex_[i-1];

    inindex_.resize(n_);   
	for (int i=0; i<m_; i++)
		inindex_[inedges_[i].v] = i;
	for (int i=1; i<n_; i++)
		if (inindex_[i] < inindex_[i-1])
			inindex_[i] = inindex_[i-1];
}

void TIC::SaveConverter(char *filename) {
    FILE* f = fopen(filename, "w");
    for (map<int, int>::iterator it = cvt_.begin(); it != cvt_.end(); it++) {
        fprintf(f, "%d\t%d\n", it->first, it->second);
    }
    fclose(f);
}

int TIC::GetOutNeighbor(int vertex) {
	if (vertex == 0)
		return outindex_[vertex]+1;
	else 
		return outindex_[vertex]-outindex_[vertex-1];
}

int TIC::GetInNeighbor(int vertex) {
	if (vertex == 0)
		return inindex_[vertex]+1;
	else 
		return inindex_[vertex]-inindex_[vertex-1];
}

Edge TIC::GetOutEdge(int vertex, int idx)
{
	if (vertex == 0)
		return outedges_[idx];
	else
		return outedges_[outindex_[vertex-1]+1+idx];
}

Edge TIC::GetInEdge(int vertex, int idx)
{
	if (vertex == 0)
		return inedges_[idx];
	else
		return inedges_[inindex_[vertex-1]+1+idx];
}

void TIC::Stats() {
    printf("n: %d\n", n_);
    printf("m: %d\n", m_);
    printf("Z: %d\n", Z_);
    double mw = 0, maxw = 0;
    for	(int i = 0; i < m_; i++) {
        mw = 0;
	for (int z = 0; z < Z_; z++) {
            mw = mw > inedges_[i].w[z] ? mw : inedges_[i].w[z];
	}
        maxw += mw;
    }
    printf("average maxium z weight(affects non-topic spread): %lg\n", maxw/m_);
}

void TIC::Print() {
    printf("==================in========================\n");
    for (unsigned int i = 0; i < inedges_.size(); i++) {
        printf("%d\t%d", inedges_[i].u, inedges_[i].v);
        for (int z = 0; z < Z_; z++) {
            printf("\t%lg", inedges_[i].w[z]);
        }
        printf("\n");
    }
    printf("==================out========================\n");
    for (unsigned int i = 0; i < outedges_.size(); i++) {
        printf("%d\t%d", outedges_[i].u, outedges_[i].v);
        for (int z = 0; z < Z_; z++) {
            printf("\t%lg", outedges_[i].w[z]);
        }
        printf("\n");
    }
}

double TIC::PProb(const vector<double>& v1, const vector<double>& v2) {
    double rst = 0;
    for (unsigned int i = 0; i < v1.size(); i++) {
        rst += v1[i]*v2[i];
    }
    rst = rst > 1 ? 1 : rst;
    return rst;
}

double TIC::MaxPProb(vector<double>& v) {
    double rst = 0;
    for (unsigned int i = 0; i < v.size(); i++) {
        if (rst < v[i]) rst = v[i];
    }
    return rst;
}

bool min_heap_cmp(const di&a, const di&b) {
    return a.first > b.first;
}

double TIC::GetOutEdgeWeight(int vertex, int idx, const vector<double>& item) {
    int index;
	if (vertex == 0)
		index = idx;
	else
		index = outindex_[vertex-1]+1+idx;

    if (outcflag_[index] == false) {
        outcflag_[index] = true;
        Edge e = outedges_[index];
//        printf("vertex %d, e.u %d\n", vertex, e.u);
        outcache_[index] = PProb(e.w, item);
    }
    return outcache_[index];
}

double TIC::GetInEdgeWeight(int vertex, int idx, const vector<double>& item) {
    int index;
	if (vertex == 0)
		index = idx;
	else
		index = inindex_[vertex-1]+1+idx;

    if (incflag_[index] == false) {
        incflag_[index] = true;
        Edge e = inedges_[index];
        incache_[index] = PProb(e.w, item);
    }
    return incache_[index];
}

vector<di> TIC::DijkstraWithStopIds(int u, const double max, const set<int>& stops) {
    vector<di> rst;
    vector<di> seen;  // distance heap
    map<int, double> dist;
    dist[u] = 0;
    seen.push_back(di(0, u));
    make_heap(seen.begin(), seen.end(), min_heap_cmp);

    // variables used in the while loop.
    di top;
    double v_dist, w_dist; int v, w;
    vector<di> neighbors;

    while (!seen.empty()) {
        top = seen.front();
        pop_heap(seen.begin(),seen.end(), min_heap_cmp); seen.pop_back();
        v_dist = top.first; v = top.second;

        if (v_dist > max) {
            break;
        }

        if (stops.count(v) == 1) continue;

        // dist must contain v already
        assert (dist.count(v) != 0);

        if (dist.at(v) < v_dist) {
            // this is an old value selected bofore.
            continue;
        }

        // shortest distance vertex for sure.
        rst.push_back(di(exp(-v_dist), v));

        for (unsigned int i = 0; i < GetOutNeighbor(v); i++) {
            Edge e = GetOutEdge(v, i);
            w_dist = -log(MaxPProb(e.w));
            w = e.v;
            if (dist.count(w) == 0 || w_dist + v_dist < dist.at(w)) {
                dist[w] = w_dist + v_dist;
                seen.push_back(di(w_dist+v_dist, w)); push_heap(seen.begin(), seen.end(), min_heap_cmp);
            }
        }
    }
//    printf("Dijkstra size :%d\n", rst.size());
    return rst;
}

vector<di> TIC::Dijkstra(int u, const double max) {
    set<int> empty;
    return DijkstraWithStopIds(u, max, empty);
}

double TIC::Spread(vector<int> seeds, vector<double> item) {
    srand (time(NULL));
    int n_times = 1000;

    int h, t, k_top = seeds.size();
    vector<int> list;
    set<int> active;
    int resultSize = 0;
    for (int it=0; it < n_times; it++) {
        list.clear();
        active.clear();
        for (int i = 0; i < k_top; i++) 
        {
            list.push_back(seeds[i]);
            active.insert(seeds[i]);
        }
        resultSize += k_top;

        h = 0;
        t = k_top;

        while (h<t) {
            int k = GetOutNeighbor(list[h]);
            for (int i=0; i<k; i++) {
                Edge e = GetOutEdge(list[h], i);
                if (active.count(e.v) == 1) continue;
                if (((double)rand()/(double)RAND_MAX) < PProb(e.w, item)) {
                    list.push_back(e.v);
                    active.insert(e.v);
                    t++;
                    resultSize++;
                }
            }
            h++;
        }
    }
    return (double)resultSize / (double)n_times;
}

void TIC::Convert2PMIA(vector<double> item, char *filename) {
    FILE *f = fopen(filename, "w");
    fprintf(f, "%d %d\n", n_, m_);
    map<pair<int, int>, double> m;
    int u, v;
    for (int i = 0; i < m_; i++) {
        u = outedges_[i].u;
        v = outedges_[i].v;
        m[make_pair(u, v)] = PProb(outedges_[i].w, item);
    }
    for (int i = 0; i < m_; i++) {
        u = outedges_[i].u;
        v = outedges_[i].v;
        pair<int, int>uv = make_pair(u, v);
        pair<int, int>vu = make_pair(v, u);
        fprintf(f, "%d %d %lg %lg\n", u, v,
                m.count(uv) == 0 ? 0 : m[uv],
                m.count(vu) == 0 ? 0 : m[vu]);
        fprintf(f, "%d %d %lg %lg\n", v, u,
                m.count(vu) == 0 ? 0 : m[vu],
                m.count(uv) == 0 ? 0 : m[uv]);
    }
    fclose(f);
}
