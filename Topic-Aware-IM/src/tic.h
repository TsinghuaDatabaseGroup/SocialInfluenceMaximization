#ifndef TIC_H
#define TIC_H

#include "global.h"
#include <vector>
#include <map>

using namespace std;

typedef std::pair<double, int> di;

struct Edge 
{
	int u,v;
	vector<double> w;
};


class TIC
{
private:
	static int n_;
	static int m_;
	static int Z_;

    static map<int, int> cvt_;  // ID converter

//	static vector<int> indegree_;
//	static vector<int> outdegree_;
	static vector<int> inindex_;
	static vector<int> outindex_;
	static vector<Edge> inedges_;
	static vector<Edge> outedges_;
	static vector<double> outcache_;
	static vector<double> incache_;
	static vector<bool> outcflag_;  // is cached ?
	static vector<bool> incflag_;

    static void init();

	static void qsort_inedges(int h, int t);
	static void qsort_outedges(int h, int t);
    static void RemoveMultiPaths();
    static void IndexNodesOnEdges();

public:
	static int	GetN();
	static int	GetM();
	static int	GetZ();

//	static int	GetOutDegree(int vertex);
//	static int	GetInDegree(int vertex);
	static int	GetOutNeighbor(int vertex);
	static int	GetInNeighbor(int vertex);
	static Edge	GetOutEdge(int vertex, int idx);
	static Edge	GetInEdge(int vertex, int idx);
	static double	GetOutEdgeWeight(int vertex, int idx, const vector<double>&);
	static double	GetInEdgeWeight(int vertex, int idx, const vector<double>&);

	static void Build2TICFromFile(char *);
    static void Print();
    static void SaveConverter(char *);
    static void Stats();

    static double MaxPProb(vector<double>&);
    static double PProb(const vector<double>&, const vector<double>&);
    static vector<di> Dijkstra(int u, const double theta);
    static vector<di> DijkstraWithStopIds(int u, const double theta, const set<int>&);

    static double Spread(vector<int>, vector<double>);
    
    static void Convert2PMIA(vector<double>, char *);
    static void resetCache();
};

#endif
