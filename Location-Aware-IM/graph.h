#ifndef GRAPH_H
#define GRAPH_H

#include <string.h>
#include "limit.h"
#include <vector>

using namespace std;

struct Edge 
{
	int u,v;	
	double w1,w2;
};

class Graph
{
private:

	static int n;
	static int m;
    static bool edgesSorted;
	static vector<int> degree;
	static vector<int> index;
	static vector<Edge> edges;

	static void qsort_edges(int h, int t);
    static void RemoveMultiPaths();
    static void IndexNodesOnEdges();

public:

	static int	GetN();
	static int	GetM();
	static int	GetDegree(int node);
	static int	GetNeighbor(int node);
	static Edge	GetEdge(int node, int idx);
	static void BuildWC();
	static void Build2WC();
    static void Build2WCFromFile(char* );
};

#endif
