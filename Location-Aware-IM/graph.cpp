#include "graph.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>

using namespace std;

int  Graph::n = 0;
int  Graph::m = 0;
bool Graph::edgesSorted = false;
vector<int> Graph::index;
vector<int> Graph::degree;
vector<Edge> Graph::edges;

void Graph::qsort_edges(int h, int t)
{
    edgesSorted = true;
	if (h<t) 
	{
		int i = h, j = t;
		Edge mid = edges[(i+j)/2];
		edges[(i+j)/2] = edges[i];

		while (i<j) 
		{
			while ((i<j) && ((edges[j].u>mid.u)||((edges[j].u==mid.u)&&(edges[j].v>mid.v))))
				j--;
			if (i<j) {
				edges[i] = edges[j];
				i++;
			}
			while ((i<j) && ((edges[i].u<mid.u)||((edges[i].u==mid.u)&&(edges[i].v<mid.v))))
				i++;
			if (i<j) {
				edges[j] = edges[i];
				j--;
			}
		}

		edges[i] = mid;
		qsort_edges(h, i-1);
		qsort_edges(i+1, t);
	}
}

int Graph::GetN()
{
	return n;
}

int Graph::GetM()
{
	return m;
}

int Graph::GetDegree(int node)
{
	return degree[node];
}

int Graph::GetNeighbor(int node)
{
	if (node == 0)
		return index[node]+1;
	else 
		return index[node]-index[node-1];
}

Edge Graph::GetEdge(int node, int idx)
{
	if (node == 0)
		return edges[idx];
	else
		return edges[index[node-1]+1+idx];
}

void Graph::RemoveMultiPaths() {
    // not allow multi-paths between any two vertices 
    //
    // must sort edges first
    if (!edgesSorted)
	    qsort_edges(0, 2*m-1);

	int m1 = 0;
	for (int i=1; i<2*m; i++)
	{
		if ((edges[i].u != edges[m1].u) || (edges[i].v != edges[m1].v))
		{
			m1++;
			edges[m1] = edges[i];
		}
	}
	if (m!=0)
		m = m1+1;
//    for (int i = 0; i < m; i++) {
//        cout << edges[i].u << "-->" << edges[i].v << ": " << exp(-edges[i].w1) << endl;
//    }
}

void Graph::IndexNodesOnEdges() {
    // index nodes on edges vector
    // edges with same source nodes are aggregated together
    //
    // must sort edges first
    if (!edgesSorted)
	    qsort_edges(0, 2*m-1);

	index.resize(n, 0);
	for (int i=0; i<m; i++)
		index[edges[i].u] = i;
	for (int i=1; i<n; i++)
		if (index[i] < index[i-1])
			index[i] = index[i-1];
}


void Graph::Build2WCFromFile(char* filename)
{
    FILE* f = fopen(filename, "r");

	fscanf(f, "%d %d", &n, &m);	
	degree.resize(n, 0);
	edges.resize(2*m);

	for (int i=0; i<2*m; i++)
	{
		fscanf(f, "%d %d %lg %lg", &edges[i].u, &edges[i].v, &edges[i].w1, &edges[i].w2);
		edges[i].u--;
		edges[i].v--;
		edges[i].w1=-log(edges[i].w1);
		edges[i].w2=-log(edges[i].w2);
		degree[edges[i].u]++;
	}

	qsort_edges(0, 2*m-1);

    RemoveMultiPaths();

    IndexNodesOnEdges();

    fclose(f);
}

void Graph::Build2WC()
{
	scanf("%d %d", &n, &m);	
	degree.resize(n, 0);
	edges.resize(2*m);

	for (int i=0; i<2*m; i++)
	{
		scanf("%d %d %lg %lg", &edges[i].u, &edges[i].v, &edges[i].w1, &edges[i].w2);
		edges[i].u--;
		edges[i].v--;
		edges[i].w1=-log(edges[i].w1);
		edges[i].w2=-log(edges[i].w2);
		degree[edges[i].u]++;
	}

	qsort_edges(0, 2*m-1);

    RemoveMultiPaths();

    IndexNodesOnEdges();
}

// void Graph::Stats()
// {
// 	printf("number of vertices:\t%d\n",n);
// 	printf("number of edges:\t%d\n",m/2);
// 	printf("density:\t%lg\n",double(m)/n/(n-1));
// 	int maxdegree=0;
// 	double tdegree=0.0;
// 	int i,j,k;
// 	for (i=0;i<n;i++)
// 	{
// 		if (degree[i]>maxdegree) maxdegree=degree[i];
// 		tdegree+=degree[i];
// 		//if (degree[i]%2) printf("%d\n", i);
// 	}
// 	printf("average degree:\t%lg\n",tdegree/n);
// 	printf("maximal degree:\t%d\n",maxdegree);
// 	int maxcmp=0,ncmp=0;
// 	bool *used=new bool[n];
// 	memset(used,0,n);
// 	while (1)
// 	{
// 		queue<int> q;
// 		for (i=0;i<n;i++)
// 			if (!used[i]) break;
// 		if (i==n) break;
// 		ncmp++;
// 		int cmpsize=0;
// 		q.push(i);
// 		used[i]=true;
// 		while (!q.empty())
// 		{
// 			k=q.front();
// 			q.pop();
// 			cmpsize++;
// 			j=GetNeighbor(k);
// 			for (i=0;i<j;i++)
// 			{
// 				Edge e=GetEdge(k,i);
// 				if (used[e.v]) continue;
// 				q.push(e.v);
// 				used[e.v]=true;
// 			}
// 		}
// 		if (cmpsize>maxcmp) maxcmp=cmpsize;
// 	}
// 	printf("# of connected component:\t%d\n",ncmp);
// 	printf("largest component size:\t%d\n",maxcmp);
// 	printf("average component size:\t%lg\n",double(n)/ncmp);
// }
