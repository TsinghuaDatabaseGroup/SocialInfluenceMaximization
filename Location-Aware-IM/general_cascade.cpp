#include "graph.h"
#include "general_cascade.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "limit.h"
//#include <time.h>

int GeneralCascade::n = 0;
int GeneralCascade::m = 0;
int	GeneralCascade::targetSize = 0;
int	GeneralCascade::resultSize = 0;
int	GeneralCascade::target[MAX_K]={0};
bool GeneralCascade::built = false;

void GeneralCascade::Build()
{
	n = Graph::GetN();
	m = Graph::GetM();
	built = true;
}

void GeneralCascade::SetTarget(int size, int set[])
{
	if (!built) Build();
	targetSize = size;
	for (int i=0; i<size; i++)
		target[i] = set[i];
	resultSize = 0;
}

double GeneralCascade::Run(int num_iter, int size, int set[])
{
	SetTarget(size, set);

	int		h, t;
	int		*list  = new int[n];
	bool	*active= new bool[n];

	for (int it=0; it<num_iter; it++)
	{
		//printf("%d\n",it);
		memset(active, 0, sizeof(bool)*n);
		for (int i=0; i<targetSize; i++) 
		{
			list[i] = target[i];
			active[target[i]] = true;
            if (Args::isNative(target[i]))
                resultSize += 1;
		}
//		resultSize += targetSize;

		h = 0;
		t = targetSize;

		while (h<t) 
		{
			int k = Graph::GetNeighbor(list[h]);
			//printf("%d \n",k);
			for (int i=0; i<k; i++)
			{
				Edge e = Graph::GetEdge(list[h], i);
				if (active[e.v]) continue;
				//printf("%d %d %lg %lg\n", e.u, e.v, exp(-e.w1), exp(-e.w2));
		//for (int j=0; j< e.c; j++)
					if (((double)rand()/(double)RAND_MAX) < exp(-e.w1))
					{
						list[t] = e.v;
						active[e.v] = true;
						t++;
                        if (Args::isNative(e.v))
						    resultSize++;
						//break;
					}
			}
			h++;
		}
	}
	delete[] active;
	delete[] list;
	return (double)resultSize / (double)num_iter;
}
