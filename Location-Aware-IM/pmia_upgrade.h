#ifndef SPT_upgrade_H
#define SPT_upgrade_H

#include "limit.h"
#include <vector>
using namespace std;

class SPT_upgrade
{
private:
	static int n;
	static int top;
	static int k;
	static double d[MAX_K];  // incremental influence of topk nodes
	static int list[MAX_K];
	static char file[STR_LEN];
	static int topk[MAX_K];
	static vector<int> dd;  // number of members in in-arborescence
	static double longest;
	static vector<double> dp;  // incremental influence
	static vector<bool> used;  // selected as topk
	static vector<double*> self;  // self[node][i] is node's ith child's activation probability in node's in-arborescence
    // r = lastupdate[i] means node i get selected as seed at r step (total k steps)
	static vector<int> lastupdate;
	static vector<double *>delta;  // alpha in the paper, delta[node][i] is node's ith child i's alpha over node
    // children consists of PMIIA members, children[node][i] = childlist[i] in in-arborescence of node
	static vector<int *>children, path;  
//    static vector<bool *>iv; // invalid children (blocked, no longer connected)
    // S[j] is j's rank in childlist + 1, then take the negative. S[j] < 0 means j is a member of in/out-arborescence.
	static int *S, *numchild, *queue;
	static double *distance, *b;
	static int *heap;
    // childlist is all the members of in/out-arborescence, ordered by their distance.
    // parent[i] is id of #i node's successor.
	static int *childlist, *oldchildlist, *parent;
	static bool *validlist[MAX_K];  // validlist[round][node] means round_th selected seed is invalid in node's in-arborescence

    /* my own enhancement to support subgraph PMIA
     */
    static int bbound;
	static vector<int *>bchildren, bpath;  //  b means "Back up"
	static vector<double*> bself;  // self[node][i] is node's ith child's activation probability in node's in-arborescence
	static vector<int> bdd;  // number of members in in-arborescence
	static vector<double *>bdelta;  // alpha in the paper, delta[node][i] is node's ith child i's alpha over node

public:
	static double Build(int num, int bound);
//	static double Build(int num, int k, int bound, double (*Run)(int num_iter, int size, int set[]), double (*RunFast)(int num_iter, int size, int set[]));
	static void BuildFromFile(int bound);
	static int GetNode(int i);
	static int GetMax(int round);
	static int GetMax0(int round);
	static int generateSPT_newfrom(int round, int node);
	static void updateSPT_to(int node);
	static int generateSPT_newto(int node);
//	static int generateSPT_newto0(int node);
//	static int count(int node);
	static char* filename(int bound);
    
    /* my own enhancement to support subgraph PMIA
     */
    static vector<bool> is_native;  // if this node is a native user.
    static int init_bound;  // bound used to initialize PMIA globally

    // initialization of PMIA algorithm globally
    static void reset_native();
    static void reset();
    static double init(int bound);
    static void gc();
    static void exit();
    // pmia selection main loop
    static void select(int num);
};

#endif

