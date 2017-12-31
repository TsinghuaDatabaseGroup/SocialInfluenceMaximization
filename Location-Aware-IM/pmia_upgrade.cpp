#include "pmia_upgrade.h"
#include "graph.h"
#include <stdio.h>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <vector>
#include <cassert>
int SPT_upgrade::n = 0;
int SPT_upgrade::top = 0;
double SPT_upgrade::d[MAX_K];  // incremental influence ?
int SPT_upgrade::list[MAX_K];  // topk
char SPT_upgrade::file[] = "SPT_upgrade_0000.txt";
static char path_log[STR_LEN] = "SPT_upgrade.log.txt";

int SPT_upgrade::k=1;
vector<int> SPT_upgrade::dd(MAX_NODE,0);
double SPT_upgrade::longest = log(100.0);
vector<double> SPT_upgrade::dp(MAX_NODE,1.0);
vector<bool> SPT_upgrade::used(MAX_NODE);
vector<double *>SPT_upgrade::self(MAX_NODE);
vector<int> SPT_upgrade::lastupdate(MAX_NODE,-1);
vector<double *>SPT_upgrade::delta(MAX_NODE);
vector<int *>SPT_upgrade::children(MAX_NODE,NULL);
vector<int *>SPT_upgrade::path(MAX_NODE,NULL);

int *SPT_upgrade::S;
double *SPT_upgrade::distance, *SPT_upgrade::b;
int *SPT_upgrade::numchild;
int *SPT_upgrade::queue;
int *SPT_upgrade::heap;
int *SPT_upgrade::childlist, *SPT_upgrade::oldchildlist, *SPT_upgrade::parent;
bool *SPT_upgrade::validlist[MAX_K]={NULL};

/* my own enhancement to support subgraph PMIA
 */
int SPT_upgrade::bbound=-1;
vector<int *>SPT_upgrade::bchildren(MAX_NODE,NULL);
vector<int *>SPT_upgrade::bpath(MAX_NODE,NULL);
vector<double *>SPT_upgrade::bself(MAX_NODE, NULL);
vector<int> SPT_upgrade::bdd(MAX_NODE,0);
vector<double *>SPT_upgrade::bdelta(MAX_NODE, NULL);

vector<bool> SPT_upgrade::is_native(MAX_NODE);


int SPT_upgrade::GetMax(int round)  // get node with maximum incremental influence. round is the number of selected seeds.
{
    double max = -1000000.0;
    int mp = -1;
    for (int j=0; j<n; j++)
        if (!used[j] && lastupdate[j]!=round)  // j is not selected as seed yet. And it has been updated in current round.
        {
            double tmp = dp[j];
            if (tmp >max)
            {
                max = tmp;
                mp = j;
            }
        }
    return mp;
}

int SPT_upgrade::generateSPT_newfrom(int round, int node){
    int top=0, bottom=0;
    distance[node]=0;
    heap[0]=node;
    top++;
    while (true){
        //pop out of heap
        if (distance[heap[0]]<longest) S[heap[0]]=-1; 
        else break;
        childlist[bottom++]=heap[0];
        for (int i=0;i<Graph::GetNeighbor(heap[0]);i++){
            Edge e=Graph::GetEdge(heap[0],i);
            if (used[e.v] || S[e.v]<0) continue;
            if (distance[e.v]>distance[heap[0]]+e.w1+EPS) {
                if (S[e.v]>=n){
                    distance[e.v]=distance[heap[0]]+e.w1;
                    heap[top]=e.v;
                    int j=top++, x=(j-1)/2;
                    double temp=distance[heap[j]];
                    while (j>0) {
                        if (distance[heap[x]]>temp){
                            heap[j]=heap[x];
                            if (S[heap[j]]<n) S[heap[j]]=j;
                            j=x;
                            x=(j-1)/2;
                        }
                        else break;
                    }
                    heap[j]=e.v;
                    S[e.v]=j;
                }
                else{  // e.v has been added into heap before?
                    distance[e.v]=distance[heap[0]]+e.w1;
                    int j=S[e.v], x=(j-1)/2;
                    double temp=distance[heap[j]];
                    while (j>0) {
                        if (distance[heap[x]]>temp){
                            heap[j]=heap[x];
                            if (S[heap[j]]<n) S[heap[j]]=j;
                            j=x;
                            x=(j-1)/2;
                        }
                        else break;
                    }
                    heap[j]=e.v;
                    S[e.v]=j;
                }
            }//endif

        }// end for
        heap[0]=heap[--top];
        if (!top) break;
        //siftdown
        int j=0, x=j*2+1;
        double temp=distance[heap[j]];
        while (x<top){
            if (x+1<top && distance[heap[x+1]]<distance[heap[x]]) x++;
            if (distance[heap[x]]<temp){
                heap[j]=heap[x];
                S[heap[j]]=j;
                j=x; x=j*2+1;
            }
            else break;
        }
        heap[j]=heap[top];
        if (S[heap[j]]<n) S[heap[j]]=j;
    }  // end while
    for (int i=0;i<bottom;i++)  // restore
    {
        int child=childlist[i];
        distance[child]=longest;
        S[child]=n;
        parent[child]=-1;
    }
    //update tree node set and heuristic
    validlist[round]=new bool[n];
    memset(validlist[round],false,n);
    for (int i=0;i<bottom;i++){
        int child=childlist[i];  // child is members in out-arborescence of node
        oldchildlist[i]=child;
        if (dd[child]>0 && is_native[child]) {  // only care about natives
            // subtract previous incremental influence
            for (int j=0;j<dd[child];j++) {
                dp[children[child][j]]-=delta[child][j]*(1-self[child][j]);
            }

            //block some nodes in child's inverse tree
            for (int j=0;j<round;j++) {
                if (validlist[j][child]) {
                    bool invalid=false;
                    int k;
                    for (k=0;k<dd[child] && children[child][k]!=list[j];k++) ;
                    if (k==dd[child]) {
                        continue;
                    }

                    for (;k > 0 && children[child][k]!=node;k=path[child][k]) {
                    }  // jth seed is not blocked by node.
                    validlist[j][child]=k==0;  //k==0 means ok
                }
            }
            validlist[round][child]=true;
        }
    }
    // 	fprintf(f_log, "PMIOA source: %d\n", node);
    //     for (int i=1; i< bottom; i++) fprintf(f_log, "%d\n", childlist[i]);
    //     fprintf(f_log, "\n");
    //    for (int i=0;i<bottom;i++) 			generateSPT_newto(oldchildlist[i]);
    for (int i=0;i<bottom;i++) {
        if (is_native[oldchildlist[i]])  // only care about natives
            generateSPT_newto(oldchildlist[i]);
    }

//    for (int i=0;i<bottom;i++) {
//        printf("node %d, ap[%d] == %f\n", node, oldchildlist[i], self[oldchildlist[i]][0]);
//    }

    return bottom;
}

int SPT_upgrade::generateSPT_newto(int node){
    int top=0, bottom=0;  // bottom is the number members in arborescence
    if (used[node]) {
        // node is already selected, set its arborescence to contain only
        // itself.
        dd[node]=1;  
        path[node][0]=0;
        self[node][0]=1;
        delta[node][0]=0;
        return 1;
    }
    distance[node]=0;
    heap[0]=node;
    top++;
    parent[node]=node;
    b[node]=1;
    while (true){
        //stack out of heap
        if (distance[heap[0]]<longest) S[heap[0]]=-bottom-1; 
        else break;
        childlist[bottom++]=heap[0];
        if (parent[heap[0]]!=heap[0]) numchild[parent[heap[0]]]++;
        if (!used[heap[0]])
            for (int i=0;i<Graph::GetNeighbor(heap[0]);i++){
                Edge e=Graph::GetEdge(heap[0],i);
                // e.v is invalid seed in node's arborescence
                if (used[e.v] && !validlist[lastupdate[e.v]][node] || S[e.v]<0) continue;  
                if (distance[e.v]>distance[heap[0]]+e.w2+EPS) {
                    parent[e.v]=heap[0];
                    b[e.v]=exp(-e.w2);
                    if (S[e.v]>=n){
                        distance[e.v]=distance[heap[0]]+e.w2;
                        heap[top]=e.v;
                        int j=top++, x=(j-1)/2;
                        double temp=distance[heap[j]];
                        while (j>0) {
                            if (distance[heap[x]]>temp){
                                heap[j]=heap[x];
                                if (S[heap[j]]<n) S[heap[j]]=j;
                                j=x;
                                x=(j-1)/2;
                            }
                            else break;
                        }
                        heap[j]=e.v;
                        S[heap[j]]=j;
                    }
                    else{
                        distance[e.v]=distance[heap[0]]+e.w2;
                        int j=S[e.v], x=(j-1)/2;
                        double temp=distance[heap[j]];
                        while (j>0) {
                            if (distance[heap[x]]>temp){
                                heap[j]=heap[x];
                                if (S[heap[j]]<n) S[heap[j]]=j;
                                j=x;
                                x=(j-1)/2;
                            }
                            else break;
                        }
                        heap[j]=e.v;
                        S[e.v]=j;
                    }
                }//endif

            }// end for
        heap[0]=heap[--top];
        if (!top) break;
        //siftdown
        int j=0, x=j*2+1;
        double temp=distance[heap[j]];
        while (x<top){
            if (x+1<top && distance[heap[x+1]]<distance[heap[x]]) x++;
            if (distance[heap[x]]<temp){
                heap[j]=heap[x];
                S[heap[j]]=j;
                j=x; x=j*2+1;
            }
            else break;
        }
        heap[j]=heap[top];
        if (S[heap[j]]<n) S[heap[j]]=j;
    }  // end while
    //update tree node set and heuristic
    if (!dd[node]){
        children[node]=new int[bottom];
        delta[node]=new double[bottom];
        self[node]=new double[bottom];
        path[node]=new int[bottom];
    }
    dd[node]=bottom; 
    int head=0, tail=0;
    for (int i=0;i<bottom;i++){
        children[node][i]=childlist[i];
        if (numchild[childlist[i]])	self[node][i]=1;  // childlist[i] is a leaf
        else {
            self[node][i]=used[childlist[i]]?1:0;
            queue[tail++]=i;  // queue is a list of leaf nodes in arborescence.
        }
        path[node][i]=-S[parent[childlist[i]]]-1;  // path is ith child's parent's idx
    }
    for (int i=0;i<bottom;i++)  // restore to default value
    {
        int child=childlist[i];
        distance[child]=longest;
        S[child]=n;
        parent[child]=-1;
    }

    int x,u;
    while (head<tail) {  // calculate activation probability
        x=queue[head++];
        u=path[node][x];
        //b[v] is the weight of edge child -> parent in an arborescence
        if (numchild[childlist[u]] > 0) {
            self[node][u]*=(1-self[node][x]*b[childlist[x]]);
        }
        if (!--numchild[childlist[u]]) {
            self[node][u]=1-self[node][u];
            queue[tail++]=u;
        }		
    }
    numchild[node]=0;  // restore
    delta[node][queue[--head]]=1;
//    printf("node %d, ap: %f\n", childlist[i], self[node][i]);
    dp[node]+=1-self[node][x];
    for (head--;head>=0;head--) {
        x=queue[head], u=path[node][x];
        delta[node][x]=(1-self[node][u])/(1-self[node][x]*b[childlist[x]])*b[childlist[x]]*delta[node][u];
        dp[childlist[x]]+=delta[node][x]*(1-self[node][x]);
    }
//    fprintf(f_log, "PMIIA source: %d\n", node);
//      for (int i=0; i< bottom; i++) fprintf(f_log, "%d ap:%f alpha:%f, pp:%f\n", childlist[i], self[node][i], delta[node][i], b[childlist[i]]);
//      fprintf(f_log, "\n");

//    printf("PMIIA source: %d\n", node);
//    for (int i=0; i< dd[node]; i++) {
//          if (path[node][i] != -1)
//            printf("%d ap:%f alpha:%f, pp:%f, parent:%d\n", childlist[i],
//                self[node][i], delta[node][i], b[childlist[i]],
//                children[node][path[node][i]]);
//          else
//            printf("%d ap:%f alpha:%f, pp:%f, parent:%d\n", childlist[i],
//                self[node][i], delta[node][i], b[childlist[i]],
//                path[node][i]);
//    }
//    for (int i=0;i<bottom;i++)
//        printf("%d ", childlist[i]);
//    printf("\n");
//    for (int i=0;i<bottom;i++) 
//        printf("%d ", path[node][i]);
//    printf("\n");
//    for (int i=0;i<bottom;i++) 
//        printf("%f ", b[childlist[i]]);
//    printf("\n");
    return bottom;
}

// double SPT_upgrade::Build(int num, int bound)  // num is topk
// {
//     clock_t start, end;
//     start = clock();
//     n = Graph::GetN();
//     longest=log(double(bound));  // influence threshold
//     top = num;
//     double treesize=0;
//     S = new int[n];
//     distance = new double[n];
//     b = new double[n];
//     heap = new int[n];
//     childlist = new int[n];
//     oldchildlist = new int[n];
//     parent = new int[n];
//     numchild = new int[n];
//     queue = new int[n];
// 
//     used.resize(n);
//     lastupdate.resize(n);
//     children.resize(n);
//     dp.resize(n);
//     self.resize(n);
//     dd.resize(n);
//     delta.resize(n);
//     path.resize(n);
//     // no use in this function
// //    int set[SET_SIZE];  // a temporal placeholder for list, actually no use here
// 
//     double old = 0.0;  // no use in this function
// 
//     int i=0;
//     for (i=0; i<n; i++)
//     {
//         lastupdate[i] = -1;
//         children[i]=NULL;
//         dp[i]=0.0;
//         self[i]=NULL;
//         used[i]=false;
//     }
//     for (i=0; i<n; i++)
//         dd[i] = 0;
//     for (int i=0;i<n;i++) distance[i]=longest;
//     for (int i=0;i<n;i++) S[i]=n;
//     for (int i=0;i<n;i++) parent[i]=-1;
//     for (int i=0;i<n;i++) numchild[i]=0;
// 
//     for (i=0;i<n;i++)
//     {
//         double size=generateSPT_newto(i);
//         treesize+=size;
//     }
//     end = clock();
//     double	timer = (double)(end - start);
// //    printf("initialization time: %lg \n", timer / CLOCKS_PER_SEC);
// 
// //    for (int _ = 0;  _ < n; _++)
// //        printf("%d's inc_inf: %f\n", _, dp[_]);
// 
//     i=0;
//     double max = -1000000.0;
//     int mp;
//     {
//         int x=GetMax(i);
// //        set[i] = x;  // no use in this function
// 
//         lastupdate[x] = i;			
//         double improve=dp[x];
//         if (improve > max) {
//             max=improve;
//             mp=x;
//         }
//     }
//     used[mp] = true;
// //    set[i] = mp;  // no use in this function
// 
//     list[i] = mp;
//     d[i] = max;
//     old+=d[i];  // no use in this function
// //    printf("selected: %d inc_inf:%f\n", mp, max);
//     generateSPT_newfrom(i, mp);
// //    for (int _ = 0;  _ < n; _++)
// //        printf("%d's inc_inf: %f\n", _, dp[_]);
// 
//     for (i=1; i<top; i++)
//     {
//         max = -1000000.0;
//         int x=GetMax(i);
// //        set[i] = x;  // no use in this function
// 
//         lastupdate[x] = i;			
//         double improve=dp[x];
//         if (improve > max) {
//             max=improve;
//             mp=x;
//         }
//         used[mp] = true;
//  //       set[i] = mp;  // no use in this function
// 
//         list[i] = mp;
//         d[i] = max;
//         old+=d[i];  // no use in this function
// //        printf("selected: %d inc_inf:%f\n", mp, max);
//         generateSPT_newfrom(i, mp);
// //        for (int _ = 0;  _ < n; _++)
// //            printf("%d's inc_inf: %f\n", _, dp[_]);
//     }
//     int ct=0;
//     delete[] childlist;
//     delete[] oldchildlist;
//     delete[] distance;
//     delete[] S;
//     delete[] heap;
//     delete[] b;
//     delete[] parent;
//     delete[] numchild;
//     delete[] queue;
// 
//     for (i=0;i<n;i++)
//         if (dd[i]) {
//             delete[] children[i];
//             delete[] delta[i];
//             delete[] self[i];
//             delete[] path[i];
//             dd[i]=0;
//         }
//     for (i=0;i<top;i++)
//         delete[] validlist[i];
// 
// //    sprintf(file,"SPT_upgrade_%04d.txt", bound);
// //    FILE *out = fopen(file, "w");
// //    fprintf(out, "%d\n", top);
// //    for (i=0; i<top; i++)
// //        fprintf(out, "%d\t%g\n", list[i], d[i]);
// //    fclose(out);
// //    printf("%lg ",treesize/n);
//     return treesize/n;
// }

double SPT_upgrade::init(int bound)
{
    clock_t start, end;
    start = clock();
    n = Graph::GetN();
    longest=log(double(bound));  // influence threshold
    double treesize=0;
    S = new int[n];
    distance = new double[n];
    b = new double[n];
    heap = new int[n];
    childlist = new int[n];
    oldchildlist = new int[n];
    parent = new int[n];
    numchild = new int[n];
    queue = new int[n];

    used.resize(n);
    lastupdate.resize(n);
    children.resize(n);
    dp.resize(n);
    self.resize(n);
    dd.resize(n);
    delta.resize(n);
    path.resize(n);

    // resize my enhancement
    bchildren.resize(n);
    bpath.resize(n);
    bself.resize(n);
    bdd.resize(n);
    bdelta.resize(n);

    int i;
    for (i=0; i<n; i++)
    {
        lastupdate[i] = -1;
        children[i]=NULL;
        dp[i]=0.0;
        self[i]=NULL;
        used[i]=false;
    }
    for (i=0; i<n; i++)
        dd[i] = 0;
    for (int i=0;i<n;i++) distance[i]=longest;
    for (int i=0;i<n;i++) S[i]=n;
    for (int i=0;i<n;i++) parent[i]=-1;
    for (int i=0;i<n;i++) numchild[i]=0;

    for (i=0; i<n; i++)
    {
        bchildren[i]=NULL;
        bself[i]=NULL;
        bdelta[i] = NULL;
        bpath[i] = NULL;
        bdd[i] = 0;
    }
    for (i=0;i<n;i++) {
        double size=generateSPT_newto(i);
        treesize+=size;
    }
    for (int i=0;i<n;i++) {
        if (dd[i]) {
            bchildren[i] = new int[dd[i]];
            bdelta[i] = new double[dd[i]];
            bself[i] = new double[dd[i]];
            bpath[i] = new int[dd[i]];
            bdd[i]=dd[i];

            for (int j = 0; j < dd[i]; j++) {
                bchildren[i][j] = children[i][j];
                bdelta[i][j] = delta[i][j];
                bself[i][j] = self[i][j];
                bpath[i][j] = path[i][j];
            }
        }
    }
    end = clock();
    double	timer = (double)(end - start);
    printf("initialization time: %lg \n", timer / CLOCKS_PER_SEC);
    return treesize/n;
}

void SPT_upgrade::reset() {
    clock_t start, end;
    start = clock();
    int i=0;
    for (i=0; i<n; i++)
    {
        lastupdate[i] = -1;
        dp[i]=0.0;
        used[i]=false;
    }
//    for (int i=0;i<n;i++) distance[i]=longest;
//    for (int i=0;i<n;i++) S[i]=n;
//    for (int i=0;i<n;i++) parent[i]=-1;
//    for (int i=0;i<n;i++) numchild[i]=0;

    for (int i=0;i<n;i++) {
        if (bdd[i]) {
            dd[i]=bdd[i];
            for (int j = 0; j < bdd[i]; j++) {
                children[i][j] = bchildren[i][j];
                delta[i][j] = bdelta[i][j];
                self[i][j] = bself[i][j];
                path[i][j] = bpath[i][j];
            }
        }
    }
    for (int i=0;i<top;i++)
        delete[] validlist[i];
    end = clock();
    double	timer = (double)(end - start);
//        printf("copy initialization data time: %lg \n", timer / CLOCKS_PER_SEC);
}

void SPT_upgrade::select(int num)
{
    double timerPMIA = 0.0;
    clock_t start;
    for (int i=0;i<n;i++) {
        if (dd[i] > 0 && is_native[i]) {
            for (int j = 0; j < dd[i]; j++) {
                dp[children[i][j]] += delta[i][j]*(1 - self[i][j]);
            }
        }
    }
//    if (PROTOTYPE) {
//        printf("init dp:\n");
//        for (int i = 0; i < n; i++) {
//            printf("dp[%d] == %f\n", i, dp[i]);
//        }
//    }
    printf("\n");
    top = num;
    int i  = 0;
    double max = -1000000.0;
    int mp;
    {
        int x=GetMax(i);
        lastupdate[x] = i;			
        double improve=dp[x];
        if (improve > max) {
            max=improve;
            mp=x;
        }
    }
    used[mp] = true;
    list[i] = mp;
    d[i] = max;
//    printf("%dth selected: %d inc_inf:%f\n", i+1, mp, max);
    start = clock();
    generateSPT_newfrom(i, mp);
    timerPMIA += (double)(clock() - start) / CLOCKS_PER_SEC;
//    if (PROTOTYPE) {
//        for (int i = 0; i < n; i++) {
//            printf("inc_infl[%d] == %f\n", i, dp[i]);
//        }
//        printf("\n");
//    }
    for (i=1; i<top; i++)
    {
        max = -1000000.0;
        int x=GetMax(i);

        lastupdate[x] = i;			
        double improve=dp[x];
        if (improve > max) {
            max=improve;
            mp=x;
        }
        used[mp] = true;

        list[i] = mp;
        d[i] = max;
//        printf("%dth selected: %d inc_inf:%f\n", i+1, mp, max);
        start = clock();
        generateSPT_newfrom(i, mp);
        timerPMIA += (double)(clock() - start) / CLOCKS_PER_SEC;
//        if (PROTOTYPE) {
//            for (int _ = 0;  _ < n; _++)
//                printf("%d's inc_inf: %f\n", _, dp[_]);
//            printf("\n");
//        }
    }
    printf("PMIA time: %f\n", timerPMIA);
}

void SPT_upgrade::gc()
{
    delete[] childlist;
    delete[] oldchildlist;
    delete[] distance;
    delete[] S;
    delete[] heap;
    delete[] b;
    delete[] parent;
    delete[] numchild;
    delete[] queue;

    for (int i=0;i<n;i++)
        if (dd[i]) {
            delete[] children[i];
            delete[] delta[i];
            delete[] self[i];
            delete[] path[i];
            dd[i]=0;
        }
    for (int i=0;i<top;i++)
        delete[] validlist[i];

    for (int i=0;i<n;i++) {
        if (bdd[i]) {
            delete[] bchildren[i];
            delete[] bdelta[i];
            delete[] bself[i];
            delete[] bpath[i];
            bdd[i]=0;
        }
    }
}

void SPT_upgrade::exit()
{
    gc();
}

void SPT_upgrade::reset_native()
{
    int n = Graph::GetN();
    is_native.resize(n);
    for (int _ = 0; _ < n; _++) {
        is_native[_] = false;
    }
}

void SPT_upgrade::BuildFromFile(int bound)
{
    n = Graph::GetN();
    sprintf(file,"SPT_upgrade_%04d.txt", bound);
    FILE* in = fopen(file, "r");
    fscanf(in, "%ld", &top);
    for (int i=0; i<top; i++)
        fscanf(in, "%ld %Lg", &list[i], &d[i]);
    fclose(in);
}

int SPT_upgrade::GetNode(int i)
{
    if (i<0)
        return -1;
    if (i>=top) 
        return -1;
    return list[i];
}


char* SPT_upgrade::filename(int bound)
{
    sprintf(file,"SPT_upgrade_%04d.txt", bound);
    return file;
}
