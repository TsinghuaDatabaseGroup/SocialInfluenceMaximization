#include<cmath>
#include<stdio.h>
#include<map>
#include<set>

#include "mip.h"
#include "Cache.h"

vector<double *> Mip::mip;
vector<int *> Mip::inflees;
vector<int> Mip::numInflees;
vector<bool> Mip::changed;

vector<double *> Mip::bmip;
vector<int *> Mip::bInflees;
vector<int> Mip::bNumInflees;

vector<int *> Mip::outPath;
vector<double *> Mip::outAp;
vector<double *> Mip::outB;

vector<vector<int> > Mip::cdds;
vector<int> Mip::numCdds;

// hints
// vector<double> Mip::hap;
// vector<double> Mip::uap;
// double Mip::lowerHint;

Heap *Mip::maxHeap;  // \Heap in paper
Heap *Mip::lMaxHeap;  // \Heap in paper
int Mip::n;
int Mip::count;
int Mip::estCount;
double Mip::incInflTimer;

vector<bool> Mip::used;  // selected as topk, global structure

vector<double> Mip::ap;
vector<double> Mip::sap;
vector<vector<int> > Mip::path;
vector<vector<int> > Mip::children;
vector<vector<double> > Mip::bb;
vector<double> Mip::b;
vector<int> Mip::numInfrs;  // total number of influencers
double Mip::longest;

// dijkstar
vector<int> Mip::S; 
vector<int> Mip::numchild; 
vector<int> Mip::parent;
vector<int> Mip::heap;
vector<int> Mip::queue;
vector<int> Mip::childlist;
vector<double> Mip::distance;

vector<int> Mip::status;
vector<bool> Mip::is_native;

void Mip::init(double _longest) {
    cout << "initialization MIPs..." << endl;
    // calculate \prob(u \leadsto \v) by dijkstar shortestPathsFrom u
    n = Graph::GetN();
    count = 0;
    estCount = 0;
    incInflTimer = 0;

    // can be restored later
    mip.resize(n, NULL);  
    inflees.resize(n, NULL);
    numInflees.resize(n, 0);

    changed.resize(n, false);
    bmip.resize(n, NULL);  
    bInflees.resize(n, NULL);
    bNumInflees.resize(n, 0);

    outPath.resize(n, NULL);  
    outAp.resize(n, NULL);  
    outB.resize(n, NULL);  

//    lowerHint = 0.0;
//    hap.resize(n, 0.0);
//    uap.resize(n, 0.0);

    cdds.resize(n);
    numCdds.resize(n, 0);
    longest = _longest;

    // re new
    maxHeap = new Heap(true, n);  // \Heap in paper
    lMaxHeap = new Heap(true, n);  // \Heap in paper
    used.resize(n, false);
    ap.resize(n, 0);
    sap.resize(n, 0);
    path.resize(n);
    children.resize(n);
    bb.resize(n);
    b.resize(n);
    numInfrs.resize(n, 0);
    S.resize(n, n);
    numchild.resize(n, 0);
    parent.resize(n, -1);
    heap.resize(n);
    queue.resize(n);
    distance.resize(n, longest);
    childlist.resize(n);
    status.resize(n, 0);
    is_native.resize(n, false);

    for (int i = 0; i < n; i++) {
        shortestPathsFrom(i);
//        if (PROTOTYPE) {
//            cout << i << " ===";
//            for (int j = 0; j < numInflees[i]; j++) {
//                cout << " " << inflees[i][j] << " : " << mip[i][j];
//            }
//            cout << endl;
//        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < numInflees[i]; j++) {
            int inflee = inflees[i][j];
            numCdds[inflee] += 1;
            cdds[inflee].push_back(i);
        }
    }
    cout << "initialization MIPs done." << endl;
}

void Mip::shortestPathsFrom(int node) {
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
    for (int i=0;i<top;i++)  // restore
        S[heap[i]] = n;

    bNumInflees[node] = bottom;
    numInflees[node] = bottom;
    mip[node] = new double[bottom];
    inflees[node] = new int[bottom];
    for (int i=0;i<bottom;i++)  // restore
    {
        int child=childlist[i];

        mip[node][i] = exp(-distance[child]);
        inflees[node][i] = child;

        distance[child]=longest;
        S[child]=n;
        parent[child]=-1;
    }
}

double Mip::incInfl(int node) {
//    cout << "incInfl -- " << node << endl;
    count++;
    double incInfl = 0.0;
    clock_t start = clock();

    int top=0, bottom=0;
    distance[node]=0;
    heap[0]=node;
    top++;
    parent[node] = node;
    b[node]=1;
    while (true){
        //pop out of heap
        if (distance[heap[0]]<longest) S[heap[0]]=-1; 
        else break;
        childlist[bottom++]=heap[0];
        for (int i=0;i<Graph::GetNeighbor(heap[0]);i++){
            Edge e=Graph::GetEdge(heap[0],i);
            if (used[e.v] || S[e.v]<0) continue;
            if (distance[e.v]>distance[heap[0]]+e.w1+EPS) {
                parent[e.v]=heap[0];
                b[e.v]=exp(-e.w1);
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

    // remember temporarily the out tree structure. 
    // avoid calculate it again later.
    outPath[node] = new int[bottom];
    outB[node] = new double[bottom];
    map<int, int> idx;
    for (int i=0;i<bottom;i++)  // restore
        idx[childlist[i]] = i;
    for (int i=0;i<bottom;i++) {
        outPath[node][i] = idx[parent[childlist[i]]];
        outB[node][i] = b[childlist[i]];
    }

    outAp[node] = new double[bottom];
    // calculate ap
    for (int i=0;i<bottom;i++)  
    {
        int child=childlist[i];
        if (!is_native[child]) continue;

        map<int, int> idx;
        int num = numInfrs[child];
        for (int j = 0; j < num; j++) {
            idx[children[child][j]] = j;
        }

        // add iter new tree
        int iter = child;
        while (true) {
            if (idx.count(iter) == 0) {  
                idx[iter] = num++;
                children[child].push_back(iter);
            }
            if (parent[iter] != iter)
                iter = parent[iter];
            else
                break;
        }

        // add iter's parentship relation
        path[child].resize(num, 0);
        path[child][0] = 0;
        bb[child].resize(num);
        bb[child][0] = 1.0;
        iter = child;
        while (true) {
            if (parent[iter] != iter) {
                path[child][idx[parent[iter]]] = idx[iter];
                bb[child][idx[parent[iter]]] = b[iter];
                iter = parent[iter];
            }
            else break;
        }

        // calculate ap
        int head = 0, tail = 0;

        for (int j = 0; j < num; j++)
            numchild[path[child][j]] += 1;
        for (int j = 0; j < num; j++) {
            int v = children[child][j];
            sap[j] = 1.0;
            if (used[v] || node == v) {
                queue[tail++] = j;
                numchild[j] = 0;
            }
        }
        numchild[idx[child]] -= 1; // avoid self loop

//        for (int j = 0; j < num; j++)
//            cout << path[child][j] << " ";
//        cout << endl;
//        for (int j = 0; j < num; j++)
//            cout << children[child][j] << " ";
//        cout << endl;
//        for (int j = 0; j < num; j++)
//            cout << bb[child][j] << " ";
//        cout << endl;
//        for (int j = 0; j < num; j++)
//            cout << numchild[j] << " ";
//        cout << endl;

        int x, u;
        while (head < tail) {
            x = queue[head++];
            u = path[child][x];
//            printf("x:%d u:%d\n", x, u);

            // avoid ineffective seed x
            int v = children[child][u];
            if (used[v] || node == v) continue;

            if (numchild[u] > 0)
                sap[u] *= (1 - sap[x]*bb[child][x]);
            if (!--numchild[u]) {
                sap[u] = 1 - sap[u];
                queue[tail++] = u;
            }
//            printf("sap[%d] = %f\n", u, sap[u]);
        }
        numchild[0] = 0;
        incInfl += sap[0] - ap[child];
        outAp[node][i] = sap[0];

        // restore
        path[child].resize(numInfrs[child]);
        children[child].resize(numInfrs[child]);
        bb[child].resize(numInfrs[child]);
        for (int j = 0; j < num; j++)
            sap[j] = 0.0;
//        cout << "sap " << child << " : " << sap[0] << endl;
    }

    for (int i=0;i<top;i++)  // restore
        S[heap[i]] = n;

    if (!changed[node]) {
        // keep the initial status
        changed[node] = true;
        int initNum = bNumInflees[node];
        bmip[node] = new double[initNum];
        bInflees[node] = new int[initNum];
        for (int j = 0; j < initNum; j++) {
            bmip[node][j] = mip[node][j];
            bInflees[node][j] = inflees[node][j];
        }
    }

    numInflees[node] = bottom;
    delete[] mip[node];
    delete[] inflees[node];
    mip[node] = new double[bottom];
    inflees[node] = new int[bottom];
    for (int i=0;i<bottom;i++)  // restore
    {
        int child=childlist[i];
        mip[node][i] = exp(-distance[child]);
        inflees[node][i] = child;
        distance[child]=longest;
        S[child]=n;
        parent[child]=-1;
    }

//    for (int i=0;i<bottom;i++) {
//        printf("node %d, ap[%d] == %f\n", node, childlist[i], outAp[node][i]);
//    }

    incInflTimer += (double)(clock() - start) / CLOCKS_PER_SEC;
    return incInfl;
}

// double Mip::spread(const unordered_set<int>& tset, const unordered_set<int>& influencees) {
//     double spread;
//     for (auto it= tset.begin(); it != tset.end(); it++) {
//         spread += incInfl(*it, seeds, influencees);
//     }
// }

void Mip::accept(int s) {
    used[s] = true;
    for (int i=0;i<numInflees[s];i++) {
        int child=inflees[s][i];
        ap[child] = outAp[s][i];
        // hint ap lower bound
//        if (ap[child] > hap[child]) {
//            hap[child] = ap[child];
//        }
//        // hint ap upper bound
//        if (ap[child] < uap[child]) {
//            upperHind += ap[child] - uap[child];
//            uap[child] = ap[child];
//        }

        map<int, int> idx;
        int num = numInfrs[child];
        for (int j = 0; j < num; j++) {
            idx[children[child][j]] = j;
        }

        // add iter new tree
        int j = i;
        while (true) {
            if (idx.count(inflees[s][j]) == 0) {  
                idx[inflees[s][j]] = num++;
                children[child].push_back(inflees[s][j]);
            }
            if (outPath[s][j] != j)
                j = outPath[s][j];
            else
                break;
        }
        numInfrs[child] = num;

        j = i;
        path[child].resize(num, 0);
        path[child][0] = 0;
        bb[child].resize(num);
        bb[child][0] = 1.0;
        while (true) {
            if (outPath[s][j] != j) {
                int parent = inflees[s][outPath[s][j]];
                path[child][idx[parent]] = idx[inflees[s][j]];
                bb[child][idx[parent]] = outB[s][j];
                j = outPath[s][j];
            }
            else break;
        }

        for (int j = 0; j < numCdds[child]; j++)
            status[cdds[child][j]] = 1;
    }

//    for (int i=0;i<numInflees[s];i++) {
//        printf("ap[%d] == %f\n", inflees[s][i], outAp[s][i]);
//    }

    delete[] outB[s];
    delete[] outPath[s];
    delete[] outAp[s];
}

int Mip::nextSeed(double *inc) {
    // status:
    // 0: init
    // 1: outdated
    // 2: estimated
    // 3: computed
    int s;
    while (true) {
        maxHeap->pop(&s, inc);

//        bool affected = false;
//        for (int j = 0; j < numInflees[s]; j++)
//            if (ap[inflees[s][j]] > 0)
//                affected = true;
//
//        if (!affected) {
//            accept(s);
//            topk.push_back(s);
//            cout << s << " not affected, inc: " << inc << endl;
//            break;
//        }
//        else if (status[s] == 0) {
        if (status[s] == 0) {
            incInfl(s);
//                cout << s << " use initial inc" << endl;
            return s;
        }
        else if (status[s] == 1) {
            double est = estimateIncInfl(s);
            estCount++;
            maxHeap->push(s, est);
            status[s] = 2;
            cout << s << " estimated inc: " << est << endl;
        }
        else if (status[s] == 2) {
            double com = incInfl(s);
 //           printf("compute: dp[%d] == %f\n", s, com);
            maxHeap->push(s, com);
            status[s] = 3;
        }
        else if (status[s] == 3) {
//            if (PROTOTYPE)
//                cout << s << " use computed inc" << endl;
            return s;
        }
    }
//    if (PROTOTYPE)
//        cout << i << "th seed: " << s <<
//            ", incremental influence: " << inc << endl;
}

vector<int> Mip::select(int k) {
    vector<int> topk;

    double timerCdd = 0.0;
    clock_t start = clock();
    maxHeap->clear();
    for (int i = 0; i < n; i++)
        used[i] = false;
    for (int i = 0; i < n; i++) {
        if (is_native[i]) {
            for (int j = 0; j < numCdds[i]; j++) {
                int u = cdds[i][j];
                if (used[u]) continue;
                used[u] = true;
                double inc = 0.0;
                for (int k = 0; k < numInflees[u]; k++) {
                    if (is_native[inflees[u][k]])
                        inc += mip[u][k];
                }
                if (inc > 0.0) {
                    maxHeap->push(u, inc);
                }
            }
        }
    }
    // careful, temporarily misuse used list. restore
    for (int i = 0; i < n; i++)
        used[i] = false;
//    for (int i = 0; i < n; i++) {
//        double inc = 0.0;
//        for (int j = 0; j < numInflees[i]; j++)
//            if (is_native[inflees[i][j]])
//                inc += mip[i][j];
//        if (inc > 0)
//            maxHeap->push(i, inc);
//    }
    timerCdd += (double)(clock() - start) / CLOCKS_PER_SEC;

//    if (PROTOTYPE)
//        cout << "1th seed: " << s <<
//            ", incremental influence: " << inc << endl;
    maxHeap->print();
    int s; double inc;
    for (int i = 1; i <= k; i++) {
        s = nextSeed(&inc);
        printf("%dth new seed: dp[%d] == %f\n", i, s, inc);
        maxHeap->print();
//        printf("%f\n", inc);
        accept(s);
        topk.push_back(s);
    }

    printf("expansion query find candidates operation %f seconds.\n", timerCdd);
    printf("expansion query heap operation %f seconds.\n", maxHeap->timer);

//    reset();
    return topk;
}

double Mip::estimateIncInfl(int s) {
    double est = 0.0;
    for (int j = 0; j < numInflees[s]; j++)
        if (is_native[inflees[s][j]])
            est += mip[s][j] * (1 - ap[inflees[s][j]]);
    return est;
}

double Mip::lowerIncInfl(int s) {
    double est = 0.0;
    for (int j = 0; j < numInflees[s]; j++)
        if (is_native[inflees[s][j]])
            if (mip[s][j] > ap[inflees[s][j]])
                est += mip[s][j] - ap[inflees[s][j]];
    return est;
}

// double Mip::estimateLowerIncInfl(const vector<int>& vertices) {
//     vector<double> hap(n, 0.0);
//     for (vector<int>::const_iterator it = vertices.begin(); it != vertices.end(); it++) {
//         double inc = 0.0;
//         for (int j = 0; j < numInflees[*it]; j++) {
//             int inflee = inflees[*it][j];
//             if (is_native[inflee]) {
//                 hap[inflee] = hap[inflee] > ap[inflee] ? hap[inflee]: ap[inflee];
//                 hap[inflee] = hap[inflee] > mip[*it][j] ? hap[inflee]: mip[*it][j];
//             }
//         }
//     }
// 
//     double est = 0.0;
//     for (vector<int>::const_iterator it = vertices.begin(); it != vertices.end(); it++) {
//         for (int j = 0; j < numInflees[*it]; j++) {
//             int inflee = inflees[*it][j];
//             if (is_native[inflee]) {
//                 est += hap[inflee];
//                 hap[inflee] = 0.0;
//             }
//         }
//     }
//     return est;
// }

double Mip::estimateLowerIncInfl(const vector<int>& vertices) {
    vector<double> hap(n, 0.0);
    double est = 0.0;
    for (vector<int>::const_iterator it = vertices.begin(); it != vertices.end(); it++) {
        double inc = 0;
        for (int j = 0; j < numInflees[*it]; j++) {
            int inflee = inflees[*it][j];
            if (is_native[inflee] && !(hap[inflee] > 0) && !(ap[inflee] > 0)) {
                est += mip[*it][j];
                inc += mip[*it][j];
                hap[inflee] = mip[*it][j];
            }
        }
        printf("%d's lowerBound really is %f\n", *it, inc);
    }
    return est;
}

// double Mip::estimateLowerIncInfl(const vector<int>& vertices) {
//     double est = 0.0;
//     for (vector<int>::const_iterator it = vertices.begin(); it != vertices.end(); it++) {
//         double inc = 0.0;
//         for (int j = 0; j < numInflees[*it]; j++) {
//             int inflee = inflees[*it][j];
//             if (is_native[inflee]) {
//                 if (mip[*it]pj]
//                 est -= hap[inflee];
//                 hap[inflee] = hap[inflee] > mip[*it][j] ? hap[inflee]: mip[*it][j];
//                 est += hap[inflee];
//                 inc += mip[*it][j] > hap[inflee] ? mip[*it][j] - hap[inflee] : 0;
//             }
//         }
//         printf("%d's lowerBound really is %f\n", *it, inc);
//     }
// 
//     return est;
// }

// double Mip::estimateUpperIncInfl(const vector<int>& vertices) {
//     double est = 0.0;
//     for (vector<int>::const_iterator it = vertices.begin(); it != vertices.end(); it++) {
//         for (int j = 0; j < numInflees[*it]; j++) {
//             int inflee = inflees[*it][j];
//             if (is_native[inflee]) {
//                 est -= uap[inflee];
//                 uap[inflee] = 1 - (1 - uap[inflee])*(1 - mip[*it][j]);
//                 est += uap[inflee];
//             }
//         }
//     }
// 
//     return est;
// }

void Mip::sortInit(const vector<int>& vertices, vector<int>& candidates,
        vector<double>& cddsDp)
{
    reset_native();
    maxHeap->clear();
    int numCandidates = 0;
//    for (vector<int>::const_iterator it = vertices.begin(); it !=
//            vertices.end(); it++)
//        is_native[*it] = true;
//    for (int i = 0; i < n; i++) {
//        double inc = 0.0;
//        for (int j = 0; j < numInflees[i]; j++)
//            if (is_native[inflees[i][j]])
//                inc += mip[i][j];
//        if (inc > 0 + EPS) {
//            numCandidates++;
//            maxHeap->push(i, inc);
//        }
//    }
//
//    for (int i = 0; i < n; i++)
//        used[i] = false;

//    vector<bool> used(n, false);
    for (vector<int>::const_iterator it = vertices.begin(); it !=
            vertices.end(); it++)
        is_native[*it] = true;

    for (vector<int>::const_iterator it = vertices.begin(); it !=
            vertices.end(); it++) {
        for (int j = 0; j < numCdds[*it]; j++) {
            int u = cdds[*it][j];
//            if (used[u]) continue;
//            used[u] = true;
            double inc = 0.0;
            for (int k = 0; k < numInflees[u]; k++) {
                if (is_native[inflees[u][k]])
                    inc += mip[u][k];
            }
            if (inc > 0.0) {
                maxHeap->push(u, inc);
//                numCandidates++;
            }
        }
    }
//    // careful, temporarily misuse used list. restore
//    for (int i = 0; i < n; i++)
//        used[i] = false;

    numCandidates = maxHeap->size();
    int s; double inc;
    int j = 0;
    candidates.resize(numCandidates);
    cddsDp.resize(numCandidates);
    while (maxHeap->pop(&s, &inc)) {
        candidates[j] = s;
        cddsDp[j] = inc;
        j++;
    }

    for (int i = 0; i < numCandidates; i++) {
        printf("list[%d] == %f\n", candidates[i], cddsDp[i]);
    }
    printf("\n");
    for (vector<int>::const_iterator it = vertices.begin(); it !=
            vertices.end(); it++)
        is_native[*it] = false;
}

vector<int> Mip::assembly(int k, const vector<vector<int> >& candidates, const
        vector<vector<double> > &initDp) {
    vector<int> topk;

    double timerCdd = 0.0;
    clock_t start = clock();

    // get random access
    vector<double> dp(n, 0.0);
    int numLists = candidates.size();
    vector<int> sizeList(numLists, 0);
    for (int i = 0; i < numLists; i++) {
        sizeList[i] = candidates[i].size();
//    // careful, temporarily misuse used list. Dangerous?
//        for (vector<int>::const_iterator it = candidates[i].begin(); it !=
//                candidates[i].end(); it++) {
//            if (used[*it]) continue;
//            used[*it] = true;
//            double init = 0.0;
//            for (int j = 0; j < numInflees[*it]; j++)
//                if (is_native[inflees[*it][j]])
//                    init += mip[*it][j];
//            dp[*it] = init;
//        }
        for (int j = 0; j < sizeList[i]; j++) {
            dp[candidates[i][j]] += initDp[i][j];
        }            
    }
    // careful, temporarily misuse used list. restore
//    for (int i = 0; i < n; i++)
//        used[i] = false;
    timerCdd += (double)(clock() - start) / CLOCKS_PER_SEC;

    cout << "number of lists: " << numLists << endl;

    // print out for prototype
//    printf("candidates init dp (TA)\n");
//    for (int i = 0; i < candidates.size(); i++) {
//        printf("region #%d\n", i);
//        for (int j = 0; j < candidates[i].size(); j++) {
//            printf("dp[%d] == %f\n", candidates[i][j], initDp[i][j]);
//        }
//    }

    /*
       TA select top 1
       */
    vector<bool> seen(n, false);
    int s;  // new seed
    double t = INFINITY;  // threshold value
    double B_H = 0.0;  // bound
    int d = 0;  // depth
    while (t > B_H) {
        // sorted access
        t = 0.0;
        for (int i = 0; i < numLists; i++) {
            // if run out of list
            if (d >= sizeList[i]) continue;
                
            int v = candidates[i][d];
            if (!seen[v]) {
                maxHeap->push(v, dp[v]);
                seen[v] = true;
                printf("push vertex: dp[%d] == %f\n", v, dp[v]);
            }
            t += initDp[i][d];
        }
        printf("%d depth, threshold value %f, bound %f\n", d, t, B_H);
        d += 1;
        maxHeap->top(&s, &B_H);
        maxHeap->print();
    }
    s = nextSeed(&B_H);
    accept(s);
    printf("new seed: dp[%d] == %f\n", s, B_H);
    topk.push_back(s);
    maxHeap->print();

    /*
       TA select top 2-k
       */
    int _k = 2;
    while (_k <= k) {
        s = nextSeed(&B_H);
        while (t > B_H) {
            // sorted access
            t = 0.0;
            for (int i = 0; i < numLists; i++) {
                // if run out of list
                if (d >= sizeList[i]) continue;
                    
                int v = candidates[i][d];
                
                if (!seen[v]) {
                    if (status[v] == 0) {
                        maxHeap->push(v, dp[v]);
                        printf("push vertex: dp[%d] == %f\n", v, dp[v]);
                    }
                    else if (status[v] == 1) {
                        double est = estimateIncInfl(v);
                        estCount++;
                        maxHeap->push(v, est);
                        status[v] = 2;
                        printf("push vertex: est[%d] == %f\n", v, est);
                    }
                    seen[v] = true;
                }
                maxHeap->print();
                t += initDp[i][d];
            }

            printf("%d depth, threshold value %f, s %d, bound %f\n", d, t, s, B_H);
            d += 1;
        }
        status[s] = 3;
        maxHeap->push(s, B_H);
        s = nextSeed(&B_H);
        printf("%dth new seed: dp[%d] == %f, threshold == %f\n", _k, s, B_H, t);
        accept(s);
        _k += 1;
        topk.push_back(s);
        maxHeap->print();
    }

    printf("assembly query find candidates operation %f seconds.\n", timerCdd);
    printf("assembly query heap operation %f seconds.\n", maxHeap->timer);

//    reset();
    return topk;
}

void Mip::cacheTau(int k, vector<int>& topTau, vector<double>& tauDp)
{
    for (int i = 0; i < n; i++) {
        double inc = 0.0;
        for (int j = 0; j < numInflees[i]; j++)
            if (is_native[inflees[i][j]])
                inc += mip[i][j];
        if (inc > 0)
            maxHeap->push(i, inc);
    }

    int s; double inc;
    for (int i = 1; i <= k; i++) {
        s = nextSeed(&inc);
//        printf("%dth new seed: dp[%d] == %f\n", i, s, inc);
        accept(s);
        topTau.push_back(s);
        tauDp.push_back(inc);
    }
//    reset();
}

double Mip::getHints(vector<int>&hint, const vector<vector<int> >& topTau, const vector<vector<double> >& dpTau, int num)
{   // return the estimated influence of topTau
    set<int> seeds;
    double infl = 0.0f;
    int numOfRegions = topTau.size();
    vector<int> positions(numOfRegions, 0);
    vector<int> tauSize(numOfRegions);
    for (int i = 0; i < numOfRegions; i++) {
        tauSize[i] = topTau[i].size();
    }
    int count = 0;
    while (count < num) {
        double max = 0;
        int mp = -1; 
        int r = -1;
        for (int i = 0; i < numOfRegions; i++) {
            // run out of this region.
            if (positions[i] > tauSize[i] -1) continue;

            int t = topTau[i][positions[i]];
            while (positions[i] < tauSize[i] && seeds.count(t) > 0) {
                infl += dpTau[i][positions[i]];
//                hintDp[t] += dpTau[i][positions[i]];
                positions[i] += 1;
                t = topTau[i][positions[i]];
            }
            if (positions[i] > tauSize[i] -1) continue;

            double tmp = dpTau[i][positions[i]];
            if (max < tmp) {
                max = tmp;
                mp = t;
                r = i;
            }
        }
//        printf("mp == %d, max == %f\n", mp, max);
        if (mp < 0)
            break;
        else {
            if (positions[r] < tauSize[r])
                positions[r] += 1;
            infl += max;
            seeds.insert(mp);
            hint.push_back(mp);
            count += 1;
        }
    }

    return infl;
}

vector<int> Mip::dumbHint(int k, double acceptRatio, const vector<vector<int> >& candidates, const vector<vector<double> > &initDp) 
{
    vector<int> topk(k);
    clock_t start = clock();

    double upperBound = 0.0;
    vector<int> holder(k);  // hold temporary seeds used by upperBound
    // hold temporary seeds' incremental influence used by upperBound
    vector<double> dpHolder(k);  
    double lowerBound = 0.0;
    vector<int> lowers(k);

    double seedInfl = 0.0;

    double timerCdd = 0.0;
    start = clock();
    // get random access
    // careful, temporarily misuse used list. Dangerous?
    vector<double> dp(n, 0.0);
    int numLists = candidates.size();
    vector<int> sizeList(numLists, 0);
    for (int i = 0; i < numLists; i++) {
        sizeList[i] = candidates[i].size();
        for (int j = 0; j < sizeList[i]; j++) {
            dp[candidates[i][j]] += initDp[i][j];
        }            
    }

    timerCdd += (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("dumb hint query find candidates operation %f seconds.\n", timerCdd);

    int _k = 1;
    /*
       TA select top 1
       */
    vector<bool> seen(n, false);
    int s;  // new seed
    double t = INFINITY;  // threshold value
    double B_H = 0.0;  // bound
    int d = 0;  // depth
    while (t > B_H) {
        // sorted access
        t = 0.0;
        for (int i = 0; i < numLists; i++) {
            // if run out of list
            if (d >= sizeList[i]) continue;
                
            int v = candidates[i][d];
            if (!seen[v]) {
                maxHeap->push(v, dp[v]);
                seen[v] = true;
 //               printf("push vertex: dp[%d] == %f\n", v, dp[v]);
            }
            t += initDp[i][d];
        }
//        printf("%d depth, threshold value %f, bound %f\n", d, t, B_H);
        d += 1;
        maxHeap->top(&s, &B_H);
    }
    s = nextSeed(&B_H);
    accept(s);
    printf("new seed: dp[%d] == %f\n", s, B_H);
    topk[0] = s;
    lowers[0] = s;
    seedInfl += B_H;

//    upperBound = seedInfl + (k - _k) * B_H;
//    if (seedInfl / upperBound > acceptRatio) {
////        reset();
//        return topk;
//    }

    /*
       TA select top 2-k
       */
    _k = 2;
    while (_k <= k) {
        s = nextSeed(&B_H);
        while (t > B_H) {
            // sorted access
            t = 0.0;
            for (int i = 0; i < numLists; i++) {
                // if run out of list
                if (d >= sizeList[i]) continue;
                    
                int v = candidates[i][d];
                
                if (!seen[v]) {
                    if (status[v] == 0) {
                        maxHeap->push(v, dp[v]);
//                        printf("push vertex: dp[%d] == %f\n", v, dp[v]);
                    }
                    else if (status[v] == 1) {
                        double est = estimateIncInfl(v);
                        estCount++;
                        maxHeap->push(v, est);
                        status[v] = 2;
 //                       printf("push vertex: est[%d] == %f\n", v, est);
                    }
                    seen[v] = true;
                }

                t += initDp[i][d];
            }
//            printf("%d depth, threshold value %f, s %d, bound %f\n", d, t, s, B_H);
            d += 1;
        }
        status[s] = 3;
        maxHeap->push(s, B_H);
        s = nextSeed(&B_H);
        printf("%dth new seed: dp[%d] == %f, threshold == %f\n", _k, s, B_H, t);
        accept(s);
        topk[_k-1] = s;
        lowers[_k-1] = s;
        seedInfl += B_H;

        // calculate lower bound
        lMaxHeap->clear();
        for (int i = 0; i < n; i++) {
            if (!used[i]) {
                lMaxHeap->push(i, lowerIncInfl(i));
            }
        }
        int j = 0;
        int l; double low;
        while (j < k - _k) {
            lMaxHeap->pop(&l, &low);
            printf("%d's lower bound %f\n", l, low);
            lowers[_k + j] = l;
            j++;
        }
        lowerBound = seedInfl + estimateLowerIncInfl(lowers);
        printf("lower bound bound %f\n", lowerBound);

        // calculate upperBound
        j = 0;
        upperBound = seedInfl; // + (k - _k) * B_H;
        while (j < k - _k) {
            int t; double tp;
            if (maxHeap->pop(&t, &tp)) {
                if (status[t] == 1) {
                    double est = estimateIncInfl(t);
                    maxHeap->push(t, est);
                    status[t] = 2;
                }
                else {
                    holder[j] = t;
                    dpHolder[j] = tp;
                    upperBound += dpHolder[j];
                    printf("compliment upper bound %d, %f\n", holder[j], dpHolder[j]);
                    j++;
                }
            }
            else {
                holder[j] = holder[j-1];
                dpHolder[j] = dpHolder[j-1];
                upperBound += dpHolder[j];
                printf("compliment upper bound %d, %f\n", holder[j], dpHolder[j]);
                j++;
            }
        }
        j = 0;
        while (j < k - _k) {  // restore
            maxHeap->push(holder[j], dpHolder[j]);
            j++;
        }

        printf("%dth lowerBound: %f, upperBound: %f, accept rate: %f\n", _k, lowerBound, upperBound, acceptRatio);
//        upperBound = seedInfl + (k - _k) * B_H;
        if (lowerBound / upperBound > acceptRatio) {
            printf("%dth lowerBound: %f, upperBound: %f, accept rate: %f\n", _k, lowerBound, upperBound, acceptRatio);
            int next; double  _;
            while(_k < k && maxHeap->pop(&next, &_)) {
                topk[_k++] = next;
            }
//            reset();
            return topk;
        }

        _k += 1;
    }

//    reset();
    return topk;
}

vector<int> Mip::hint(int k, double acceptRatio, const vector<vector<int> >& candidates, const vector<vector<double> > &initDp, const vector<vector<int> >& topTau, const vector<vector<double> > &dpTau) 
{
    vector<int> topk(k);
//    int step = k / 100 == 0 ? 1 : k / 100;
    int step = 1;

    vector<int> hints;
    vector<int> holder(k);  // hold temporary seeds used by upperBound
    // hold temporary seeds' incremental influence used by upperBound
    vector<double> dpHolder(k);  

    clock_t start = clock();
    double timerHints = 0.0;
    double lowerHint = getHints(hints, topTau, dpTau, k);
//    int hintUsed = hints.size();
    timerHints = (double)(clock() - start) / CLOCKS_PER_SEC;

    printf("hints are ");
    for (int i = 0; i < hints.size(); i++)
        printf("%d ", hints[i]);
    printf("\n");

    for (int i = 0; i < candidates.size(); i++) {
        printf("Region #%d\n", i);
        for (int j = 0; j < candidates[i].size(); j++) {
            printf("init-dp[%d] = %f\n", candidates[i][j], initDp[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    for (int i = 0; i < topTau.size(); i++) {
        printf("top tau #%d\n", i);
        for (int j = 0; j < topTau[i].size(); j++) {
            printf("dp[%d] = %f\n", topTau[i][j], dpTau[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    double upperBound = 0.0;
    double seedInfl = 0.0;
//    lowerHint = estimateLowerIncInfl(hints);

    double timerCdd = 0.0;
    start = clock();
    // get random access
    // careful, temporarily misuse used list. Dangerous?
    vector<double> dp(n, 0.0);
    int numLists = candidates.size();
    vector<int> sizeList(numLists, 0);
    for (int i = 0; i < numLists; i++) {
        sizeList[i] = candidates[i].size();
        for (int j = 0; j < sizeList[i]; j++) {
            dp[candidates[i][j]] += initDp[i][j];
        }            
    }

    timerCdd += (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("hint query find candidates operation %f seconds.\n", timerCdd);

    int _k = 1;
    /*
       TA select top 1
       */
    vector<bool> seen(n, false);
    int s;  // new seed
    double t = INFINITY;  // threshold value
    double B_H = 0.0;  // bound
    int d = 0;  // depth
    while (t > B_H) {
        // sorted access
        t = 0.0;
        for (int i = 0; i < numLists; i++) {
            // if run out of list
            if (d >= sizeList[i]) continue;
                
            int v = candidates[i][d];
            if (!seen[v]) {
                maxHeap->push(v, dp[v]);
                seen[v] = true;
 //               printf("push vertex: dp[%d] == %f\n", v, dp[v]);
            }
            t += initDp[i][d];
        }
//        printf("%d depth, threshold value %f, bound %f\n", d, t, B_H);
        d += 1;
        maxHeap->top(&s, &B_H);
    }
    s = nextSeed(&B_H);
    accept(s);
    printf("new seed: dp[%d] == %f\n", s, B_H);
    topk[0] = s;
    seedInfl += B_H;

    int j = 0;
    // calculate upperBound
    upperBound = seedInfl; // + (k - _k) * B_H;
    maxHeap->print();
    j = 0;
    while (j < k - _k) {
        int t; double tp;
        if (maxHeap->pop(&t, &tp)) {
            if (status[t] == 1) {
                double est = estimateIncInfl(t);
                maxHeap->push(t, est);
                status[t] = 2;
            }
            else {
                holder[j] = t;
                dpHolder[j] = tp;
                upperBound += dpHolder[j];
                printf("compliment upper bound %d, %f\n", holder[j], dpHolder[j]);
                j++;
            }
        }
        else {
            holder[j] = holder[j-1];
            dpHolder[j] = dpHolder[j-1];
            upperBound += dpHolder[j];
            printf("compliment upper bound %d, %f\n", holder[j], dpHolder[j]);
            j++;
        }
    }
    j = 0;
    while (j < k - _k) {  // restore
        maxHeap->push(holder[j], dpHolder[j]);
        j++;
    }
    maxHeap->print();

    timerHints += (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("%dth hint lowerBound: %f, upperBound: %f, accept rate: %f\n", _k, lowerHint, upperBound, acceptRatio);

//    start = clock();
//    // calculate upperBound
//    upperBound = seedInfl; // + (k - _k) * B_H;
//    while (j < k - _k) {
//        if (!maxHeap->pop(&(holder[j]), &(dpHolder[j]))) {
//            holder[j] = holder[j-1];
//            dpHolder[j] = dpHolder[j-1];
//        }
//        upperBound += dpHolder[j];
//        j++;
//    }
//    j = 0;
//    while (j < k - _k) {  // restore
//        maxHeap->push(holder[j], dpHolder[j]);
//        j++;
//    }
//    // calculate lowerBound
    int countHint = 0;
    int iterHint = 0;
//    for (iterHint = 0; countHint < k - _k && iterHint < k; iterHint++) {
//        if (!used[hints[iterHint]]) {
//            topk[_k + countHint] = hints[iterHint];
//            countHint++;
//        }
//    }
//    while (iterHint < hintUsed) {  // new seed is not in current hints.
//        int aborted = hints[hintUsed--];
//        for (int j = 0; j < numInflees[aborted]; j++) {
//            int inflee = inflees[aborted][j];
//            if (is_native[inflee]) {
//                if (fabs(hap[inflee] - mip[aborted][j]) < EPS) {
//                    lowerHint -= hap[inflee];
//                    hap[inflee] = 0.0;
//                }
//            }
//        }
//    }
//    timerHints += (double)(clock() - start) / CLOCKS_PER_SEC;
//
//    printf("1th lowerBound: %f, upperBound: %f\n", lowerHint, upperBound);
//    
//    if (lowerHint / upperBound > acceptRatio) {
////        reset();
//        printf("collect hints %f seconds\n", timerHints);
//        return topk;
//    }

    /*
       TA select top 2-k
       */
    _k = 2;
    while (_k <= k) {
        s = nextSeed(&B_H);
        while (t > B_H) {
            // sorted access
            t = 0.0;
            for (int i = 0; i < numLists; i++) {
                // if run out of list
                if (d >= sizeList[i]) continue;
                    
                int v = candidates[i][d];
                
                if (!seen[v]) {
                    if (status[v] == 0) {
                        maxHeap->push(v, dp[v]);
//                        printf("push vertex: dp[%d] == %f\n", v, dp[v]);
                    }
                    else if (status[v] == 1) {
                        double est = estimateIncInfl(v);
                        estCount++;
                        maxHeap->push(v, est);
                        status[v] = 2;
 //                       printf("push vertex: est[%d] == %f\n", v, est);
                    }
                    seen[v] = true;
                }

                t += initDp[i][d];
            }
//            printf("%d depth, threshold value %f, s %d, bound %f\n", d, t, s, B_H);
            d += 1;
        }
        status[s] = 3;
        maxHeap->push(s, B_H);
        s = nextSeed(&B_H);
        printf("%dth new seed: dp[%d] == %f, threshold == %f\n", _k, s, B_H, t);
        accept(s);
        topk[_k-1] = s;
        seedInfl += B_H;

        if (_k % step == 0) {
            start = clock();
            // calculate upperBound
            upperBound = seedInfl; // + (k - _k) * B_H;
            maxHeap->print();
            j = 0;
            while (j < k - _k) {
                int t; double tp;
                if (maxHeap->pop(&t, &tp)) {
                    if (status[t] == 1) {
                        double est = estimateIncInfl(t);
                        maxHeap->push(t, est);
                        status[t] = 2;
                    }
                    else {
                        holder[j] = t;
                        dpHolder[j] = tp;
                        upperBound += dpHolder[j];
                        printf("compliment upper bound %d, %f\n", holder[j], dpHolder[j]);
                        j++;
                    }
                }
                else {
                    holder[j] = holder[j-1];
                    dpHolder[j] = dpHolder[j-1];
                    upperBound += dpHolder[j];
                    printf("compliment upper bound %d, %f\n", holder[j], dpHolder[j]);
                    j++;
                }
            }
            j = 0;
            while (j < k - _k) {  // restore
                maxHeap->push(holder[j], dpHolder[j]);
                j++;
            }
            maxHeap->print();

//            // calculate lowerBound
//            countHint = 0;
//            iterHint = 0;
//            for (iterHint = 0; countHint < k - _k && iterHint < k; iterHint++) {
//                if (!used[hints[iterHint]]) {
//                    topk[_k + countHint] = hints[iterHint];
//                    countHint++;
//                }
//            }
//            while (hintUsed > 0 && iterHint < hintUsed) {  // new seed is not in current hints.
//                int aborted = hints[--hintUsed];
//                for (int j = 0; j < numInflees[aborted]; j++) {
//                    int inflee = inflees[aborted][j];
//                    if (is_native[inflee]) {
//                        if (fabs(hap[inflee] - mip[aborted][j]) < EPS) {
//                            lowerHint -= hap[inflee];
//                            hap[inflee] = 0.0;
//                        }
//                    }
//                }
//            }
            timerHints += (double)(clock() - start) / CLOCKS_PER_SEC;
            printf("%dth hint lowerBound: %f, upperBound: %f, accept rate: %f\n", _k, lowerHint, upperBound, acceptRatio);

            if (lowerHint / upperBound > acceptRatio) {
                printf("%dth hint lowerBound: %f, upperBound: %f, accept rate: %f\n", _k, lowerHint, upperBound, acceptRatio);
                printf("collect hints %f seconds\n", timerHints);
//                reset();
                return hints;
            }
        }

        _k += 1;
    }

    printf("collect hints %f seconds\n", timerHints);
//    reset();
    return topk;
}

void Mip::reset() {
    clock_t start = clock();
    count = 0;
    incInflTimer = 0;
    estCount = 0;
    n = Graph::GetN();
    maxHeap->clear();
//    lowerHint = 0.0;

    for (int i = 0; i < n; i++) {
        used[i] = false;
        ap[i] = 0.0;
//        hap[i] = 0.0;
//        uap[i] = 0.0;

        // reset influencees's influence paths from seeds
        path[i].clear();
        children[i].clear();
        bb[i].clear();
        b[i] = 0.0;
        numInfrs[i] = 0;

        status[i] = 0;

        if (changed[i]) {
            numInflees[i] = bNumInflees[i];
            delete[] mip[i];
            delete[] inflees[i];
            mip[i] = new double[numInflees[i]];
            inflees[i] = new int[numInflees[i]];
            for (int j = 0; j < bNumInflees[i]; j++) {
                mip[i][j] = bmip[i][j];
                inflees[i][j] = bInflees[i][j];
            }
            delete[] bmip[i];
            delete[] bInflees[i];
        }
        changed[i] = false;
//      // memory access error to delete not exist 
//        if (outPath[i] != NULL) {
//            delete[] outPath[i];
//            delete[] outAp[i];
//            delete[] outB[i];
//        }
    }
    printf("reset takes %f seconds\n", (double)(clock() - start) / CLOCKS_PER_SEC);
}

void Mip::gc() {
    delete maxHeap;
    for (int i = 0; i < n; i++) {
        delete[] mip[i];
        delete[] inflees[i];
//        if (outPath[i] != NULL) {
//            delete[] outPath[i];
//            delete[] outAp[i];
//            delete[] outB[i];
 //       }
    }
}

void Mip::exit() {
    gc();
}

void Mip::reset_native()
{
    is_native.resize(n);
    for (int _ = 0; _ < n; _++) {
        is_native[_] = false;
    }
}
//    maxHeap->clear();
//    minHeap->clear();
