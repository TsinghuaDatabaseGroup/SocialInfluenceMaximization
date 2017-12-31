#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>    // std::max
#include <time.h>    // std::max

#include "tim.h"

#define INF 2147483647

double* TIM::ap_;
double* TIM::tap_;
vector<HNode> TIM::H_;
Q TIM::q_;
double TIM::theta_;
double TIM::theta2_;
int TIM::round_;
double *TIM::dist_;
int *TIM::seen_;
int *TIM::seen_idx_;
int *TIM::children_;
double *TIM::infl2_;
bool *TIM::is_border_;
int *TIM::pred_;

vector<vector<di> > TIM::iees_cache_;
vector<di> TIM::L_;
int TIM::cursorL_;
int TIM::n_;
bool *TIM::used_;
double *TIM::lap_;
vector<UBSample> TIM::ubsamples_;
vector<LBSample> TIM::lbsamples_;
// vector<vector<di> > TIM::Iee_;
// vector<int> TIM::pre_;
// vector<set<int> > TIM::In_;

double exactTimer; clock_t eStart;
double boundTimer; clock_t bStart;
double dTimer; clock_t dStart; double dSize;
double apTimer; clock_t apStart;
double pTimer; clock_t pStart;
int boundCounter;
int exactCounter;
int dCounter;
clock_t start;
double tTimer; clock_t tStart;

bool tim_min_heap_cmp(const di&a, const di&b) {
    return a.first > b.first;
}

bool max_candidate_cmp(const HNode&a, const HNode&b) {
    return a.infl < b.infl;
}

void TIM::Reset() {
#ifdef VERBOSE
    printf("Reseting TIM workspace ... Reset is needed for each query.\n");
#endif

    TIC::resetCache();

    round_ = 0;

    H_.clear();
    cursorL_ = 0;

//    Iee_.resize(n_);
    iees_cache_.resize(n_);

//    pre_.resize(n_);
//    In_.resize(n_);
    for (unsigned int i = 0; i < n_; i++) ap_[i] = 0;
    for (unsigned int i = 0; i < n_; i++) tap_[i] = 0;

    // reset timer 
    exactCounter = 0;
    boundCounter = 0;
    dCounter = 0;
    exactTimer = 0;
//    boundTimer = 0;
    dTimer = 0;
    apTimer = 0;
    pTimer = 0;
    dSize = 0;
    tTimer = 0;  // test
}

void TIM::Init() {
    theta_ = 0.01;
//    theta2_ = 0.1;
    L_.clear();
    n_ = TIC::GetN();
    infl2_ = new double[n_];

    ap_ = new double[n_];
    tap_ = new double[n_];
    used_ = new bool[n_];
    for (unsigned int i = 0; i < n_; i++) used_[i] = false;

    dist_ = new double[n_];
    seen_ = new int[n_];
    seen_idx_ = new int[n_];
    children_ = new int[n_];
    pred_ = new int[n_];
    is_border_ = new bool[n_];
    for (unsigned int i = 0; i < n_; i++) dist_[i] = INF;
    for (unsigned int i = 0; i < n_; i++) seen_idx_[i] = n_;
    for (unsigned int i = 0; i < n_; i++) is_border_[i] = true;

#ifdef VERBOSE
    PrintArguments();
#endif
}

void TIM::PrintArguments() {
    printf("N:%d M:%d Z:%d theta:%lg theta2 %lg\n", n_, TIC::GetM(), TIC::GetZ(), theta_, theta2_);
}

void TIM::InitGreedy(char* filename) {
    Init();
#ifdef VERBOSE
    printf("Initializing Greedy BestFirst framework...\n");
#endif

    // Load \List
    FILE *f = fopen(filename, "r");
    int n, id; double d, d2;
    fscanf(f, "%d", &n);
    assert(n == n_);
    
    for (int i = 0; i < n; i++) {
        fscanf(f, "%d %lg %lg", &id, &d, &d2);
        L_.push_back(di(d, id));
        infl2_[id] = d2;
    }

    fclose(f);
//    printf("|L_| = %d\n", L_.size());
}

void TIM::InitApprox(char *gfile, char *lfile, char *ufile) {
    InitGreedy(gfile);
#ifdef VERBOSE
    printf("Initializing Approximation framework...\n");
#endif
    UBLoad(lfile);
    LBLoad(ufile);
#ifdef VERBOSE
    printf("Initializing Approximation framework done.\n");
#endif
}

// void TIM::InitGreedyUpApprox() {
// }

void PrintQ(Q q) {
    printf("Q = (%d, {", q.k);
    for (unsigned int i = 0; i < q.item.size(); i++) {
        printf("%lg, ", q.item[i]);
    }
    printf("})\n");
}

S TIM::Approx(Q q, double epsilon) {
    Reset();
    q_ = q;
//    assert(q.item.size() == TIC::GetZ());
#ifdef VERBOSE
    PrintQ(q_);
#endif

    S rsp;
    di s;

    lap_ = new double[n_];
    for (int i = 0; i < n_; i++) lap_[i] = 0;

    vector<int> candidates = LBSelect();
    double UB = UBCalculate(), LB = 0, Real = 0;
    int select_idx = -1;  /* the index of candidate that are selected by Greedy*/
    int step = q_.k/10;
    for (round_ = 1; round_ <= q_.k; round_++) {
        if (round_%step == 0) {
            LB = Real + LBCalculate(candidates);
#ifdef VERBOSE
            printf("round %d, Ubound %lg, Lbound %lg\n", round_, UB, LB);
#endif
        }
//        printf("candidates:\n");
//        for (unsigned int i = 0; i < candidates.size(); i++) {
//            printf("%d ", candidates[i]);
//        }
//        printf("\n");
        if (LB > UB*epsilon) {
            for (unsigned int i = 0; i < candidates.size(); i++) {
                rsp.ids.push_back(candidates[i]);
                rsp.minflu.push_back(0);
            }

            for (int i = 0; i < rsp.ids.size(); i++) {
                // restore
                used_[rsp.ids[i]] = false;
            }
            return rsp;
        }

        s = BestFirst();
        rsp.ids.push_back(s.second);
        rsp.minflu.push_back(s.first);
        Real += s.first;
        UpdateAP(s.second);
        used_[s.second] = true;
//        printf("seed %d, marginal %lg\n", s.second, s.first);

        for (unsigned int i = 0; i < candidates.size(); i++)
            if (candidates[i] == s.second) select_idx = i;
        if (!(select_idx < 0)) {
            for (unsigned int i = select_idx; i < candidates.size()-1; i++) {
                candidates[i] = candidates[i+1];
            }
        }
        candidates.pop_back();
        select_idx = -1;  /* restore */
    }

    for (int i = 0; i < rsp.ids.size(); i++) {
        // restore
        used_[rsp.ids[i]] = false;
    }

    return rsp;
}

S TIM::Greedy(Q q) {
    return GreedyBound(q, true);
}

S TIM::GreedyNoBound(Q q) {
    return GreedyBound(q, false);
}

S TIM::GreedyBound(Q q, bool bounded) {
    start = clock();
    Reset();
#ifdef VERBOSE
    printf("reset takes time %lg s.\n", double(clock() - start) / CLOCKS_PER_SEC);
#endif
//    assert(q.item.size() == TIC::GetZ());
    q_ = q;
#ifdef VERBOSE
    PrintQ(q_);
#endif

    S rsp;
    di s;
    for (round_ = 1; round_ <= q_.k; round_++) {
        if (bounded) s = BestFirst();
        else s = BestFirstNoBound();
        rsp.ids.push_back(s.second);
        rsp.minflu.push_back(s.first);
//        tStart = clock();  // not here
        UpdateAP(s.second);
        used_[s.second] = true;
#ifdef VERBOSE
        printf("round %d, seed %d, marginal %lg\n", round_, s.second, s.first);
#endif
//        tTimer += double(clock()-tStart);
    }

#ifdef VERBOSE
//    printf("bound takes time %lg s, %d times.\n", boundTimer/CLOCKS_PER_SEC, boundCounter);
//    printf("exact takes time %lg s, %d times.\n", exactTimer/CLOCKS_PER_SEC, exactCounter);
    printf("bound takes %d times.\n", boundCounter);
    printf("exact takes %d times.\n", exactCounter);
//    printf("Dijkstra time %lg s, %d times, average size %lg.\n", dTimer/CLOCKS_PER_SEC, dCounter, dSize/dCounter);
//    printf("Dijkstra push time %lg s.\n", tTimer/CLOCKS_PER_SEC);
//    printf("marginal ap takes time %lg s.\n", apTimer/CLOCKS_PER_SEC);
//    printf("PProb takes time %lg s.\n", pTimer/CLOCKS_PER_SEC);
//    printf("UpdateAP takes time %lg s.\n", tTimer/CLOCKS_PER_SEC);
//    printf("H_ operation takes time %lg s.\n", tTimer/CLOCKS_PER_SEC);
#endif
#ifdef EXP
     printf("%d, %d, %d ", boundCounter, exactCounter, round_);
#endif
    
    for (int i = 0; i < rsp.ids.size(); i++) {
        // restore
        used_[rsp.ids[i]] = false;
    }

    return rsp;
}

void TIM::UpdateAP(const int u) {
    int v;
    vector<di>& rst = iees_cache_[u];
    for (unsigned int i = 0; i < rst.size(); i++) {
        v = rst[i].second;
        ap_[v] += rst[i].first;
        tap_[v] = ap_[v];
//        In_[v].insert(rst[pre[i]].second);
    }
}

double TIM::MarginalAPOf(const int v) {
//    apStart = clock();

    if (used_[v]) return 0;
    double ap_v = 1.0;
    for (unsigned int i = 0; i < TIC::GetInNeighbor(v); i++) {
        Edge e = TIC::GetInEdge(v, i);
        if (tap_[e.u] > 0) {
            ap_v *= (1 - TIC::GetInEdgeWeight(v, i, q_.item)*tap_[e.u]);
        }
    }
//    printf("In_[%d]={", v);
//    for (set<int>::iterator it = In_[v].begin(); it != In_[v].end(); it++) {
//        printf("%d ", *it);
//    }
//    printf("}\n");
    ap_v = 1.0 - ap_v;
    tap_[v] = ap_v;
#ifdef VERBOSE
    if (ap_v < ap_[v]) {
        printf("ap(%d)=%lg, tap=%lg\n", v, ap_[v], tap_[v]);
    }
    assert (ap_v >= ap_[v]);
#endif
    
//    apTimer += double(clock() - apStart);

    return ap_v - ap_[v];
}

// S TIM::GreedyUpApprox(Q q) {
//     Reset();
//    assert(q.item.size() == TIC::GetZ());
// }

int TIM::Dijkstra(int u, double max) {
//    printf("Dijkstra(%d, %lg)\n", u, max);
//    dCounter += 1;
//    dStart = clock();

//    vector<di> rst;
//    vector<int> pre;
    dist_[u] = 0;
    int top = 0;
    int bottom = 0;
    seen_[top++] = u;
    children_[bottom++] = u;
    seen_idx_[u] = -1;

//    rst.push_back(di(1.0f, u));
//    pre.push_back(0);

    // variables used in the while loop.
    double w_dist; int v, w;
    double tmp;

    while (top > 0) {
        v = seen_[0];
//        printf("Dijkstra pop (%d, %lg)\n", v, dist_[v]);

        if (dist_[v] < max) {
            seen_idx_[v] = -1;
        } else break;

        // shortest distance vertex for sure.
        if (v != u) {
//            rst.push_back(di(exp(-dist_[v]), v));
            children_[bottom++] = v;
            is_border_[pred_[v]] = false;
        }

        for (unsigned int i = 0; i < TIC::GetOutNeighbor(v); i++) {
            Edge e = TIC::GetOutEdge(v, i);
//            printf("e(%d, %d)\n", e.u, e.v);
            w = e.v;

            if (used_[w] || seen_idx_[w] < 0) continue;

//            pStart = clock();
            w_dist = -log(TIC::GetOutEdgeWeight(v, i, q_.item));
//            pTimer += double(clock() - pStart);

//            tStart = clock();
            if (w_dist + dist_[v] < dist_[w]) {
                dist_[w] = w_dist + dist_[v];
                pred_[w] = v;
                // push
                int j;
                if (seen_idx_[w] >= n_) {  // not seen before
                    seen_[top] = w;
                    j = top++;
                } else {
                    j = seen_idx_[w];
                }
                int x = (j-1)/2;
                tmp = dist_[seen_[j]];
                while (j > 0) {
                    if (dist_[seen_[x]] > tmp) {
                        seen_[j] = seen_[x];
                        if (seen_idx_[seen_[j]] < n_) seen_idx_[seen_[j]] = j;
                        j = x; x = (j-1)/2;
                    } else break;
                }
                seen_[j] = w;
                seen_idx_[w] = j;
            }
 //           tTimer += double(clock() - tStart);
        }

        // pop
        seen_[0] = seen_[--top];
        if (!top) break;
        int j = 0, x = j*2 + 1;
        tmp = dist_[seen_[j]];
        while (x < top) {
            if (x+1 < top && dist_[seen_[x+1]] < dist_[seen_[x]]) x=x+1;
            if (dist_[seen_[x]] < tmp) {
                seen_[j] = seen_[x];
                seen_idx_[seen_[j]] = j;
                j = x; x = j*2+1;
            } else break;
        }
        seen_[j] = seen_[top];
        if (seen_idx_[seen_[j]] < n_) seen_idx_[seen_[j]] = j;
    }

    // restore
    for (unsigned int i = 0; i < bottom; i++) {
//        dist_[rst[i].second] = INF;
//        seen_idx_[rst[i].second] = n_;
        dist_[children_[i]] = INF;
        seen_idx_[children_[i]] = n_;
    }

    for (unsigned int i = 0; i < top; i++) {
        seen_idx_[seen_[i]] = n_;
    }

//    dTimer += double(clock() - dStart);

//    return rst;
//    dSize += bottom;
    return bottom;
}

void TIM::LoadTIC(char *filename) {
    TIC::Build2TICFromFile(filename);
}

void TIM::UpdateHeap(bool bounded) {
    if (cursorL_ < n_) {
        double maxL = L_[cursorL_].first;
        double maxH = H_.empty() ? 0 : H_.front().infl;
        while (maxH < maxL && cursorL_ < n_) {
//            printf("maxL: (%d, %lg)\n", L_[cursorL_].id, L_[cursorL_].infl);
//            printf("maxH: (%d, %lg)\n", H_.empty() ? -1 : H_.front().id, maxH);
            di ln = L_[cursorL_++];

            HNode c;
            c.id = ln.second;
            c.infl = ln.first;
            if (bounded) Bounded(c);
            else Exact(c);
            H_.push_back(c); push_heap(H_.begin(), H_.end(), max_candidate_cmp);

            maxL = L_[cursorL_].first;
            maxH = H_.front().infl;
        }

//        printf("|H_| = %d\n", H_.size());
//        printf("maxL: (%d, %lg), maxH: (%d, %lg)\n", 
//            L_[cursorL_].second, L_[cursorL_].first,
//            H_.empty() ? -1 : H_.front().id, maxH);
    }
}

// double TIM::DijkstraV(int u, double max) {
//    vector<di> rst = Dijkstra(u, max);
//    double rt = 0;
//    for (unsigned int i = 0; i < rst.size(); i++) {
//        rt += rst[i].first;
//    }
//
//    return rt;
// }

di TIM::BestFirst() {
    UpdateHeap(true);
    HNode seed = H_.front();
//    tStart = clock();  // not here
    pop_heap(H_.begin(),H_.end(), max_candidate_cmp); H_.pop_back();
//    tTimer += double(clock()-tStart);
//    while (seed.round < round_) {
    while (seed.round < round_ || seed.status != COMPUTED) {
//        printf("check seed (%d, %lg, %d, %d)\n", seed.id, seed.infl, seed.round, seed.status);
        if (seed.round < round_) {  // estimate
            Bounded(seed);
            H_.push_back(seed); push_heap(H_.begin(), H_.end(), max_candidate_cmp);
        } else if (seed.status == ESTIMATED) {  // compute
            Exact(seed);
            H_.push_back(seed); push_heap(H_.begin(), H_.end(), max_candidate_cmp);
        }
//        Exact(seed);
//        tStart = clock();
//        H_.push_back(seed); push_heap(H_.begin(), H_.end(), max_candidate_cmp);
//        tTimer += double(clock()-tStart);

        UpdateHeap(true);
        seed = H_.front();
//        tStart = clock();
        pop_heap(H_.begin(),H_.end(), max_candidate_cmp); H_.pop_back();
//        tTimer += double(clock()-tStart);
    }

//    printf("BestFirst round %d: s = (%d, %lg)\n", round_, seed.id, seed.infl);

    return di(seed.infl, seed.id);
}

di TIM::BestFirstNoBound() {
    UpdateHeap(false);
    HNode seed = H_.front();
//    tStart = clock();  // not here
    pop_heap(H_.begin(),H_.end(), max_candidate_cmp); H_.pop_back();
//    tTimer += double(clock()-tStart);
    while (seed.round < round_) {
//        printf("check seed (%d, %lg, %d, %d)\n", seed.id, seed.infl, seed.round, seed.status);
        Exact(seed);
//        tStart = clock();
        H_.push_back(seed); push_heap(H_.begin(), H_.end(), max_candidate_cmp);
//        tTimer += double(clock()-tStart);

        UpdateHeap(false);
        seed = H_.front();
//        tStart = clock();
        pop_heap(H_.begin(),H_.end(), max_candidate_cmp); H_.pop_back();
//        tTimer += double(clock()-tStart);
    }

//    printf("BestFirst round %d: s = (%d, %lg)\n", round_, seed.id, seed.infl);

    return di(seed.infl, seed.id);
}

void TIM::Exact(HNode& seed) {
    exactCounter += 1;
//    eStart = clock();

//   printf("Exact (%d)\n", u);
    int bottom = Dijkstra(seed.id, -log(theta_));

    // calculate marginal influence of seed
    tap_[seed.id] = 1.0f;
    double rt = tap_[seed.id] - ap_[seed.id], minflu;
//    iees_cache_[seed.id].resize(rst.size());
    iees_cache_[seed.id].resize(bottom);
    iees_cache_[seed.id][0] = di(tap_[seed.id] - ap_[seed.id], seed.id);

    int id;
//    for (unsigned int i = 1; i < rst.size(); i++) {
    for (unsigned int i = 1; i < bottom; i++) {
//        id = rst[i].second;
        id = children_[i];
        minflu = MarginalAPOf(id);
        rt += minflu;
//        printf("ier %d, iee %d, ap: %lg, minflu: %lg, %lg\n", rst[pre[i]].second, id, ap_[id], minflu, rst[i].first);
//        rst[i] = di(minflu, id);
        iees_cache_[seed.id][i] = di(minflu, id);
    }

    for (int i = 0; i < bottom; i++) {
//    for (int i = 0; i < rst.size(); i++) {
        // restore
        tap_[children_[i]] = ap_[children_[i]];
//        tap_[rst[i].second] = ap_[rst[i].second];
    }

//    printf("Exact (%d, %lg, %lg, %d)\n", seed.id, rt, seed.infl, rst.size());

//    assert( rt <= seed.infl );  // marginal influence submodularity.

    seed.infl = rt;
    seed.round = round_;
    seed.status = COMPUTED;

//    exactTimer += double(clock() - eStart);
}

void TIM::Bounded(HNode& seed) {
//    printf("Bounded (%d)\n", u);
    boundCounter += 1;
//    bStart = clock();

//    if (Iee_[seed.id].empty()) {
//        Exact(seed);
//        return;
//    }

    int bottom = Dijkstra(seed.id, -log(theta2_));
    tap_[seed.id] = 1.0f;
    double rt = tap_[seed.id] - ap_[seed.id], minflu;
    int id;
    for (unsigned int i = 1; i < bottom; i++) {
        id = children_[i];
        minflu = MarginalAPOf(id);
        rt += minflu;
    }
    double rt2 = 0; int rt2s = 0;
    for (unsigned int i = 0; i < bottom; i++) {
        id = children_[i];
        if (is_border_[id]) {
            rt2 += (tap_[id]-ap_[id]) * infl2_[id];
            rt2s += 1;
            is_border_[id] = false;  /* restore */
        } 
    }
    /* restore */
    for (int i = 0; i < bottom; i++) tap_[children_[i]] = ap_[children_[i]];
//    printf("Bounded (%d, %lg, %lg, %d)\n", seed.id, rt + rt2, rt, rt2s);
    
//    double infl = 0;
//    vector<di> iees = Iee_[seed.id];
//    double minflu;
//    for (unsigned int i = 0; i < iees.size(); i++) {
//        minflu = (1.0 - ap_[iees[i].second]) * iees[i].first;
//        infl += minflu;
//        printf("iee %d, ap: %lg, minflu: %lg, %lg\n", iees[i].second, ap_[iees[i].second], minflu, iees[i].first);
//    }

//    boundTimer += double(clock() - bStart);

//    printf("Bounded (%d, %lg), (%d, %lg)\n", seed.id, infl, seed.id, seed.infl);
    seed.infl = min(rt + rt2, seed.infl);  /* lazy forward */
    seed.round = round_;
    seed.status = ESTIMATED;
}

// vector<vector<int> > TIM::LBGenerate(vector<int> seeds) {
//     vector<vector<int> > rst;
//     set<int> sset;
//     for (unsigned int i = 0; i < seeds.size(); i++) {
//         vector<di> r = TIC::DijkstraWithStopIds(seeds[i], -log(theta_), sset);
//         sset.insert(seeds[i]);
//         vector<int> t(r.size());
//         for (unsigned int j = 0; j < r.size(); j++) {
//             t[j] = r[j].second;
//         }
//         rst.push_back(t);
//     }
//     return rst;
// }

void TIM::UBLoad(char *filename) {
#ifdef VERBOSE
    printf("Loading Upper Bound samples ...\n");
#endif
    int n_samples, Z, k;
    FILE *f = fopen(filename, "r");
    fscanf(f, "%d\t%d\t%d", &n_samples, &Z, &k);
    assert(Z == TIC::GetZ());
    ubsamples_.resize(n_samples);

    for (int i = 0; i < n_samples; i++) {
        fscanf(f, "%lg", &(ubsamples_[i].ub));
        ubsamples_[i].item.resize(Z);
        for (int z = 0; z < Z; z++) {
            fscanf(f, "%lg", &(ubsamples_[i].item[z]));
        }
    }
    fclose(f);
}

void TIM::LBLoad(char *filename) {
#ifdef VERBOSE
    printf("Loading Lower Bound samples ...\n");
#endif
    int n_samples, Z, k;
    FILE *f = fopen(filename, "r");
    fscanf(f, "%d %d %d", &n_samples, &Z, &k);
    assert(Z == TIC::GetZ());
    lbsamples_.resize(n_samples);

//    vector<vector<int> >lbouts;
//    vector<int> lbout; int size;
    for (int i = 0; i < n_samples; i++) {
//        lbouts.clear();
        lbsamples_[i].item.resize(Z);
        for (int z = 0; z < Z; z++) {
            fscanf(f, "%lg", &(lbsamples_[i].item[z]));
        }
        lbsamples_[i].seeds.resize(k);
        for (int j = 0; j < k; j++) {
            fscanf(f, "%d", &(lbsamples_[i].seeds[j]));
//            lbout.clear();
//            fscanf(f, "%d ", &size);
//            for (int v = 0; v < size; v++) {
//                fscanf(f, "%d", &v);
//                lbout.push_back(v);
//            }
//            lbouts.push_back(lbout);
        }
//        lbsamples_[i].lbouts = lbouts;
    }
    fclose(f);
}

double TIM::LBMarginalAPOf(const int v) {
    if (used_[v]) return 0;
    double ap_v = 1.0;
    for (unsigned int i = 0; i < TIC::GetInNeighbor(v); i++) {
        Edge e = TIC::GetInEdge(v, i);
        if (lap_[e.u] > 0) {
//            printf("e(%d, %d, %lg), ap(%d)=%lg\n", e.u, e.v, TIC::PProb(q_.item, e.w), e.u, tap_[e.u]);
            ap_v *= (1 - TIC::GetInEdgeWeight(v, i, q_.item)*lap_[e.u]);
        }
    }
    ap_v = 1.0 - ap_v;
    lap_[v] = ap_v;

    assert (ap_v >= ap_[v]);
//    printf("ap(%d)=%lg, tap=%lg\n", v, ap_[v], tap_[v]);
    return ap_v - ap_[v];
}

// double TIM::LBCalculate(const vector<vector<int> >& lbouts, vector<di>& rst) {
//     rst.clear();
//     for (unsigned int i = 0; i < n_; i++) lap_[i] = ap_[i];
//     double rt = 0, minflu, m; int s;
//     for (unsigned int i = 0; i < lbouts.size(); i++) {
//         const vector<int>& lbout = lbouts[i];
//         s = lbout[0];
//         minflu = 1.0f - lap_[s];
//         rt += minflu;
//         lap_[s] = 1.0f;  // seed
//         for (unsigned int j = 1; j < lbout.size(); j++) {
//             m = LBMarginalAPOf(lbout[j]);
//             rt += m;
//             minflu += m;
//         }
//         rst.push_back(di(minflu, s));
//     }
// 
//     return rt;
// }

double TIM::LBCalculate(const vector<int>& seeds) {
    for (unsigned int i = 0; i < n_; i++) lap_[i] = ap_[i];
    double rt = 0;
    for (unsigned int i = 0; i < seeds.size(); i++) {
        rt += 1 - ap_[seeds[i]];
        lap_[seeds[i]] = 1.0f;  // seed
        used_[seeds[i]] = true;
    }
    for (unsigned int i = 0; i < n_; i++) {
        LBMarginalAPOf(n_-i-1);
    }
    for (unsigned int i = 0; i < n_; i++) {
        rt += LBMarginalAPOf(i);
    }

    for (unsigned int i = 0; i < seeds.size(); i++) {
        /* restore */
        used_[seeds[i]] = false;
    }

    return rt;
}

double KLD(const vector<double>& P, const vector<double>& Q) {
//    assert(P.size() == Q.size());
//    printf("P:\n");
//    for (unsigned int i = 0; i < P.size(); i++) {
//        printf("%lg ", P[i]);
//    }
//    printf("\nQ:\n");
//    for (unsigned int i = 0; i < Q.size(); i++) {
//        printf("%lg ", Q[i]);
//    }
//    printf("\n");
    double rt = 0;
    double p, q;
    for (unsigned int i = 0; i < P.size(); i++) {
        p = P[i] == 0 ? 0.0000001 : P[i];
	q = Q[i] == 0 ? 0.0000001 : Q[i];
        rt += log(p/q) * p;
    }
//    printf("rt: %lg\n", rt);
    return rt;
}

double TIM::UBCalculate() {
    double rt;

    // KLD
    double m = 1.0f, val; int idx = 0;
    for (unsigned int i = 0; i < ubsamples_.size(); i++) {
        UBSample& ub = ubsamples_[i];
        val = KLD(q_.item, ub.item);
//        printf("KLD %lg\n", val);
        if (val < m) {
            m = val;
            rt = ub.ub;
            idx = i;
        }
    }

#ifdef VERBOSE
    printf("chosen %dth upper bound item sample:\n", idx);
    for (unsigned int i = 0; i < ubsamples_[idx].item.size(); i++) {
        printf("%lg ", ubsamples_[idx].item[i]);
    }
    printf("bound: %lg", ubsamples_[idx].ub);
    printf("\n");
#endif

    // Cosin
    return rt;
}

// vector<vector<int> > TIM::LBSelect() {
//     // KLD
//     double m = 1.0f, val; int key = 0;
//     for (unsigned int i = 0; i < lbsamples_.size(); i++) {
//         LBSample& lb = lbsamples_[i];
// //        printf("lb item %d", lb.item.size());
//         val = KLD(q_.item, lb.item);
//         if (val < m) {
//             m = val;
//             key = i;
//         }
//     }
// 
//     printf("chosen lower bound item sample:\n");
//     for (unsigned int i = 0; i < lbsamples_[key].item.size(); i++) {
//         printf("%lg ", lbsamples_[key].item[i]);
//     }
//     printf("\n");
// 
//     // Cosin
//     return lbsamples_[key].lbouts;
// }

vector<int> TIM::LBSelect() {
    // KLD
    double m = 1.0f, val; int key = 0;
    for (unsigned int i = 0; i < lbsamples_.size(); i++) {
        LBSample& lb = lbsamples_[i];
//        printf("lb item %d", lb.item.size());
        val = KLD(q_.item, lb.item);
        if (val < m) {
            m = val;
            key = i;
        }
    }

#ifdef VERBOSE
    printf("chosen lower bound item sample:\n");
    for (unsigned int i = 0; i < lbsamples_[key].item.size(); i++) {
        printf("%lg ", lbsamples_[key].item[i]);
    }
    printf("\n");
#endif

    // Cosin
    return lbsamples_[key].seeds;
}
