#include<iostream>
#include<set>

#include "limit.h"
#include "graph.h"
#include "Quadtree.h"
#include "pmia_upgrade.h"
#include "general_cascade.h"
#include "Cache.h"
#include "mip.h"

using namespace std;

clock_t start;

void readLocations(Quadtree *qt, char* arg) {
    // read loc
    printf("Quadtree insertion begin...\n");
    int numLoc = 0;
    char filename_loc[]="very_very_long_loc_file_name.loc";
    sprintf(filename_loc, "data/%s.loc",  arg);
    FILE *loc = fopen(filename_loc, "r");
    fscanf(loc, "%d", &numLoc);
    int id;
    float x, y; 
    for (int _ = 0; _ < numLoc; _++) {
        fscanf(loc, "%d %f %f", &id, &y, &x);
        if (!qt->insert(id, x, y)) {
            printf("Quadtree insertion failed: %d %f %f\n", id, y, x);
        }
    }
    fclose(loc);

    printf("Quadtree insertion done. %d population, %d regions total.\n", qt->population, qt->count);
}

void readRelations(char* arg) {
    // read influence
    char filename_inf[]="very_very_long_inf_file_name.inf";
    sprintf(filename_inf, "data/%s.inf",  arg);
    printf("Build WC model...\n");
    Graph::Build2WCFromFile(filename_inf);
    GeneralCascade::Build();
    cout << "Build WC model done."<<endl;
}

double toSimulateOnce(int setsize, int (*GetNode)(int i), double (*Run)(int num_iter, int size, int set[]))  // simulate influence spread at one step
{
    int set[MAX_NODE];
    int t;
    for (t=0; t<setsize; t++)
    {
        set[t] = GetNode(t);
    }
    return Run(NUM_ITER, t, set);
}

double toSimulateOnce2(int setsize, vector<int>topk, double (*Run)(int num_iter, int size, int set[]))  // simulate influence spread at one step
{
    int set[MAX_NODE];
    int t;
    for (t=0; t<topk.size() && t<setsize; t++)
    {
        set[t] = topk[t];
    }
    return Run(NUM_ITER, t, set);
}

bool baseline(Quadtree *qt, AABB boundary, int qK, set<int> &dTopkSet, int
        *numNatives, double *dTimer, double *dSpread)
{
    start = clock();

    vector<int> natives; 
    qt->queryRange(natives, boundary);
    printf("%d natives in this query.\n", int(natives.size()));
    *numNatives = int(natives.size());

    Args::reset(Graph::GetN());
    // set new natives and candidates
    for (int i = 0; i < natives.size(); i++) {
        Args::setNative(natives[i]);
    }

    if (qK * Args::CACHE_RATIO > natives.size()) {
        printf("k == %d, only %d natives, abort this query.\n", 
                qK, int(natives.size()));
        return false;
    }

    // reset natives and candidates
    SPT_upgrade::reset_native();
    for (int i = 0; i < natives.size(); i++) {
        SPT_upgrade::is_native[natives[i]] = true;
    }

    *dTimer = (double)(clock() - start) / CLOCKS_PER_SEC;

    SPT_upgrade::reset();

    start = clock();
    SPT_upgrade::select(qK);
    *dTimer += (double)(clock() - start) / CLOCKS_PER_SEC;
    *dSpread=toSimulateOnce(
            qK, SPT_upgrade::GetNode, GeneralCascade::Run);
    for (int i = 0; i < qK; i++) {
        dTopkSet.insert(SPT_upgrade::GetNode(i));
    }
    return true;
}

bool expansion(Quadtree* qt, AABB boundary, int qK, vector<int> &eTopk, int
        *numNatives, double *eTimer, double *eSpread)
{
    start = clock();

    double timerNatives = 0.0;
    vector<int> natives; 
    qt->queryRange(natives, boundary);
//    printf("%d natives in this query.\n", int(natives.size()));
    *numNatives = int(natives.size());
    timerNatives = (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("expansion query find natives operation %f seconds.\n", timerNatives);

    Args::reset(Graph::GetN());
    // set new natives and candidates
    for (int i = 0; i < natives.size(); i++) {
        Args::setNative(natives[i]);
    }

    if (qK * Args::CACHE_RATIO > natives.size()) {
//        printf("k == %d, only %d natives, abort this query.\n", 
//                qK, int(natives.size()));
        return false;
    }

    // reset natives and candidates
    Mip::reset_native();
    for (int i = 0; i < natives.size(); i++) {
        Mip::is_native[natives[i]] = true;
    }
    // select
    eTopk = Mip::select(qK);

    *eTimer = (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("expansion query compute %d times.\n", Mip::count);
    printf("expansion query estimate %d times.\n", Mip::estCount);
    printf("expansion query compute %f seconds.\n", Mip::incInflTimer);

    Mip::reset();

    *eSpread=toSimulateOnce2(qK, eTopk, GeneralCascade::Run);
    return true;
}

void sortInitOneRegion(Quadtree *qt, int rid, AABB boundary)
{
    vector<int> list; 
    qt->queryRange(list, boundary);
    for (int i = 0; i < list.size(); i++) {
        printf("%d ", list[i]);
    }
    printf("\n");
    vector<int> cdds;
    vector<double> cddsDp; 
    Mip::sortInit(list, cdds, cddsDp);
    Cache::setCandidates(rid, cdds);
    Cache::setDP(rid, cddsDp);
}

void sortInitAllRegion(Quadtree *qt, Quadtree* lqt)
{
    printf("XY:(%f %f), halfDimension:(%f %f)\n", (lqt->boundary).center.x, (lqt->boundary).center.y, (lqt->boundary).halfDimension.x, (lqt->boundary).halfDimension.y);
    sortInitOneRegion(qt, lqt->id, lqt->boundary);
    if (lqt->NW == NULL) return;
    sortInitAllRegion(qt, lqt->NW);
    sortInitAllRegion(qt, lqt->NE);
    sortInitAllRegion(qt, lqt->SW);
    sortInitAllRegion(qt, lqt->SE);
}

bool assembly(Quadtree* qt, AABB boundary, int qK, vector<int> &aTopk, int
        *numNatives, double *aTimer, double *aSpread) 
{
    start = clock();

    double timerNatives = 0.0;
    double timerAssembly = 0.0;
    double timerLists = 0.0;

    // find natives
    vector<int> rids;
    vector<AABB> boundaries;
    qt->queryRangeRegion(rids, boundaries, boundary);

    vector<int> natives; 
    vector<int> nativesNotCovered;
    qt->queryRangeFO(natives, nativesNotCovered, boundaries, boundary);
//    printf("%d natives in this query.\n", int(natives.size()));
    *numNatives = int(natives.size());
    timerNatives += (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("assembly query find natives operation %f seconds.\n", timerNatives);

    Args::reset(Graph::GetN());
    // set new natives and candidates
    for (int i = 0; i < natives.size(); i++) {
        Args::setNative(natives[i]);
    }
    if (qK * Args::CACHE_RATIO > natives.size()) {
//        printf("k == %d, only %d natives, abort this query.\n", 
//                qK, int(natives.size()));
        return false;
    }

    // if region not cached, cache it.
    for (int i = 0; i < rids.size(); i++) {
        if (!Cache::isInitCached(rids[i])) {
            printf("Opps, need to cache region #%d first...\n", rids[i]);
            sortInitOneRegion(qt, rids[i], boundaries[i]);
        }
    }

    start = clock();
    int numOfRegions = rids.size();
    vector<vector<int> > candidates;
    vector<vector<double> > initDp;
    candidates.resize(numOfRegions + 1);
    initDp.resize(numOfRegions + 1);
    // collect all candidates, and their initial dp
    for (int i = 0; i < numOfRegions; i++) {
        // collect all possible candidates
        Cache::getCandidates(candidates[i], rids[i], -1); 
        Cache::getDP(initDp[i], rids[i], -1);
    }
    // fringe other users
    Mip::sortInit(nativesNotCovered, candidates[numOfRegions], initDp[numOfRegions]);
    timerLists = (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("assembly query construct lists operation %f seconds.\n", timerLists);

    start = clock();
    // reset natives and candidates
    Mip::reset_native();
    for (int i = 0; i < natives.size(); i++) {
        Mip::is_native[natives[i]] = true;
    }
    aTopk = Mip::assembly(qK, candidates, initDp);
    *aTimer = (double)(clock() - start) / CLOCKS_PER_SEC;
    *aTimer += timerNatives;
    *aTimer += timerLists;

    printf("assembly query compute %d times.\n", Mip::count);
    printf("assembly query estimate %d times.\n", Mip::estCount);
    printf("assembly query compute %f seconds.\n", Mip::incInflTimer);

    Mip::reset();

    *aSpread=toSimulateOnce2(qK, aTopk, GeneralCascade::Run);
    return true;
}

void findTauOneRegion(Quadtree *qt, int rid, AABB boundary, int k)
{
    vector<int> natives; 
    qt->queryRange(natives, boundary);
    vector<int> topTau;
    vector<double> tauDp; 
    // reset natives and candidates
    Mip::reset_native();
    for (int i = 0; i < natives.size(); i++) {
        Mip::is_native[natives[i]] = true;
    }
    k = k > natives.size()/Args::CACHE_RATIO ? natives.size()/Args::CACHE_RATIO : k;

    Mip::cacheTau(k, topTau, tauDp);
    Mip::reset();

    Cache::setTopTau(rid, topTau);
    Cache::setDPTau(rid, tauDp);
}

bool hint(Quadtree* qt, AABB boundary, int qK, double acceptRatio, vector<int> &hTopk, int *numNatives, double *hTimer, double *hSpread) {
    start = clock();

    double timerNatives = 0.0;
    double timerHint = 0.0;
    double timerLists = 0.0;

    // find natives
    vector<int> rids;
    vector<AABB> boundaries;
    qt->queryRangeRegion(rids, boundaries, boundary);

    vector<int> natives; 
    vector<int> nativesNotCovered;
    qt->queryRangeFO(natives, nativesNotCovered, boundaries, boundary);
//    printf("%d natives in this query.\n", int(natives.size()));
    *numNatives = int(natives.size());
    timerNatives += (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("hint query find natives operation %f seconds.\n", timerNatives);

    Args::reset(Graph::GetN());
    // set new natives and candidates
    for (int i = 0; i < natives.size(); i++) {
        Args::setNative(natives[i]);
    }
    if (qK * Args::CACHE_RATIO > natives.size()) {
//        printf("k == %d, only %d natives, abort this query.\n", 
//                qK, int(natives.size()));
        return false;
    }

    // if region not cached, cache it.
    for (int i = 0; i < rids.size(); i++) {
        if (!Cache::isInitCached(rids[i])) {
//            printf("Opps, need to cache region #%d first...\n", rids[i]);
            sortInitOneRegion(qt, rids[i], boundaries[i]);
        }
    }
    for (int i = 0; i < rids.size(); i++) {
        if (!Cache::isTauCached(rids[i])) {
//            printf("Opps, need to cache tau region #%d first...\n", rids[i]);
            findTauOneRegion(qt, rids[i], boundaries[i], qK);
        }
    }

    start = clock();
    int numOfRegions = rids.size();
    vector<vector<int> > candidates(numOfRegions + 1);
    vector<vector<double> > initDp(numOfRegions + 1);
    vector<vector<int> > topTau(numOfRegions);
    vector<vector<double> > dpTau(numOfRegions);
    // collect all candidates, and their initial dp
    for (int i = 0; i < numOfRegions; i++) {
        // collect all possible candidates
        Cache::getCandidates(candidates[i], rids[i], -1); 
        Cache::getDP(initDp[i], rids[i], -1);
        Cache::getTopTau(topTau[i], rids[i], qK); 
        Cache::getDPTau(dpTau[i], rids[i], qK);
    }
    Mip::sortInit(nativesNotCovered, candidates[numOfRegions], initDp[numOfRegions]);
    timerLists = (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("hint query construct lists operation %f seconds.\n", timerLists);
    // fringe other users

    start = clock();
    // reset natives and candidates
    Mip::reset_native();
    for (int i = 0; i < natives.size(); i++) {
        Mip::is_native[natives[i]] = true;
    }
    hTopk = Mip::hint(qK, acceptRatio, candidates, initDp, topTau, dpTau);

    *hTimer = (double)(clock() - start) / CLOCKS_PER_SEC;
    *hTimer += timerNatives;
    *hTimer += timerLists;

    printf("hint query compute %d times.\n", Mip::count);
    printf("hint query estimate %d times.\n", Mip::estCount);
    printf("hint query compute %f seconds.\n", Mip::incInflTimer);

    Mip::reset();

    *hSpread=toSimulateOnce2(qK, hTopk, GeneralCascade::Run);
    return true;
}

bool dumbHint(Quadtree* qt, AABB boundary, int qK, double acceptRatio, vector<int> &hTopk, int *numNatives, double *hTimer, double *hSpread) {
    start = clock();

    double timerNatives = 0.0;
    double timerHint = 0.0;
    double timerLists = 0.0;

    // find natives
    vector<int> rids;
    vector<AABB> boundaries;
    qt->queryRangeRegion(rids, boundaries, boundary);

    vector<int> natives; 
    vector<int> nativesNotCovered;
    qt->queryRangeFO(natives, nativesNotCovered, boundaries, boundary);
//    printf("%d natives in this query.\n", int(natives.size()));
    *numNatives = int(natives.size());
    timerNatives += (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("dumb hint query find natives operation %f seconds.\n", timerNatives);

    Args::reset(Graph::GetN());
    // set new natives and candidates
    for (int i = 0; i < natives.size(); i++) {
        Args::setNative(natives[i]);
    }
    if (qK * Args::CACHE_RATIO > natives.size()) {
//        printf("k == %d, only %d natives, abort this query.\n", 
//                qK, int(natives.size()));
        return false;
    }

    // if region not cached, cache it.
    for (int i = 0; i < rids.size(); i++) {
        if (!Cache::isInitCached(rids[i])) {
//            printf("Opps, need to cache region #%d first...\n", rids[i]);
            sortInitOneRegion(qt, rids[i], boundaries[i]);
        }
    }

    start = clock();
    int numOfRegions = rids.size();
    vector<vector<int> > candidates(numOfRegions + 1);
    vector<vector<double> > initDp(numOfRegions + 1);
    // collect all candidates, and their initial dp
    for (int i = 0; i < numOfRegions; i++) {
        // collect all possible candidates
        Cache::getCandidates(candidates[i], rids[i], -1); 
        Cache::getDP(initDp[i], rids[i], -1);
    }
    // fringe other users
    Mip::sortInit(nativesNotCovered, candidates[numOfRegions], initDp[numOfRegions]);
    timerLists = (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("dumb hint query construct lists operation %f seconds.\n", timerLists);

    start = clock();
    // reset natives and candidates
    Mip::reset_native();
    for (int i = 0; i < natives.size(); i++) {
        Mip::is_native[natives[i]] = true;
    }
    hTopk = Mip::dumbHint(qK, acceptRatio, candidates, initDp);

    *hTimer = (double)(clock() - start) / CLOCKS_PER_SEC;
    *hTimer += timerNatives;
    *hTimer += timerLists;

    printf("dumb hint query compute %d times.\n", Mip::count);
    printf("dumb hint query estimate %d times.\n", Mip::estCount);
    printf("dumb hint query compute %f seconds.\n", Mip::incInflTimer);

    Mip::reset();

    *hSpread=toSimulateOnce2(qK, hTopk, GeneralCascade::Run);
    return true;
}

void cacheOneRegion(Quadtree *qt, int rid, AABB boundary)
{
//    Cache::setIsCached(rid, true);
//    SPT_cache::reset_native();
//    vector<int> list; 
//    vector<float> nativesX;
//    vector<float> nativesY;
//    qt->queryRange(list, nativesX, nativesY, boundary);
//    for (int j = 0; j < list.size(); j++) {
//        SPT_cache::is_native[list[j]] = true;
//    }
////    int ck = list.size() / Args::CACHE_RATIO;
////    ck = ck > Args::MAX_CK ? Args::MAX_CK : ck;
////    ck = ck < Args::MIN_CK ? Args::MIN_CK : ck;
//    int ck = list.size();
//    ck = ck > Args::MAX_CK ? Args::MAX_CK : ck;
//    vector<int> cdds; 
//    vector<double> cddsDp; 
//    vector<int> topTau;
//    vector<double> dpTau;
//    SPT_cache::cache(ck, bound, true, cdds, cddsDp, topTau, dpTau);
//    printf("still have noncached region! rid == %d, population == %d, caching %d candidates \n", rid, int(list.size()), int(topTau.size()));
//    Cache::setCandidates(rid, cdds);
//    Cache::setDP(rid, cddsDp);
//    Cache::setTopTau(rid, topTau);
//    Cache::setDPTau(rid, dpTau);
}

int main(int argc, char * argv[])
{
    system("cd tmp");

    if (argc<=1) { 
        printf("Apply PMIA algorithm on local topk query problem");
        return 0;
    }
    // set program settings
    string s;
    s = "prototype";
    Quadtree *qt;
    if (!s.compare(argv[3])) {
        Args::MAX_QK = 20;
        Args::MAX_CK = 20;
        Args::MIN_CK = 0;
        Args::CACHE_RATIO = 1;
        Args::CAPACITY = 4;
        qt = new Quadtree(XY(0.0, 0.0), XY(1.0, 1.0), Args::CAPACITY);
    }
    else {
        Args::MAX_QK = 5000;
        Args::MAX_CK = 5000;
        Args::MIN_CK = 10;
        Args::CACHE_RATIO = 100;
        Args::CAPACITY = 500;
        qt = new Quadtree(XY(0.0, 0.0), XY(180.0, 90.0), Args::CAPACITY);
    }

    readLocations(qt, argv[3]);

    readRelations(argv[3]);

//    qt->print();

    // global initialization
    int bound;
    sscanf(argv[2], "%d", &bound);

    // load queries
    int numQueries;
    scanf("%d", &numQueries);
    float xlim, ylim;
    float x, y; 
    int qK;
    double accept;

    printf("PMIA upgrade initialization...\n");

    s="-m";
    if (!s.compare(argv[1])) {
        // record topk result
        char filename_topk[]="very_very_long_topk_file_name.topk";
        sprintf(filename_topk, "tmp/%s_%d_exp", argv[3], bound);
        FILE *out = fopen(filename_topk, "w");

        sprintf(filename_topk, "tmp/%s_%d_base", argv[3], bound);
        FILE *bout = fopen(filename_topk, "w");

        SPT_upgrade::init(bound);
        Mip::init(-log(1.0/double(bound)));
        for (int _ = 0; _ < numQueries; _++) {
            // get query region natives
            scanf("%f %f %f %f %d", &y, &x, &ylim, &xlim, &qK);

            if (qK > Args::MAX_QK) {
                printf("can not answer more than %d topks !\n", Args::MAX_QK);
                continue;
            }
            int numNatives;

            bool succ;
            printf("baseline query begin...\n");
            double dTimer, dSpread;
            set<int> dTopkSet;
            succ = baseline(qt, AABB(XY(x,y), XY(xlim,ylim)), qK, dTopkSet, &numNatives, &dTimer, &dSpread);
            if (!succ) continue;
            printf("baseline query end.\n");

            printf("expansion based query begin...\n");
            double eTimer, eSpread;
            vector<int> eTopk;
            succ = expansion(qt, AABB(XY(x,y), XY(xlim,ylim)), qK, eTopk, &numNatives, &eTimer, &eSpread);
            if (!succ) continue;
            printf("expansion query end.\n");

            fprintf(out, "%lg %d %lg %d ",
                    eTimer, int(eTopk.size()), eSpread, numNatives);
            fprintf(bout, "%lg %d %lg %d ",
                    dTimer, int(dTopkSet.size()), dSpread, numNatives);
            printf("%lg %lg %d %d %lg %lg %lg %d ",
                    eTimer, dTimer,
                    int(eTopk.size()), qK,
                    eSpread, dSpread, eSpread/dSpread,
                    numNatives
                    );
            // accuracy
            int count = 0;
            for (int i = 0; i < eTopk.size(); i++) {
                if (dTopkSet.count(eTopk[i]) != 0)
                    count += 1;
            }
            double accuracy = double(count) / eTopk.size();
            printf("%f ", accuracy);
            printf("\n");
            fprintf(out, "%f ", accuracy);
            fprintf(out, "\n");
            fprintf(bout, "\n");
            fflush(stdout);
            fflush(out);
            fflush(bout);
        }
        Mip::exit();
        SPT_upgrade::exit();
        fclose(out);
        fclose(bout);
    }

    s="-a";
    if (!s.compare(argv[1])) {
        // record topk result
        char filename_topk[]="very_very_long_topk_file_name.topk";
        sprintf(filename_topk, "tmp/%s_%d_ta", argv[3], bound);
        FILE *out = fopen(filename_topk, "w");

        Mip::init(-log(1.0/double(bound)));
        Cache::init(Quadtree::count);

        sortInitAllRegion(qt, qt);

        for (int _ = 0; _ < numQueries; _++) {
            // get query region natives
            scanf("%f %f %f %f %d", &y, &x, &ylim, &xlim, &qK);

            if (qK > Args::MAX_QK) {
                printf("can not answer more than %d topks !\n", Args::MAX_QK);
                continue;
            }
            int numNatives;

            double aTimer, aSpread;
            vector<int> aTopk;
            bool succ;
            printf("assembly based query begin...\n");
            succ = assembly(qt, AABB(XY(x,y), XY(xlim,ylim)), qK, aTopk, &numNatives, &aTimer, &aSpread);
            if (!succ) continue;
            printf("assembly query end.\n");

//            printf("expansion based query begin...\n");
//            double eTimer, eSpread;
//            vector<int> eTopk;
//            succ = expansion(qt, AABB(XY(x,y), XY(xlim,ylim)), qK, eTopk, &numNatives, &eTimer, &eSpread);
//            if (!succ) continue;
//            printf("expansion query end.\n");
//            printf("expansion query compute %d times.\n", Mip::count);
//            printf("expansion query compute %f seconds.\n", Mip::incInflTimer);

            fprintf(out, "%lg %d %lg %d ",
                    aTimer, int(aTopk.size()), aSpread, numNatives);
            printf("%lg %d %lg %d ",
                    aTimer, int(aTopk.size()), aSpread, numNatives);
//            printf("%lg %lg %d %d %lg %lg %lg %d ",
//                    aTimer, eTimer,
//                    int(aTopk.size()), int(aTopk.size()),
//                    aSpread, eSpread, aSpread/eSpread,
//                    numNatives
//                    );
            fprintf(out, "\n");
            printf("\n");
            fflush(stdout);
            fflush(out);
        }
        Mip::exit();
        fclose(out);
    }

    s="-h";
    if (!s.compare(argv[1])) {
        // record topk result
        char filename_topk[]="very_very_long_topk_file_name.topk";
        sprintf(filename_topk, "tmp/%s_%d_hint", argv[3], bound);
        FILE *out = fopen(filename_topk, "w");

        sprintf(filename_topk, "tmp/%s_%d_dumb_hint", argv[3], bound);
        FILE *dhout = fopen(filename_topk, "w");


        Mip::init(-log(1.0/double(bound)));
        Cache::init(Quadtree::count);
        for (int _ = 0; _ < numQueries; _++) {
            // get query region natives
            scanf("%f %f %f %f %d", &y, &x, &ylim, &xlim, &qK);

            if (qK > Args::MAX_QK) {
                printf("can not answer more than %d topks !\n", Args::MAX_QK);
                continue;
            }
            int numNatives;

            bool succ;
//            printf("expansion based query begin...\n");
//            double eTimer, eSpread;
//            vector<int> eTopk;
//            succ = expansion(qt, AABB(XY(x,y), XY(xlim,ylim)), qK, eTopk, &numNatives, &eTimer, &eSpread);
//            if (!succ) continue;
//            printf("expansion query end.\n");
//            printf("expansion query compute %d times.\n", Mip::count);
//            printf("expansion query compute %f seconds.\n", Mip::incInflTimer);

            printf("hint based query begin...\n");
            double acceptRatio = 0.9;
            for (int i = 3; i <= 5; i++) {
                double hTimer, hSpread;
                vector<int> hTopk;
                succ = hint(qt, AABB(XY(x,y), XY(xlim,ylim)), qK, acceptRatio, hTopk, &numNatives, &hTimer, &hSpread);
                if (!succ) continue;

                double dhTimer, dhSpread;
                vector<int> dhTopk;
                succ = dumbHint(qt, AABB(XY(x,y), XY(xlim,ylim)), qK, acceptRatio, dhTopk, &numNatives, &dhTimer, &dhSpread);
                if (!succ) continue;

                fprintf(out, "%lg %d %lg %d ",
                        hTimer, int(hTopk.size()), hSpread, numNatives);
                fprintf(dhout, "%lg %d %lg %d ",
                        dhTimer, int(dhTopk.size()), dhSpread, numNatives);
                printf("%lg %d %lg %d ",
                        hTimer, int(hTopk.size()), hSpread, numNatives);
                printf("%lg %d %lg %d ",
                        dhTimer, int(dhTopk.size()), dhSpread, numNatives);
//                printf("%lg %lg %lg %d %d %d %lg %lg %lg %lg %lg %d ",
//                        hTimer, dhTimer, eTimer, int(hTopk.size()),
//                        int(dhTopk.size()), int(hTopk.size()), hSpread,
//                        dhSpread, eSpread, hSpread/eSpread, dhSpread/eSpread,
//                        numNatives
//                        );
                printf("\n");
                fprintf(out, "\n");
                fprintf(dhout, "\n");
                fflush(stdout);
                fflush(out);
                fflush(dhout);

                acceptRatio -= 0.05;
            }
            printf("hint query end.\n");
        }
        Mip::exit();
        fclose(out);
        fclose(dhout);
    }

    qt->clear();
}
