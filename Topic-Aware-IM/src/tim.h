#ifndef TIM_H
#define TIM_H

#include <set>
#include <vector>

#include "tic.h"

using namespace std;


typedef struct
{
    vector<double> item;
    int k;
} Q;

typedef struct
{
    vector<int> ids;
    vector<double> minflu;
} S;

// typedef struct
// {
//     vector<di> iees;
//     vector<int> pre;
// } Out;  // out arborescence

#define ESTIMATED 0
#define COMPUTED 1
typedef struct
{
    int id;
    double infl;
    int round;
    int status;  // ESTIMATED, COMPUTED
} HNode;

// typedef struct
// {
//     int id;
//     double infl;
//     double infl2;  // theta / theta2_
// } LNode;

// Lower Bound
typedef struct
{
//    vector<vector<int> > lbouts;
    vector<int> seeds;
    vector<double> item;
} LBSample;

typedef struct
{
    vector<double> item;
    double ub;
} UBSample;

class TIM 
{
private:
    // needed for all methods
    static int n_;
    static double *ap_;
    static Q q_;
    static int round_;
    static bool *used_;

    // not needed for every method
    // Greedy
    static vector<HNode> H_;
    static vector<vector<di> > iees_cache_;
    static double* tap_;  // temporary ap
//    static vector<int> pre_;
//    static vector<set<int> > In_;
    // List
    static vector<di> L_;  // \List  <influ, number_of_iee>
    static double *infl2_;  /* sqrt(theta_) */
    static int cursorL_;  // cursor on \List
    // Dijkstra
    static double *dist_;
    static int *seen_;
    static int *seen_idx_;
    static int *children_;
    static int *pred_;  /* predecessor id */
    static bool *is_border_;  /* if on the border, upper bounded needed */

//    static vector<vector<di> > Iee_;  // exact calculated out arborescence
    // Lower
    static double *lap_;
    static vector<LBSample> lbsamples_;
    // upper
    static vector<UBSample> ubsamples_;

    static int Dijkstra(int, double);
//    static double DijkstraV(int, double);
    static void UpdateHeap(bool);
    static void Exact(HNode&);
    static void Bounded(HNode&);
    static di BestFirst(); 
    static di BestFirstNoBound(); 
    static void UpdateAP(const int);
    static double MarginalAPOf(const int);

    static void Init();

    // Lower
    static void LBLoad(char *);
//    static double LBCalculate(const vector<vector<int> >&, vector<di>&);
    static double LBCalculate(const vector<int>&);
    static vector<int> LBSelect();
    static double LBMarginalAPOf(const int);

    // Upper Approx
    static void UBLoad(char *);
    static double UBCalculate();

    static S GreedyBound(Q, bool);

public:
    // load
    static void LoadTIC(char *);
    static double theta_;
    static double theta2_;

    // query
    static void InitGreedy(char *);
    static S Greedy(Q);
    static S GreedyNoBound(Q);
//    static void InitGreedyUpApprox();
//    static S GreedyUpApprox(Q);
    static void InitApprox(char *, char *, char *);
    static S Approx(Q, double epsilon);

//    static vector<vector<int> > LBGenerate(vector<int>);

    static void Reset();
    static void PrintArguments();
};

#endif
