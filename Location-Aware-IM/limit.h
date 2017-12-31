#ifndef LIMIT_H
#define LIMIT_H

#define MAX_NODE	40000
// #define	MAX_EDGE	500000
#define MAX_K		5000
// // #define MAX_POPULATION  100000
// #define MAX_REGION	100000
// #define POPULATION_COVERAGE 0.9
// 
// #define PROTOTYPE true
#define PROTOTYPE false
// 
#define STR_LEN		200
// 
// #define NUM_ITER	200
#define NUM_ITER	10
// #define SET_SIZE	50
#define EPS 1e-10

#include <vector>

using namespace std;

class Args
{
    private:
        // if this is a native
        static vector<bool> is_native;
        // if this is a candidate
        static vector<bool> is_cdd;

    public:
        static int MAX_QK;  // 5000
        static int MAX_CK;  // 5000
        static int MIN_CK;  // 10
        static int CACHE_RATIO;  // 10
        static int CAPACITY;    // 500
        static void reset(int n)
        {
            is_native.resize(n);
            is_cdd.resize(n);
            for (int _ = 0; _ < n; _++) {
                is_native[_] = false;
                is_cdd[_] = true;
            }
        }

        static void setNative(int i)
        {
            is_native[i] = true;
        }
        
        static void setCandidate(int i)
        {
            is_cdd[i] = true;
        }
        
        static bool isCandidate(int i)
        {
            return is_cdd[i];
        }
        
        static bool isNative(int i)
        {
            return is_native[i];
        }
};

#endif
