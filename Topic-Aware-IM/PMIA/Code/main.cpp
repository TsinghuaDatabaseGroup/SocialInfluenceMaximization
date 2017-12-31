#include "limit.h"
#include "graph.h"

#include "random.h"
#include "degree.h"
#include "greedy.h"
#include "degreediscount_ic.h"
#include "weighteddegree.h"

#include "SPM_gc.h"
#include "SP1M_gc.h"
#include "pmia.h"
#include "pagerank.h"
#include "general_cascade.h"
#include "mia.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>

using namespace std;

FILE* timetmpfile;
double timer;
clock_t start, ended;

void toSimulate(char *file, int (*GetNode)(int i), double (*Run)(int num_iter, int size, int set[]))
{
    FILE *out = fopen(file, "w");
    int set[MAX_NODE];
    for (int t=0; t<SET_SIZE; t++)
    {
        set[t] = GetNode(t);
        fprintf(out, "%02d \t %g\n", t+1, Run(NUM_ITER, t+1, set));
    }
    fclose(out);
}

double toSimulateOnce(int setsize, int (*GetNode)(int i), double (*Run)(int num_iter, int size, int set[]))
{
    int set[MAX_NODE];
    int t;
    for (t=0; t<setsize; t++)
    {
        set[t] = GetNode(t);
    }
    return Run(NUM_ITER, t, set);
}


int main(int argc, char * argv[])
{
    srand(time_t(NULL));
    //srand(100);
    //system("mkdir tmp");
    system("cd tmp");

    if (argc<=1) { 
        printf("-st : statistics of the weighted cascade graph\n");
        printf("-stg : statistics of the general independedent cascade graph\n");
        printf("-b : baseline(random, degree, degreediscount, weighteddegree, pagerank) for general ic\n");
        printf("-g : greedy, SPM and SP1M for general ic\n");
        printf("-p bound1 bound2 : PMIA with 1/theta from bound1 to bound 2\n");
        printf("-m bound1 bound2 : MIA with 1/theta from bound1 to bound 2\n");
        printf("example: max_influence -p 20 2000 < hep.inf \n");
        return 0;
    }
    string s;
    s="-st";
    if (!s.compare(argv[1])) {
        Graph::Build();
        Graph::Stats();
    }
    s="-stg";
    if (!s.compare(argv[1])) {
        Graph::Build2GC();
        Graph::Stats();
    }
    s="-b";
    if (!s.compare(argv[1])) {
        Graph::Build2WC();
        GeneralCascade::Build();

        // Random
        start = clock();
        Random::Build();
        ended = clock();
        timer = (double)(ended - start);
        timetmpfile = fopen("time_random.txt", "w");
        fprintf(timetmpfile,"%lg\n", timer);
        fclose(timetmpfile);
        toSimulate("GC_Random.txt", Random::GetNode, GeneralCascade::Run);
        // Weighted Degree
        start = clock();
        WeightedDegree::Build();
        ended = clock();
        timer = (double)(ended - start);
        timetmpfile = fopen("time_weighteddegree.txt", "w");
        fprintf(timetmpfile,"%lg\n", timer);
        fclose(timetmpfile);
        toSimulate("GC_WeightedDegree.txt", WeightedDegree::GetNode, GeneralCascade::Run);

        // Highest Degree
        start = clock();
        Degree::Build();
        ended = clock();
        timer = (double)(ended - start);
        timetmpfile = fopen("time_degree.txt", "w");
        fprintf(timetmpfile,"%lg\n", timer);
        fclose(timetmpfile);
        toSimulate("GC_Degree.txt", Degree::GetNode, GeneralCascade::Run);

        // DegreeDiscountIC
        start = clock();
        DegreeDiscount_IC::Build(0.01);
        ended = clock();
        timer = (double)(ended - start);
        timetmpfile = fopen("time_degreediscount_ic.txt", "w");
        fprintf(timetmpfile,"%lg\n", timer);
        fclose(timetmpfile);
        toSimulate("GC_DiscountIC.txt", DegreeDiscount_IC::GetNode, GeneralCascade::Run);

        // Pagerank
        GeneralCascade::Build();
        start = clock();
        pagerank::Build(SET_SIZE);
        ended = clock();
        timer = (double)(ended - start);
        timetmpfile = fopen("time_pagerank.txt", "w");
        fprintf(timetmpfile,"%lg\n", timer);
        fclose(timetmpfile);
        toSimulate("GC_pagerank.txt", pagerank::GetNode, GeneralCascade::Run);
    }

    // GreedyGC (improved by CELF)
    s="-g";
    if (!s.compare(argv[1])) {
        Graph::Build2WC();
        GeneralCascade::Build();	
        start = clock();
        Greedy::Build(SET_SIZE,GeneralCascade::Run);
        ended = clock();
        timer = (double)(ended - start);
        timetmpfile = fopen("time_greedy_gc.txt", "w");
        fprintf(timetmpfile,"%lg\n", timer);
        fclose(timetmpfile);
        system("copy greedy.txt greedy_gc.txt");
        system("del /Q tmp\\*");
        toSimulate("GC_Greedy.txt", Greedy::GetNode, GeneralCascade::Run);

        // GreedyGC_SPM (improved by CELF)
        start = clock();
        Greedy::Build(SET_SIZE,SPM_gc::Run);
        ended = clock();
        timer = (double)(ended - start);
        timetmpfile = fopen("time_greedy_gc_spm.txt", "w");
        fprintf(timetmpfile,"%lg\n", timer);
        fclose(timetmpfile);
        system("copy greedy.txt greedy_gc_spm.txt");
        system("del /Q tmp\\*");
        toSimulate("GC_SPM.txt", Greedy::GetNode, GeneralCascade::Run);

        // GreedyGC_SP1M (improved by CELF)
        start = clock();
        Greedy::Build(SET_SIZE,SP1M_gc::Run);
        ended = clock();
        timer = (double)(ended - start);
        timetmpfile = fopen("time_greedy_gc_sp1m.txt", "w");
        fprintf(timetmpfile,"%lg\n", timer);
        fclose(timetmpfile);
        system("copy greedy.txt greedy_gc_sp1m.txt");
        system("del /Q tmp\\*");
        toSimulate("GC_SP1M.txt", Greedy::GetNode, GeneralCascade::Run);
    }

    //control bound to test PMIA_GC
    s="-p";
    if (!s.compare(argv[1])) {
        Graph::Build2WC();
        GeneralCascade::Build();
        int bound1=10, bound2=2000;
        if (argc>=3) sscanf(argv[2],"%d",&bound1);
        if (argc>=4) sscanf(argv[3],"%d",&bound2);
        char SPTfilename[]="PMIA_control.txt";
        FILE *out = fopen(SPTfilename, "w");
        char timefilename[]="time_PMIA_0000.txt";
        char SPT_new_WC[]="GC_PMIA_0000.txt";
        for (int bound=bound1;bound<bound2; bound+=bound){
            printf("%d ",bound);
            double spread, treesize=0;
#ifdef COUNT
            {
                spread=SPT_new::Build(SET_SIZE, bound);
                printf("%lg\n",spread);
                continue;
            }
#endif
            sprintf(timefilename,"time_PMIA_%04d.txt", bound);
            sprintf(SPT_new_WC,"GC_PMIA_%04d.txt",bound);
            {
                start = clock();
                treesize=SPT_new::Build(SET_SIZE, bound);
                ended = clock();
                printf("\n");
                timer = (double)(ended - start);
                timetmpfile = fopen(timefilename, "w");
                fprintf(timetmpfile,"%lg\n", timer);
                fclose(timetmpfile);
                spread=toSimulateOnce(SET_SIZE, SPT_new::GetNode, GeneralCascade::Run);
                toSimulate(SPT_new_WC, SPT_new::GetNode, GeneralCascade::Run);
            }
            fprintf(out,"%lg, %lg, %d, %lg, %d\n", timer/CLOCKS_PER_SEC, spread, bound, treesize, SET_SIZE);
        }
        fclose(out);
    }

    //control bound to test MIA_GC
    s="-m";
    if (!s.compare(argv[1])) {
        Graph::Build2WC();
        GeneralCascade::Build();
        int bound1=10, bound2=2000;
        if (argc>=3) sscanf(argv[2],"%d",&bound1);
        if (argc>=4) sscanf(argv[3],"%d",&bound2);
        char SPTfilename[]="MIA_control.txt";
        FILE *out = fopen(SPTfilename, "w");
        char timefilename[]="time_MIA_0000.txt";
        char SPT_new_WC[]="GC_MIA_0000.txt";
        for (int bound=bound1;bound<bound2; bound+=bound){
            printf("%d ",bound);
            double spread, treesize=0;
#ifdef COUNT
            {
                spread=SPT_new::Build(SET_SIZE, bound);
                printf("%lg\n",spread);
                continue;
            }
#endif
            sprintf(timefilename,"time_MIA_%04d.txt", bound);
            sprintf(SPT_new_WC,"GC_MIA_%04d.txt",bound);
            {
                start = clock();
                treesize=MIA::Build(SET_SIZE, bound);
                ended = clock();
                printf("\n");
                timer = (double)(ended - start);
                timetmpfile = fopen(timefilename, "w");
                fprintf(timetmpfile,"%lg\n", timer);
                fclose(timetmpfile);
                spread=toSimulateOnce(SET_SIZE, MIA::GetNode, GeneralCascade::Run);
                toSimulate(SPT_new_WC, MIA::GetNode, GeneralCascade::Run);
            }
            fprintf(out,"%lg %lg %d %lg\n", timer, spread, bound, treesize);
        }
        fclose(out);
    }

}

