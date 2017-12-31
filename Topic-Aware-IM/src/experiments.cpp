#include "stdio.h"
#include "time.h"
#include <string>
#include <math.h>

#include "tim.h"

void PrintS(S s) {
    printf("S = \n");
    for (unsigned int i = 0; i < s.ids.size(); i++) {
        printf("%d\t%lg\n", s.ids[i], s.minflu[i]);
    }
}

int
main(int argc, char *argv[]) {
    clock_t start, end;
    double timer;

    string cmd;
    cmd = "-p";
    if (!cmd.compare(argv[1])) {
        TIC::Build2TICFromFile(argv[2]);
        // load query items
        char dm[] = "very_very_very_very_very_long_dm_file_name";
        FILE *f = fopen(argv[3], "r");
        int nq = 0, Z = 0, k;
        int q_id;
        sscanf(argv[5], "%d", &q_id);
        vector<double> item;
        fscanf(f, "%d %d", &nq, &Z);
        item.resize(Z);
        for (int i = 0; i < nq; i++) {
            fscanf(f, "%d", &k);
            for (int z = 0; z < Z; z++) {
                fscanf(f, "%lg", &(item[z]));
            }
            if (q_id == i) {
                sprintf(dm, "%s", argv[4]);
//                printf("%s\n", dm);
                start = clock();
                TIC::Convert2PMIA(item, dm);
#ifdef VERBO    SE
                printf("convert time %lg\n", double(clock()-start)/CLOCKS_PER_SEC);
#endif
#ifdef EXP
                printf("%lg\n", double(clock()-start)/CLOCKS_PER_SEC);
#endif
            }
        }
        fclose(f);
    }

    bool bounded = false;
    cmd = "-b";
    string cmd2 = "-f";
    if (!cmd.compare(argv[1]) || !cmd2.compare(argv[1])) {
        if (!cmd.compare(argv[1])) bounded = true;
        else bounded = false;
//        printf("bounded? %d", bounded ? 1 : 0);
        TIM::LoadTIC(argv[2]);
        TIM::InitGreedy(argv[3]);
        int bound1, bound2;
        sscanf(argv[4], "%d", &bound1);
        sscanf(argv[5], "%d", &bound2);
        TIM::theta_ = 1.0/bound1;
        TIM::theta2_ = 1.0/bound2;

        // load query items
        FILE *f = fopen(argv[6], "r");
        int nq = 0, Z = 0, k;
        sscanf(argv[7], "%d", &k);
        vector<double> item;
        fscanf(f, "%d %d", &nq, &Z);
        item.resize(Z);
        for (int i = 0; i < nq; i++) {
            start = clock();

            for (int z = 0; z < Z; z++) {
                fscanf(f, "%lg", &(item[z]));
            }
            Q q = {.item=item, .k=k};
            S s;
            if (bounded) s = TIM::Greedy(q);
            else s = TIM::GreedyNoBound(q);
#ifdef VERBOSE
            PrintS(s);
#endif
            end = clock();
            timer = double(end - start);

            double spread = TIC::Spread(s.ids, q.item);
            
            double rt = 0;
            for	(unsigned int i = 0; i < s.ids.size(); i++) {
                rt += s.minflu[i];
            }

#ifdef VERBOSE
            printf("spread : %lg, time: %lg\n", spread, timer/CLOCKS_PER_SEC);
#endif
#ifdef EXP
            for (int z = 0; z < Z; z++) {
                printf("%lg", item[z]);
            }
            printf("\n");
            printf("%lg, %lg, %lg, %d, %d, %d\n", timer/CLOCKS_PER_SEC, spread, rt, bound1, bound2, k);
#endif
        }
        fclose(f);
    }

    cmd = "-a";
    if (!cmd.compare(argv[1])) {
        TIM::LoadTIC(argv[2]);
        TIM::InitApprox(argv[3], argv[4], argv[5]);
        int bound1, bound2; double epsilon;
        sscanf(argv[6], "%d", &bound1);
        sscanf(argv[7], "%d", &bound2);
        TIM::theta_ = 1.0/bound1;
        TIM::theta2_ = 1.0/bound2;

        // load query items
        FILE *f = fopen(argv[8], "r");
        sscanf(argv[9], "%lg", &epsilon);
        int nq = 0, Z = 0, k;
        sscanf(argv[10], "%d", &k);
        vector<double> item;
        fscanf(f, "%d %d", &nq, &Z);
        item.resize(Z);
        double placeholder;
        for (int i = 0; i < nq; i++) {
            start = clock();
            fscanf(f, "%d", &placeholder);
            for (int z = 0; z < Z; z++) {
                fscanf(f, "%lg", &(item[z]));
            }
            Q q = {.item=item, .k=k};
            S s = TIM::Approx(q, epsilon);

#ifdef VERBOSE
            PrintS(s);
#endif

            end = clock();
            timer = double(end - start);

            double spread = TIC::Spread(s.ids, q.item);

#ifdef VERBOSE
            printf("spread : %lg, time: %lg\n", spread, timer/CLOCKS_PER_SEC);
#endif
#ifdef EXP
            for (int z = 0; z < Z; z++) {
                printf("%lg", item[z]);
            }
            printf("\n");
            printf("%lg, %lg, %d, %d, %d, %lg\n", timer/CLOCKS_PER_SEC, spread, bound1, bound2, k, epsilon);
#endif
        }
        fclose(f);
    }
}
