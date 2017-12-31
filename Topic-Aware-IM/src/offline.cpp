// offline methods

// #define VERBOSE

#include "stdio.h"
#include <string>
#include <math.h>
#include <algorithm>    // std::sort

#include "tim.h"

typedef struct
{
    int id;
    double infl;
    double infl2;  // theta / theta2_
} LNode;

// bool max_lnode_cmp(const di&a, const di&b) {
// //    return a.infl > b.infl;
//     return a.first > b.first;
// }

bool max_lnode_cmp(const LNode&a, const LNode&b) {
    return a.infl > b.infl;
}

void List(double theta, double theta2, char *filename) {
    double delta = theta / theta2;  // 0.01 / 0.1 = 0.1
    vector<LNode> l;
    for (int i = 0; i < TIC::GetN(); i++) {
        if (i % 1000 == 0) printf("processing vertex %d...\n", i);
        vector<di> rst = TIC::Dijkstra(i, -log(theta));
        double rt = 0, rt2 = 0;
        for (unsigned int i = 0; i < rst.size(); i++) {
            if (rst[i].first > delta) rt2 += rst[i].first;
            rt += rst[i].first;
        }
        LNode ln = {.id = i, .infl=rt, .infl2 = rt2};
        l.push_back(ln);
    }
    sort(l.begin(), l.end(), max_lnode_cmp);
    FILE *f = fopen(filename, "w");
    fprintf(f, "%d\n", TIC::GetN());
    for (int i = 0; i < TIC::GetN(); i++) {
        fprintf(f, "%d\t%lg\t%lg\n", l[i].id, l[i].infl, l[i].infl2);
    }
    fclose(f);
}

void UpSample(int k, char *filename, char *outname) {
    printf("up sampling...\n");
    FILE* f = fopen(filename, "r");
    FILE *fo = fopen(outname, "w");
    int nq = 0, Z = 0;
    vector<double> item;
    fscanf(f, "%d %d", &nq, &Z);
    fprintf(fo, "%d\t%d\t%d\n", nq, Z, k);
    item.resize(Z);
    double tmp;
    for (int i = 0; i < nq; i++) {
#ifdef VERBOSE
        printf("processing sample %d...\n", i);
#endif
        fscanf(f, "%lg", &tmp);
        for (int z = 0; z < Z; z++) {
            fscanf(f, "%lg", &(item[z]));
        }
        Q q = {.item=item, .k=k};
        S s = TIM::Greedy(q);
        double spread = 0;
        for (unsigned int i = 0; i < s.minflu.size(); i++) {
            spread += s.minflu[i];
        }
        fprintf(fo, "%lg\n", spread);
        for (int z = 0; z < Z; z++) {
            fprintf(fo, "%lg ", item[z]);
        }
        fprintf(fo, "\n");
    }
    fclose(f);
    fclose(fo);
}

void LowSample(int k, char *filename, char *outname) {
    printf("low sampling...\n");
    FILE* f = fopen(filename, "r");
    FILE *fo = fopen(outname, "w");
    int nq = 0, Z = 0;
    vector<double> item;
    fscanf(f, "%d %d", &nq, &Z);
    fprintf(fo, "%d\t%d\t%d\n", nq, Z, k);
    item.resize(Z);
    for (int i = 0; i < nq; i++) {
        if (i % 100 == 0) printf("processing sample %d...\n", i);
        for (int z = 0; z < Z; z++) {
            fscanf(f, "%lg", &(item[z]));
        }
        Q q = {.item=item, .k=k};
        S s = TIM::Greedy(q);
//        vector<int> seeds = TIM::LBGenerate(s.ids);

        for (int z = 0; z < Z; z++) {
            fprintf(fo, "%lg ", item[z]);
        }
        fprintf(fo, "\n");
        for (int j = 0; j < k; j++) {
            fprintf(fo, "%d ", s.ids[j]);
        }
        fprintf(fo, "\n");
//        for (int i = 0; i < k; i++) {
//            fprintf(fo, "%d ", lbouts[i].size());
//            for (unsigned int j = 0; j < lbouts[i].size(); j++) {
//                fprintf(fo, "%d ", lbouts[i][j]);
//            }
//            fprintf(fo, "\n");
//        }
    }
    fclose(f);
    fclose(fo);
}

void UpHeapSample() {
}

int
main(int argc, char *argv[]) {
    string cmd;
    cmd = "-l";
    if (!cmd.compare(argv[1])) {
        TIC::Build2TICFromFile(argv[2]);
        int bound1 = 100, bound2 = 10;
        sscanf(argv[3], "%d", &bound1);
        sscanf(argv[4], "%d", &bound2);
        List(1.0/bound1, 1.0/bound2, argv[5]);
    }

    cmd = "-ls";
    if (!cmd.compare(argv[1])) {
        TIM::LoadTIC(argv[2]);
        TIM::InitGreedy(argv[3]);
        int bound1 = 1000, bound2=10, k = 1000;
        sscanf(argv[4], "%d", &bound1);
        sscanf(argv[5], "%d", &bound2);
        sscanf(argv[6], "%d", &k);
        TIM::theta_ = 1.0/bound1;
        TIM::theta2_ = 1.0/bound2;
        LowSample(k, argv[7], argv[8]);
    }

    cmd = "-us";
    if (!cmd.compare(argv[1])) {
        TIM::LoadTIC(argv[2]);
        TIM::InitGreedy(argv[3]);
        int bound1 = 1000, bound2=10, k = 1000;
        sscanf(argv[4], "%d", &bound1);
        sscanf(argv[5], "%d", &bound2);
        sscanf(argv[6], "%d", &k);
        TIM::theta_ = 1.0/bound1;
        TIM::theta2_ = 1.0/bound2;
        UpSample(k, argv[7], argv[8]);
    }

    cmd = "-uhs";
    if (!cmd.compare(argv[1])) {
        UpHeapSample();
    }
}
