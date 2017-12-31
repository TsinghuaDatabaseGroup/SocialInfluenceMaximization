#include "tim.h"

int main() {
    TIC::Build2TICFromFile("./lastfm.tedge");
    TIC::Stats();
    TIC::Print();
    TIC::SaveConverter("./lastfm.cvt");

    int vertex = 23;
    printf("%d's out neighbor\n", vertex);
    for (int i = 0; i < TIC::GetOutNeighbor(vertex); i++) {
        Edge e = TIC::GetOutEdge(vertex, i);
        printf("%d\t%d", e.u, e.v);
        for (int z = 0; z < TIC::GetZ(); z++) {
            printf("\t%lg", e.w[z]);
        }
        printf("\n");
    }
    printf("%d's in neighbor\n", vertex);
    for (int i = 0; i < TIC::GetInNeighbor(vertex); i++) {
        Edge e = TIC::GetInEdge(vertex, i);
        printf("%d\t%d", e.u, e.v);
        for (int z = 0; z < TIC::GetZ(); z++) {
            printf("\t%lg", e.w[z]);
        }
        printf("\n");
    }
}
