#include "user_defined.h"
#include <stdio.h>
#define N 2

void sum (Full_particle x[N]) {
    for (int i=0; i<N;++i){
        x[i].epot += (double)i + 1;
    }
}
void more_sum (Full_particle x[N]) {
    for (int i=0; i<N;++i){
        x[i].epot += (double)i + 1;
    }
}
int main () {
    Full_particle ptlc[N] = {10};
    // sum(ptlc);
    // more_sum(ptlc);
    for (int i=0; i<N;++i){
        printf("%f\n", ptlc[i].epot);
    }
}