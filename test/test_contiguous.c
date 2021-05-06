#include <stdio.h>
#include <stdlib.h>

int main(){

    int n_dims=3;
    int n_points=4;
    double ** pts;
    pts=(double**)malloc(n_points*sizeof(double**));
    for(int i=0;i<n_points;i++){
        pts[i]=(double*)malloc(n_dims*sizeof(double));
        printf("%x\n",pts[i]);
    }
}