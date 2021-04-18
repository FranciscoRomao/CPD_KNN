#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <time.h>
#include <omp.h>

void printArray(double* arr,int size){
    printf("\n");
    for(int i=0;i<size;i++){
        printf(" %f", arr[i]);
    }
    printf("\n");
}

void setUp(double* a, int size){
    int i;

    srand(time(NULL));
    for (i = 0; i<size; ++i) {
        a[i] = rand();
    }
    return;
}

int check_min(double* a, int size){
    double min=DBL_MAX;
    for(int i=0;i<size-1;i++){
        if(a[i]<min)
            min=a[i];
    }
    printf("\n\nMINIMUM IS: %f",min);
}

int main(){
    double arr[1000];
    setUp(arr,1000);
    //printArray(arr,1000);

    double min = DBL_MAX;

    #pragma omp parallel for reduction(min:min)
    for(int i=0; i<1000; i++)
    {
        if(arr[i] < min)
            min = arr[i];
    }
    check_min(arr,1000);
    printf("\n\nPARALLEL MINIMUM IS %f",min);
    return min;
}