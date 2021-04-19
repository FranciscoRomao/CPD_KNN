#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(){
    omp_set_nested(2);
    #pragma omp parallel
    {
        #pragma omp for
        for(int i=0;i<10;i++){
                printf("omp_get_thread_num(): %d\n",omp_get_thread_num());
        }
    }

}