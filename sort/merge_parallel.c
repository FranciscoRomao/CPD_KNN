#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>

// Merges two subarrays of arr[].
// First subarray is arr[l..m]
// Second subarray is arr[m+1..r]
void merge(double arr[], int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;
 
    /* create temp arrays */
    double* L =(double*)malloc(n1*sizeof(double));
    double* R =(double*)malloc(n2*sizeof(double));
 
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];
 
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }
 
    /* Copy the remaining elements of L[], if there
    are any */
    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }
 
    /* Copy the remaining elements of R[], if there
    are any */
    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
    free(L);
    free(R);
}
 
/* l is for left index and r is right index of the
sub-array of arr to be sorted */
void mergeSort(double arr[], int l, int r)
{
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        int m = l + (r - l) / 2;
        #pragma omp parallel
        {
            #pragma omp sections
            {       
                
                // Sort first and second halves
                #pragma omp section
                mergeSort(arr, l, m);
                #pragma omp section
                mergeSort(arr, m + 1, r);
            }
            #pragma omp single nowait
            merge(arr, l, m, r); 
        }        
    }
}

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
        a[i] = rand() % size;
    }
    return;
}

int check_sorted(double* a, int size){
    for(int i=0;i<size-1;i++){
        if(a[i+1]<a[i]){
            printf("\n\n\na[i+1]:%f a[i]:%f\n",a[i+1],a[i]);
            return 0;
        }
    }
    return 1;
}

/* Driver code */
int main()
{      
    double exec_time;
    int arr_size=1000000;
    double* arr=(double*)malloc(arr_size*sizeof(double));
    if(arr==NULL){
        printf("Error allocating memory");
        return -1;
    }
    setUp(arr,arr_size);
    //int arr_size = sizeof(arr) / sizeof(arr[0]);
 
    printf("Given array is \n");
    printArray(arr, arr_size);
    
    printf("\n\n\n\n");
    exec_time = -omp_get_wtime();
    mergeSort(arr, 0, arr_size - 1);
    exec_time += omp_get_wtime();
    
    printf("\nSorted array is \n");
    printArray(arr, arr_size);

    if(check_sorted(arr,arr_size))
        printf("Array is sorted!\n\n");
    printf("%.1lf\n", exec_time); 
    free(arr);
    return 0;
}