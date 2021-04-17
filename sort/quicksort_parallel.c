#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>

// A utility function to swap two elements
void swap(double *pts, long a, long b)
{   
    
    double temp1 = pts[a];

    pts[a] = pts[b];
    pts[b] = temp1;

}

/* This function takes last element as pivot, places 
the pivot element at its correct position in sorted 
array, and places all smaller (smaller than pivot) 
to left of pivot and all greater elements to right 
of pivot */
long partition(double *pts, long low, long high)
{
    double pivot = pts[high]; // pivot
    long i = (low - 1);               // Index of smaller element and indicates the right position of pivot found so far

    for (long j = low; j <= high - 1; j++)
    {
        // If current element is smaller than the pivot
        if (pts[j] < pivot)
        {
            i++; // increment index of smaller element
            swap(pts, i, j);
        }
    }
    swap(pts, i + 1, high);
    return (i + 1);
}

/* The main function that implements QuickSort 
arr[] --> Array to be sorted, 
low --> Starting index, 
high --> Ending index */
void quick_sort(double *pts, long low, long high)
{   
    int cutoff=10000;
    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now 
        at right place */
        long pi = partition(pts, low, high);

        // Separately sort elements before
        // partition and after partition
        if(high-low<cutoff){
            quick_sort(pts, low, pi - 1);
            quick_sort(pts, pi + 1, high);
        }
        else
        {   
            #pragma omp task
            {quick_sort(pts, low, pi - 1);}
            #pragma omp task
            {quick_sort(pts,pi + 1, high);}
        }
    }
}

void quick_sort_parallel(double *pts, int n_points){
    #pragma omp parallel
    {
        #pragma omp single nowait
        {
            quick_sort(pts, 0, n_points);
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
    quick_sort_parallel(arr, arr_size - 1);
    exec_time += omp_get_wtime();
    printf("\nSorted array is \n");
    printArray(arr, arr_size);

    if(check_sorted(arr,arr_size))
        printf("Array is sorted!");
    free(arr);
    printf("%.1lf\n", exec_time); 
    return 0;
}