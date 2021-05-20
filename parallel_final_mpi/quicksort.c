#include <stdlib.h>
#include <omp.h>
#include "quicksort.h"

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

void quick_sort_serial(double *pts, long low, long high)
{   
    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now 
        at right place */
        long pi = partition(pts, low, high);

        // Separately sort elements before
        // partition and after partition
        quick_sort_serial(pts, low, pi - 1);
        quick_sort_serial(pts,pi + 1, high);
        
    }
}

/* The main function that implements QuickSort 
arr[] --> Array to be sorted, 
low --> Starting index, 
high --> Ending index */
void quick_sort(double *pts, long low, long high, int cutoff)
{   
    cutoff=10000;
    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now 
        at right place */
        long pi = partition(pts, low, high);

        // Separately sort elements before
        // partition and after partition
        if(high-low<cutoff){
            quick_sort(pts, low, pi - 1,cutoff);
            quick_sort(pts, pi + 1, high,cutoff);
        }
        else
        {   
            #pragma omp task
            {quick_sort(pts, low, pi - 1,cutoff);}
            #pragma omp task
            {quick_sort(pts,pi + 1, high,cutoff);}
        }
    }
}

void quick_sort_parallel(double *pts, int n_points, int cutoff){
    #pragma omp parallel
    {
        #pragma omp single nowait
        {
            quick_sort(pts, 0, n_points,cutoff);
        }
    }
}

