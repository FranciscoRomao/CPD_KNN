#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "geometry.h"

void flag(int n)
{
    printf("Flag %d\n", n);
    fflush(stdout);
}
/*
void swap(double **pts, double *projections, long a, long b)
{   
    int n_dims=2;
    double* copy=(double*)malloc(n_dims*sizeof(double));
    for(int i=0;i<n_dims;i++){
        copy[i]=pts[a][i];
    }
    for(int i=0;i<n_dims;i++){
        pts[a][i]=pts[b][i];
    }
    for(int i=0;i<n_dims;i++){
        pts[b][i]=copy[i];
    }
    double temp2 = projections[a];
    projections[a] = projections[b];
    projections[b] = temp2;

}
*/

// A utility function to swap two elements
void swap(double **pts, double *projections, long a, long b)
{   
    
    double *temp1 = pts[a];
    double temp2 = projections[a];

    pts[a] = pts[b];
    pts[b] = temp1;

    projections[a] = projections[b];
    projections[b] = temp2;

}

/* This function takes last element as pivot, places 
the pivot element at its correct position in sorted 
array, and places all smaller (smaller than pivot) 
to left of pivot and all greater elements to right 
of pivot */
long partition(double **pts, double *projections, long low, long high)
{
    double pivot = projections[high]; // pivot
    long i = (low - 1);               // Index of smaller element and indicates the right position of pivot found so far

    for (long j = low; j <= high - 1; j++)
    {
        // If current element is smaller than the pivot
        if (projections[j] < pivot)
        {
            i++; // increment index of smaller element
            swap(pts, projections, i, j);
        }
    }
    swap(pts, projections, i + 1, high);
    return (i + 1);
}

/* The main function that implements QuickSort 
arr[] --> Array to be sorted, 
low --> Starting index, 
high --> Ending index */
void quick_sort(double **pts, double *projections, long low, long high)
{
    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now 
        at right place */
        long pi = partition(pts, projections, low, high);

        // Separately sort elements before
        // partition and after partition
        quick_sort(pts, projections, low, pi - 1);
        quick_sort(pts, projections, pi + 1, high);
    }
}

long getMedian(double **pts, long n_points, int n_dims, double *center)
{
    long idx;
    if (n_points % 2 != 0) //odd number of points
    {
        idx = (n_points - 1) / 2;
        orthogonal_projection(n_dims, pts[idx], pts[0], pts[n_points-1], center);
    }
    else //even number of points
    {
        double* center1=(double*)malloc(n_dims*sizeof(double));
        double* center2=(double*)malloc(n_dims*sizeof(double));
        
        idx = (n_points - 2) / 2; //point to the left of median coordinates
        orthogonal_projection(n_dims, pts[idx], pts[0], pts[n_points-1], center1);
        orthogonal_projection(n_dims, pts[idx+1], pts[0], pts[n_points-1], center2);
        
        for(int i=0; i<n_dims; i++)
        {
            center[i] = (center1[i] + center2[i]) / 2;
        }

        idx = idx + 1; //point to the right of median coordinates
        
        free(center1);
        free(center2);    
    }

    
    return idx;
}

double *calc_distances_to_left_limit(double *left_limmit, double **projections, long n_points, int n_dims)
{
    double *dists = malloc(sizeof(double) * n_points);

    for (long i = 0; i < n_points; i++)
    {
        dists[i] = distance(n_dims, left_limmit, projections[i]);
    }

    return dists;
}

double distance(int n_dims, double *a, double *b)
{
    double distance = 0;
    double aux = 0;
    for (int i = 0; i < n_dims; i++)
    {   
        aux = a[i] - b[i];
        distance += aux * aux;
    }
    distance = sqrt(distance);
    return distance;
}

double squared_distance(int n_dims, double *a, double *b)
{
    double distance = 0;
    double aux = 0;
    for (int i = 0; i < n_dims; i++)
    {   
        aux = a[i] - b[i];
        distance += aux * aux;
    }
    return distance;
}

long furthest_point_from_coords(int n_dims, long n_points, double **pts, double *base_coords)
{
    double max_dist = -1, curr_dist = 0;
    long idx_newpt = 0;

    for (long i = 0; i < n_points; i++)
    {
        if ((curr_dist = distance(n_dims, base_coords, pts[i])) > max_dist)
        {
            max_dist = curr_dist;
            idx_newpt = i;
        }
    }
    return idx_newpt;
}

void recursive_furthest_apart(int n_dims, long n_points, double **pts, long *idx_fp)
{
    long idx_new_fp = 0;

    idx_new_fp = furthest_point_from_coords(n_dims, n_points, pts, pts[idx_fp[0]]);
    if (idx_new_fp != idx_fp[1])
    {
        idx_fp[1] = idx_fp[0];
        idx_fp[0] = idx_new_fp;
        recursive_furthest_apart(n_dims, n_points, pts, idx_fp);
    }
    return;
}


double *subtraction(int n_dims, double *a, double *b, double* result)
{
    for (int i = 0; i < n_dims; i++)
    {
        result[i] = b[i] - a[i];
    }
    return result;
}

double *sum(int n_dims, double *a, double *b, double * result)
{
    for (int i = 0; i < n_dims; i++)
    {
        result[i] = a[i] + b[i];
    }
    return result;
}

double inner_product(int n_dims, double *a, double *b)
{
    double result = 0;
    for (int i = 0; i < n_dims; i++)
    {
        result += a[i] * b[i];
    }
    return result;
}

double *multiply(int n_dims, double *a, double constant, double* result)
{
    for (int i = 0; i < n_dims; i++)
    {
        result[i] = constant * a[i];
    }
    return result;
}

double orthogonal_projection_reduced(int n_dims, double *p, double *a, double* p_minus_a, double* b_minus_a)
{
    subtraction(n_dims, a, p, p_minus_a);
    double result=inner_product(n_dims, p_minus_a, b_minus_a);
    return result;
}

void orthogonal_projection(int n_dims, double *p, double *a, double *b,double* result_sum)
{
    double numerator = 0;
    double denominator = 0;
    double result_div = 0;
    double *result_mult = (double *)malloc(n_dims * sizeof(double));
    double *b_minus_a = (double *)malloc(n_dims * sizeof(double));
    double *p_minus_a = (double *)malloc(n_dims * sizeof(double));

    
    subtraction(n_dims, a, b, b_minus_a);
    subtraction(n_dims, a, p, p_minus_a);
    numerator = inner_product(n_dims, p_minus_a, b_minus_a);
    denominator = inner_product(n_dims, b_minus_a, b_minus_a);
    result_div = numerator / denominator;
    multiply(n_dims, b_minus_a, result_div, result_mult);
    sum(n_dims, result_mult, a, result_sum);
    free(b_minus_a);
    free(p_minus_a);
    free(result_mult);
    return;
}

void project_pts2line(int n_dims, double* projections, double *a, double *b, double **pts, long n_points)
{
    double *p_minus_a = (double *)malloc(n_dims * sizeof(double));
    double *b_minus_a = (double *)malloc(n_dims * sizeof(double));
    
    if(a[0]>b[0])
        subtraction(n_dims, b, a, b_minus_a);
    else
        subtraction(n_dims, a, b, b_minus_a);

    for (int i = 0; i < n_points; i++)
    {   
        projections[i] = orthogonal_projection_reduced(n_dims, pts[i],a,p_minus_a,b_minus_a);
    }

    free(p_minus_a);
    free(b_minus_a);
    
    return;
}