#include <math.h>
#include <stdlib.h>
#include "geometry.h"

// A utility function to swap two elements
void swap(double **pts, double **projections, double **distances2a, long a, long b)
{
    double *temp1 = pts[a];
    double *temp2 = projections[a];
    double *temp3 = distances2a[a];

    pts[a] = pts[b];
    pts[b] = temp1;

    projections[a] = projections[b];
    projections[b] = temp2;

    distances2a[a] = distances2a[b];
    distances2a[b] = temp3;
}

/* This function takes last element as pivot, places 
the pivot element at its correct position in sorted 
array, and places all smaller (smaller than pivot) 
to left of pivot and all greater elements to right 
of pivot */
long partition(double **pts, double **projections, double **distances2a, long low, long high)
{
    long pivot = distances2a[high]; // pivot
    long i = (low - 1);             // Index of smaller element and indicates the right position of pivot found so far

    for (long j = low; j <= high - 1; j++)
    {
        // If current element is smaller than the pivot
        if (distances2a[j] < pivot)
        {
            i++; // increment index of smaller element
            swap(pts, projections, distances2a, i, j);
        }
    }
    swap(pts, projections, distances2a, i + 1, high);
    return (i + 1);
}

/* The main function that implements QuickSort 
arr[] --> Array to be sorted, 
low --> Starting index, 
high --> Ending index */
void quick_sort(double **pts, double **projections, double **distances2a, long low, long high)
{
    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now 
        at right place */
        long pi = partition(pts, projections, distances2a, low, high);

        // Separately sort elements before
        // partition and after partition
        quick_sort(pts, projections, distances2a, low, pi - 1);
        quick_sort(pts, projections, distances2a, pi + 1, high);
    }
}

long getMedian(double **projections, long npts, int ndims, double *center)
{
    long idx;

    if (npts % 2 != 0) //odd number of points
    {
        idx = (npts - 1) / 2;
        center = projections[idx];
    }
    else //even number of points
    {
        idx = (npts - 2) / 2; //point to the left of median coordinates

        for (long i = 0; i < ndims; i++)
        {
            center[i] = (projections[idx][i] + projections[idx + 1][i]) / 2;
        }
        idx = idx + 1; //point to the right of median coordinates
    }
    return idx;
}

double *calc_distances_to_left_limit(double *left_limmit, double **projections, long npts, int ndims)
{
    double *dists = malloc(sizeof(double) * npts);

    for (long i = 0; i < npts; i++)
    {
        dists[i] = distance(ndims, left_limmit, projections[i]);
    }

    return dists;
}

double distance(int n_dims, double *a, double *b)
{
    double distance = 0;
    for (int i = 0; i < n_dims; i++)
    {
        distance += pow(a[i] - b[i], 2);
    }
    distance = sqrt(distance);
    return distance;
}

long furthest_point(int n_dims, long np, double **pts, long base_idx)
{
    double max_dist = -1, curr_dist = 0;
    long idx_newpt = 0;

    for (long i = 0; i < np; i++)
    {
        if ((curr_dist = distance_(n_dims, pts[base_idx], pts[i])) > max_dist)
        {
            max_dist = curr_dist;
            idx_newpt = i;
        }
    }
    return idx_newpt;
}

long *recursive_furthest_apart(int n_dims, long np, double **pts, long *idx_fp)
{
    long idx_new_fp = 0;
    idx_new_fp = furthest_point(n_dims, np, **pts, idx_fp[0]);
    if (idx_new_fp != idx_fp[1])
    {
        idx_fp[1] = idx_fp[0];
        idx_fp[0] = idx_new_fp;
        idx_fp = recursive_furthest_apart(n_dims, long np, double **pts, idx_fp)
    }
    return idx_fp;
}

double *subtraction(int n_dims, double *a, double *b)
{
    double *result = (double *)malloc(n_dims * sizeof(double));
    for (int i = 0; i < n_dims; i++)
    {
        result[i] = b[i] - a[i];
    }
    return result;
}

double *sum(int n_dims, double *a, double *b)
{
    double *result = (double *)malloc(n_dims * sizeof(double));
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

double *multiply(int n_dims, double *a, int constant)
{
    double *result = (double *)malloc(n_dims * sizeof(double));
    for (int i = 0; i < n_dims; i++)
    {
        result[i] = constant * a[i];
    }
    return result;
}

double *orthogonal_projection(int n_dims, double *p, double *a, double *b)
{
    double numerator = 0, denominator = 0, result_div = 0;
    double *result = (double *)malloc(n_dims * sizeof(double));
    double *b_minus_a = (double *)malloc(n_dims * sizeof(double));

    b_minus_a = subtraction(n_dims, a, b);
    numerator = inner_product(n_dims, subtraction(n_dims, a, p), b_minus_a);
    denominator = inner_product(n_dims, b_minus_a, b_minus_a);
    result_div = numerator / denominator;
    result = multiply(n_dims, b_minus_a, result);
    result = sum(n_dims, result, a);
    return result;
}

double **project_pts2line(int n_dims, /*[numero de pontos]*/
                          double *a,
                          double *b,
                          double **pts,
                          long *idx2project,
                          long np2project)
{
    double **projected_points = (double *)malloc(n_dims * np2project * sizeof(double));
    for (int i = 0; i < np2project; i++)
    {
        if (idx2project[0] == -1)
            projected_points[i] = orthogonal_projection(n_dims, pts[i], a, b);
        else
            projected_points[i] = orthogonal_projection(n_dims, pts[idx2project[i]], a, b);
    }
    return projected_points;
}



/*
long *recursive_furthest_apart(double **pts, int n_dims, long np, long base_index,
                               double max_distance, long *furthest_pts)
{
    if (base_index == np - 2)
    { //only one point remaining
        return furthest_pts;
    }
    long furthest_pt = furthest_point(n_dims, np, pts, base_index);
    long new_distance = 0;
    if ((new_distance = distance(n_dims, pts[base_index], pts[furthest_pt])) > max_distance)
    {
        max_distance = new_distance;
        furthest_pts[0] = base_index;
        furthest_pts[1] = furthest_pt;
    }
    return recursive_furthest_apart(pts, n_dims, np, base_index + 1, max_distance, furthest_pts);
}

long *furthest_apart(int n_dims, long np, double **pts)
{
    long *furthest_pts = (long *)malloc(2 * sizeof(long));
    double max_distance = -1;
    double curr_distance = 0;
    for (long i = 0; i < np; i++)
    {
        for (long j = i + 1; j < np; j++)
        {
            if ((curr_distance = distance(n_dims, pts[i], pts[j])) > max_distance)
            {
                furthest_pts[0] = i;
                furthest_pts[1] = j;
            }
        }
    }
    return furthest_pts;
}
*/