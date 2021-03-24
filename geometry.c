#include <math.h>
#include <stdlib.h>
#include "geometry.h"

// A utility function to swap two elements
void swap(double **pts, double **projections, double *distances2a, long a, long b)
{
    double *temp1 = pts[a];
    double *temp2 = projections[a];
    double temp3 = distances2a[a];

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
long partition(double **pts, double **projections, double *distances2a, long low, long high)
{
    double pivot = distances2a[high]; // pivot
    long i = (low - 1);               // Index of smaller element and indicates the right position of pivot found so far

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
void quick_sort(double **pts, double **projections, double *distances2a, long low, long high)
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

long getMedian(double **projections, long n_points, int n_dims, double *center)
{
    long idx;

    if (n_points % 2 != 0) //odd number of points
    {
        idx = (n_points - 1) / 2;
        center = projections[idx];
    }
    else //even number of points
    {
        idx = (n_points - 2) / 2; //point to the left of median coordinates

        for (long i = 0; i < n_dims; i++)
        {
            center[i] = (projections[idx][i] + projections[idx + 1][i]) / 2;
        }
        idx = idx + 1; //point to the right of median coordinates
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
    for (int i = 0; i < n_dims; i++)
    {
        distance += pow(a[i] - b[i], 2);
    }
    distance = sqrt(distance);
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

double *multiply(int n_dims, double *a, double constant)
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
    double numerator = 0;
    double denominator = 0;
    double result_div = 0;
    double *result_mult;
    double *result_sum;
    double *b_minus_a;
    double *p_minus_a;

    b_minus_a = subtraction(n_dims, a, b);
    p_minus_a = subtraction(n_dims, a, p);
    numerator = inner_product(n_dims, p_minus_a, b_minus_a);
    denominator = inner_product(n_dims, b_minus_a, b_minus_a);
    result_div = numerator / denominator;
    result_mult = multiply(n_dims, b_minus_a, result_div);
    result_sum = sum(n_dims, result_mult, a);
    free(b_minus_a);
    free(p_minus_a);
    free(result_mult);
    return result_sum;
}

double **project_pts2line(int n_dims, double *a, double *b, double **pts, long n_points)
{
    double **projected_points = (double **)malloc(n_dims * n_points * sizeof(double));
    for (int i = 0; i < n_points; i++)
    {
        projected_points[i] = orthogonal_projection(n_dims, pts[i], a, b);
    }
    return projected_points;
}
