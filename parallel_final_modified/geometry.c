#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "geometry.h"

typedef struct max_n_idx{
    double maximum;
    long int index; 
}max_n_idx;

max_n_idx max_n_idx_max(max_n_idx a, max_n_idx b) {
    return a.maximum > b.maximum ? a : b; 
}

#pragma omp declare reduction(max_n_idx_max_: struct max_n_idx: omp_out=max_n_idx_max(omp_out, omp_in))


/**
* Computes euclidean distance
 * @param *a point a
 * @param *b point b
 * @return euclidean distance between a and b
 */ 
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

/**
 * Computes euclidean distance without sqrt in the end for performancwe purposes
 * @param *a point a
 * @param *b point b
 * @return squared euclidean distance between a and b
 */
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

/**
 * Computes the furthest point in the set from a defined point
 * @param n_dims number of dimensions
 * @param n_point number of points in the set
 * @param **pts set of the points
 * @param coordinates of the point to search from
 */
long furthest_point_from_coords(int n_dims, long n_points, double **pts, double *base_coords, int threads_available)
{
    double curr_dist = 0;
    max_n_idx max_point={-1.0,-1};

    if (threads_available > 1)
    {
        //#pragma omp taskloop reduction(max_n_idx_max_:max_point)
        #pragma omp parallel for reduction(max_n_idx_max_:max_point) 
        for (long i = 0; i < n_points; i++)
        {   
            if ((curr_dist = squared_distance(n_dims, base_coords, pts[i])) > max_point.maximum)
            {
                max_point.maximum = curr_dist;
                max_point.index = i;
            }
        }
    }
    else 
    {
        for (long i = 0; i < n_points; i++)
        {   
            if ((curr_dist = squared_distance(n_dims, base_coords, pts[i])) > max_point.maximum)
            {
                max_point.maximum = curr_dist;
                max_point.index = i;
            }
        }
    }
    return max_point.index;
}


/**
 * Computes the two points furthest apart in the set 
 * (might get local solution but it is not relevant for the problem)
 * @param n_dims   # of dimensions
 * @param n_points # of points 
 * @param **pts array with the points
 * @param *idx_fp 
 */
void recursive_furthest_apart(int n_dims, long n_points, double **pts, long *idx_fp, int threads_available)
{
    long idx_new_fp = 0;

    idx_new_fp = furthest_point_from_coords(n_dims, n_points, pts, pts[idx_fp[0]], threads_available);
    if (idx_new_fp != idx_fp[1])
    {
        idx_fp[1] = idx_fp[0];
        idx_fp[0] = idx_new_fp;
        recursive_furthest_apart(n_dims, n_points, pts, idx_fp, threads_available);
    }
    return;
}

/**Subtracts every dimensions of two points
 * n_dims: # of dimensions
 * a : first point
 * b : second point
 * result : array with the proper dimensions
 * [return] : result
 */
double *subtraction(int n_dims, double *a, double *b, double* result)
{
    for (int i = 0; i < n_dims; i++)
    {
        result[i] = b[i] - a[i];
    }
    return result;
}

/**Sums every dimensions of two points
 * n_dims: # of dimensions
 * a : first point
 * b : second point
 * result : array with the proper dimensions
 * [return] : result
 */
double *sum(int n_dims, double *a, double *b, double *result)
{
    for (int i = 0; i < n_dims; i++)
    {
        result[i] = a[i] + b[i];
    }
    return result;
}

/**Computes the inner product of two points
 * n_dims: # of dimensions
 * a : first point
 * b : second point
 * [return] : inner product
 */
double inner_product(int n_dims, double *a, double *b)
{
    double result = 0;
    for (int i = 0; i < n_dims; i++)
    {
        result += a[i] * b[i];
    }
    return result;
}


/**
 * Multiplies every dimensions of one points with a constant
 * @param   n_dims: # of dimensions
 * @param   a : first point
 * @param   constant : constant
 * @param   result : array with the proper dimensions
 * @param   [return] : result
 */
double *multiply(int n_dims, double *a, double constant, double* result)
{
    for (int i = 0; i < n_dims; i++)
    {
        result[i] = constant * a[i];
    }
    return result;
}

/**
 * Computes the reduced orthogonal projection(used to order the points without computing the full projections)
 * @param   n_dims  # of dimensions
 * @param   p   point to project
 * @param   first   points of the line
 * @param   p_minus_a   point p subtracted by point a
 * @param   b_minus_a   point b (end of line) subtracted by point a
 * @return  result of the projection
 */
double orthogonal_projection_reduced(int n_dims, double *p, double *a, double* p_minus_a, double* b_minus_a)
{
    
    //printf("%lf %x %lf %x %x %d\n", *p, p, a, *a, p_minus_a, omp_get_thread_num());
    
    subtraction(n_dims, a, p, p_minus_a);
    
    double result=inner_product(n_dims, p_minus_a, b_minus_a);
    return result;
}


/**
 * Computes the full orthogonal projection
 * @param   n_dims  # of dimensions
 * @param   p   point to project
 * @param   a   first points of the line
 * @param   b   point of the end of the line
 * @param   result_sum  its used to return by reference the reduced orthogonal projection
*/
void orthogonal_projection(int n_dims, double *p, double *a, double *b, double* result_sum)
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

/**Given a list of points with multiple dimensions this function computes the reduced orthogonal projections
 * for all the points 
 * n_dims: # of dimensions
 * p : point to project
 * a : first points of the line
 * p_minus_a : point p subtracted by point a
 * b_minus_a : point b (end of line) subtracted by point a
 * [return] : result
 */
void project_pts2line(int n_dims, double* projections, double *a, double *b, double **pts, long n_points, int threads_available)
{   
    //printf("Projection\n");
    double *b_minus_a = (double *)malloc(n_dims * sizeof(double));
    double *p_minus_a = NULL;
    if(a[0]>b[0])
        subtraction(n_dims, b, a, b_minus_a);
    else
        subtraction(n_dims, a, b, b_minus_a);
    
    if (threads_available > 1)
    { 
        #pragma omp parallel firstprivate(p_minus_a)
        {
            p_minus_a = (double *)malloc(n_dims * sizeof(double));
            #pragma omp for
            for (int i = 0; i < n_points; i++)
            {   
                //printf("i:%d  id: %d\n\n",i,omp_get_thread_num());

                // if (p_minus_a == NULL)
                //     exit(-1);
                projections[i] = orthogonal_projection_reduced(n_dims, pts[i],a,p_minus_a,b_minus_a);
            }
            free(p_minus_a); 
        }
    }
    else
    {
        p_minus_a = (double *)malloc(n_dims * sizeof(double));
        for (int i = 0; i < n_points; i++)
        {   
            projections[i] = orthogonal_projection_reduced(n_dims, pts[i],a,p_minus_a,b_minus_a);
        }
        free(p_minus_a);
    }
   
    free(b_minus_a);
    
    return;
}