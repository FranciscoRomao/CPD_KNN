#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "geometry.h"
#include <mpi.h>
#include "gen_points.h"
#include "unsorted_median.h"

typedef struct max_n_idx
{
    double maximum;
    long int index; 
}max_n_idx;

max_n_idx max_n_idx_max(max_n_idx a, max_n_idx b)
{
    return a.maximum > b.maximum ? a : b; 
}


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
 * @param   n_dims              number of dimensions
 * @param   n_point             number of points in the set
 * @param   pts                 **pts set of the points
 * @param   base_coords         of the point to search from
 */
long furthest_point_from_coords(int n_dims, long start, long end, double **pts, double *base_coords)
{
    //printf("Coords: %lf %lf\n", base_coords[0], base_coords[1]);
    //fflush(stdout);

    double curr_dist = 0;
    max_n_idx max_point={-1.0,-1};
    
    for (long i = start; i < end; i++)
    {   
        if ((curr_dist = squared_distance(n_dims, base_coords, pts[i])) > max_point.maximum)
        {
            max_point.maximum = curr_dist;
            max_point.index = i;
        }
    }
    //printf("Index mais distante: %ld\n", max_point.index);
    //fflush(stdout);

    return max_point.index;
}

void printArray(double** vector, int n_dims, int start, int end)
{
    for(int i=start; i<=end; i++)
    {
        for(int j=0; j<n_dims; j++)
            printf("%lf ", vector[i][j]);
        
        printf("\n");
    }
    printf("\n");
}


/**
 * Computes the two points furthest apart in the set in two g
 * (might get local solution but won't affect the objective of the solution)
 * @param n_dims   # of dimensions
 * @param n_points # of points 
 * @param **pts array with the points
 * @param *idx_fp 
 */
void furthest_apart(int n_dims, long n_points, double **pts, long *idx_fp, MPI_Comm comm)
{
    int rank;
    long start;
    long end;
    long final_idx;
    int n_procs;
    int max_dist = 0;
    int save_idx=-1;
    double curr_dist = 0;
    max_n_idx max_point={-1.0,-1};
    double* recvValues;
    int winner = -1;
    //int winner2 = -1;
    long i=0;
    double* valueFromRoot = (double*)malloc(sizeof(double)*n_dims);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &n_procs);

    if(!rank)
    {
        recvValues = (double*)malloc(n_dims*n_procs*sizeof(double));
        arraycpy(valueFromRoot, pts[0], n_dims);
    }

    if(rank!=n_procs-1)
    {
        start = (int)(n_points/n_procs)*rank;
        end = (int)(n_points/n_procs)*(rank+1)-1;
    }
    else
    {
        start = (long)(n_points/n_procs)*rank;
        end = n_points-1;
    }
    
    //printf("Rank %d start %d, end %d, tem:\n", rank, start, end);
    //printArray(pts, n_dims, start, end);
    //fflush(stdout);
    
    MPI_Bcast(&(valueFromRoot[0]), n_dims, MPI_DOUBLE, 0, comm);

    save_idx = furthest_point_from_coords(n_dims, start, end, pts, valueFromRoot);

    MPI_Gather(&(pts[save_idx][0]), n_dims, MPI_DOUBLE, &(recvValues[0]), n_dims, MPI_DOUBLE, 0, comm);

    if(!rank)
    {
        max_dist = 0;
        winner = -1;

        for(i = 0; i < n_procs; i++)
        {
            //printf("Recebido %lf %lf\n", *(recvValues+i*n_dims), *(recvValues+i*n_dims+1));
            //fflush(stdout);
            curr_dist = squared_distance(n_dims, recvValues+i*n_dims, valueFromRoot);
            
            if(curr_dist > max_dist)
            {
                max_dist = curr_dist;
                winner = i;
            }
        }

        arraycpy(valueFromRoot, recvValues+winner*n_dims, n_dims);

        //printf("fromRoot antes:%lf %lf\n", valueFromRoot[0], valueFromRoot[1]);
        //fflush(stdout);
//
        //printf("winner1: %d\n", winner);
        //fflush(stdout);
    }

    MPI_Bcast(&winner, 1, MPI_INT, 0, comm);

    if(winner == rank)
    {
        //printf("A enviar: %d\n", save_idx);
        MPI_Send(&save_idx, 1, MPI_INT, 0, 0, comm);
    }

    //printf("rank %d, aqui1\n", rank);
    fflush(stdout);

    if(!rank)
    {
        MPI_Recv(&save_idx, 1, MPI_INT, winner, 0, comm, MPI_STATUS_IGNORE);
        //printf("Recebido index %d\n", save_idx);
        //fflush(stdout);
        idx_fp[0] = save_idx;
    }

    //printf("rank %d, aqui2\n", rank);
    //fflush(stdout);

    MPI_Bcast(&(valueFromRoot[0]), n_dims, MPI_DOUBLE, 0, comm);

    save_idx = furthest_point_from_coords(n_dims, start, end, pts, valueFromRoot);

    MPI_Gather(&(pts[save_idx][0]), n_dims, MPI_DOUBLE, &(recvValues[0]), n_dims, MPI_DOUBLE, 0, comm);

    //printf("rank %d, aqui3\n", rank);
    //fflush(stdout);

    if(!rank)
    {
        max_dist = 0;
        winner = -1;

        for (i = 0; i < n_procs; i++)
        {   
            //printf("Recebido %lf %lf\n", *(recvValues+i*n_dims), *(recvValues+i*n_dims+1));
            //fflush(stdout);

            curr_dist = squared_distance(n_dims, recvValues+i*n_dims, valueFromRoot);
            if(curr_dist > max_dist)
            {
                max_dist = curr_dist;
                winner = i;
            }
        }
        //printf("winner1: %d\n", winner);
        //fflush(stdout);
    }

    MPI_Bcast(&winner, 1, MPI_INT, 0, comm);

    //printf("rank %d, aqui4\n", rank);
    //fflush(stdout);

    if(winner == rank)
    {
        //printf("A enviar: %d\n", save_idx);
        //fflush(stdout);
        MPI_Send(&save_idx, 1, MPI_INT, 0, 0, comm);
    }

    //printf("rank %d, aqui5\n", rank);
    //fflush(stdout);

    if(!rank)
    {
        MPI_Recv(&save_idx, 1, MPI_INT, winner, 0, comm, MPI_STATUS_IGNORE);
        //printf("Recebido index %d\n", save_idx);
        idx_fp[1] = save_idx;
        //printf("Broadcast valores: %ld %ld\n",  idx_fp[0], idx_fp[1]);
        //fflush(stdout);
    }

    //printf("rank %d, aqui6\n", rank);
    //fflush(stdout);

    MPI_Bcast(&(idx_fp[0]), 2, MPI_LONG, 0, comm);

    //printf("rank %d, aqui7\n", rank);
    //fflush(stdout);

    //printf("rank %d - Valores: %ld %ld\n",rank,  idx_fp[0], idx_fp[1]);
    //fflush(stdout);

    //printf("----------------------------------\n");
    //fflush(stdout);

    //MPI_Barrier(comm);

    free(valueFromRoot);
    if(!rank)
        free(recvValues);

    return;
}


/**
 * Subtracts every dimensions of two points
 * @param   n_dims: # of dimensions
 * @param   a : first point
 * @param   b : second point
 * @param   result : array with the proper dimensions
 * @return  result
 */
double *subtraction(int n_dims, double *a, double *b, double* result)
{
    for (int i = 0; i < n_dims; i++)
        result[i] = b[i] - a[i];

    return result;
}


/**
 * Sums every dimensions of two points
 * @param   n_dims: # of dimensions
 * @param   a : first point
 * @param   b : second point
 * @param   result : array with the proper dimensions
 * @return  result
 */
double *sum(int n_dims, double *a, double *b, double *result)
{
    for (int i = 0; i < n_dims; i++)
        result[i] = a[i] + b[i];

    return result;
}


/**
 * Computes the inner product of two points
 * @param   n_dims: # of dimensions
 * @param   a : first point
 * @param   b : second point
 * @return  inner product
 */
double inner_product(int n_dims, double *a, double *b)
{
    double result = 0;
    for (int i = 0; i < n_dims; i++)
        result += a[i] * b[i];

    return result;
}


/**
 * Multiplies every dimensions of one points with a constant
 * @param   n_dims: # of dimensions
 * @param   a : first point
 * @param   constant : constant
 * @param   result : array with the proper dimensions
 * @param   result
 */
double *multiply(int n_dims, double *a, double constant, double* result)
{
    for (int i = 0; i < n_dims; i++)
        result[i] = constant * a[i];

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


/**
 * Given a list of points this function computes the reduced projections of every point and saves it on the projections array
 * @param   n_dims              # of dimensions
 * @param   projections         # of dimensions
 * @param   a                   first point of the line
 * @param   a                   second point of the line
 * @param   pts                 list of points
 * @param   n_points            # number of points
 * @return  result
 */
void project_pts2line(int n_dims, double* projections, double *a, double *b, double **pts, long n_points)
{
    double *b_minus_a = (double *)malloc(n_dims * sizeof(double));
    
    if(a[0]>b[0])
        subtraction(n_dims, b, a, b_minus_a);
    else
        subtraction(n_dims, a, b, b_minus_a);
    
    double* p_minus_a = (double *)malloc(n_dims * sizeof(double));

    for (int i = 0; i < n_points; i++)
        projections[i] = orthogonal_projection_reduced(n_dims, pts[i],a,p_minus_a,b_minus_a);

    free(p_minus_a);
    free(b_minus_a);
    
    return;
}