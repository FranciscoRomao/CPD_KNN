#ifndef _DISTANCE_H_
#define _DISTANCE_H_

#include <mpi.h>
MPI_Datatype MAXNIDX;
MPI_Op REDUCEMAXOP;

typedef struct node {
    int n_nodes;
    double radius;
    double* center;
    long L;
    long R;
    long id;
    struct node* nextNode;
} node;


typedef struct max_n_idx
{
    double maximum;
    long int index; 
}max_n_idx;

max_n_idx max_n_idx_max(max_n_idx a, max_n_idx b)
{
    return a.maximum > b.maximum ? a : b; 
}

void reducemaxop(MAXNIDX *in, MAXNIDX *inout, int *len, MPI_Datatype *datatype)
{   
    if(in->maximum > inout->maximum)
    {
        inout->maximum=in->maximum;
        inout->index = inout->index;
    }
}

double distance(int n_dims, double *a, double *b);
double squared_distance(int n_dims, double *a, double *b);
long furthest_point_from_coords(int n_dims, long n_points, double **pts, double *base_coords, int size_world, MPI_Comm comm, int rank);
void furthest_apart(int n_dims, long n_points, double **pts, long *idx_fp, int size_world, MPI_Comm comm ,int rank);
double *subtraction(int n_dims, double *a, double *b, double* result);
double *sum(int n_dims, double *a, double *b, double * result);
double inner_product(int n_dims, double *a, double *b);
double *multiply(int n_dims, double *a, double constant, double* result);
double orthogonal_projection_reduced(int n_dims, double *p, double *a, double* p_minus_a, double* b_minus_a);
void orthogonal_projection(int n_dims, double *p, double *a, double *b,double* result_sum);
void project_pts2line(int n_dims, double* projections, double *a, double *b, double **pts, long n_points);


#endif