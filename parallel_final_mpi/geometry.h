#ifndef _DISTANCE_H_
#define _DISTANCE_H_

typedef struct node {
    double radius;
    double * center;
    long L;
    long R;
    long id;
} node;

double distance(int n_dims, double *a, double *b);
double squared_distance(int n_dims, double *a, double *b);
long furthest_point_from_coords(int n_dims, long n_points, double **pts, double *base_coords, int threads_available);
void furthest_apart(int n_dims, long n_points, double **pts, long *idx_fp, int threads_available);
double *subtraction(int n_dims, double *a, double *b, double* result);
double *sum(int n_dims, double *a, double *b, double * result);
double inner_product(int n_dims, double *a, double *b);
double *multiply(int n_dims, double *a, double constant, double* result);
double orthogonal_projection_reduced(int n_dims, double *p, double *a, double* p_minus_a, double* b_minus_a);
void orthogonal_projection(int n_dims, double *p, double *a, double *b,double* result_sum);
void project_pts2line(int n_dims, double* projections, double *a, double *b, double **pts, long n_points, int threads_available);


#endif