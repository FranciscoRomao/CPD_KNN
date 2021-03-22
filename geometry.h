#ifndef _DISTANCE_H_
#define _DISTANCE_H_

typedef struct
{
    long id;
    long idx; //será que é necessário??
    node *lnode;
    node *rnode;
    int radius;
    double *center;
} node;

double distance(int dimensions, double *x_coords, double *y_coords);                                  //compute the distance between two points
int *furthest_apart(int numb_points, int dimensions, point *pts, int initial_index, int start_index); //comput the points that are furthest apart from each other
double **project_on_ab(double **limits_ab, double **pts);
double *subtraction(int n_dims, double *a, double *b);
double *sum(int n_dims, double *a, double *b);
double *calc_distances2a(pts);
double *multiply(int n_dims, double *a, int constant);
double inner_product(int n_dims, double *a, double *b);
long furthest_point(int n_dims, long np, double **pts, long base_idx);
double *orthogonal_projection(int n_dims, double *p, double *a, double *b);
double **project_pts2line(int n_dims, double *a, double *b, double **pts, long *idx2project, long np2project);

#endif