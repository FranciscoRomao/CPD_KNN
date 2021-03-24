#ifndef _DISTANCE_H_
#define _DISTANCE_H_

typedef struct node
{
    long id;
    struct node *lnode;
    struct node *rnode;
    double radius;
    double *center;
} node;

void swap(double **pts, double **projections, double *distances2a, long a, long b);
long partition(double **pts, double **projections, double *distances2a, long low, long high);
void quick_sort(double **pts, double **projections, double *distances2a, long low, long high);
long getMedian(double **projections, long n_points, int n_dims, double *center);
double *calc_distances_to_left_limit(double *left_limmit, double **projections, long n_points, int n_dims);
double distance(int n_dims, double *a, double *b);
long furthest_point_from_coords(int n_dims, long n_points, double **pts, double *base_coords);
void recursive_furthest_apart(int n_dims, long n_points, double **pts, long *idx_fp);
double *subtraction(int n_dims, double *a, double *b);
double *sum(int n_dims, double *a, double *b);
double inner_product(int n_dims, double *a, double *b);
double *multiply(int n_dims, double *a, double constant);
double *orthogonal_projection(int n_dims, double *p, double *a, double *b);
double **project_pts2line(int n_dims, double *a, double *b, double **pts, long n_points);
void dump_tree(node tree_root, int n_dims, long n_points);
void flag(int n); //Tirar no fim
#endif