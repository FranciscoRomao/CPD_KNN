#ifndef _DISTANCE_H_
#define _DISTANCE_H_

typedef struct
{
    int id;
    int lid;
    int rid;
    int radius;
    double *center;
} node;


double distance(int dimensions, double* x_coords, double* y_coords);//compute the distance between two points
int* furthest_apart(int numb_points,int dimensions, point* pts, int initial_index, int start_index);//comput the points that are furthest apart from each other
double **project_on_ab(double **limits_ab, double **pts);
double *calc_distances2a(pts);
long furthest_point(int n_dims, long np, double **pts, long base_idx);
double calc_Radius(double *distances, int npoints);

#endif