#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "gen_points.h"
#include "geometry.h"

void flag(int n)
{
    printf("Flag %d\n", n);
    fflush(stdout);
}

void build_tree(node *newNode, double **pts, long n_points, int n_dims)
{
    if (n_points == 1)
    {
        newNode->rnode = NULL;
        newNode->lnode = NULL;
        newNode->radius = 0;
        newNode->center = pts[0];
        return;
    }

    double **projections; //projections from all points on the line ab (including a and b)
    double radius;
    double *center;
    double *distances2a;
    long center_idx;
    long fapart_idx;
    long idx_fp[2] = {0, 0}; //pontos a e b

    node *lnode = (node *)malloc(sizeof(node));
    lnode->id = 2 * newNode->id + 1;

    node *rnode = (node *)malloc(sizeof(node));
    rnode->id = 2 * newNode->id + 2;

    //Calculate furthest points and make projections on their line
    recursive_furthest_apart(n_dims, n_points, pts, idx_fp);
    projections = project_pts2line(n_dims, pts[idx_fp[0]], pts[idx_fp[1]], pts, n_points);
    //Sort points and calculate Median
    distances2a = calc_distances_to_left_limit(pts[idx_fp[0]], projections, n_points, n_dims);
    quick_sort(pts, projections, distances2a, 0, n_points);
    center_idx = getMedian(projections, n_points, n_dims, newNode->center);
    //Calculate the radius
    fapart_idx = furthest_point_from_coords(n_dims, n_points, pts, newNode->center);
    newNode->radius = distance(n_dims, pts[fapart_idx], newNode->center);

    free(distances2a);
    for (int i = 0; i < n_points; i++)
    {
        free(projections[i]);
    }
    free(projections);

    build_tree(lnode, pts, center_idx, n_dims); //center_idx happens to be the number of points in the set
    build_tree(rnode, pts + center_idx, n_points - center_idx, n_dims);

    return;
}

void dump_tree(node *root)
{
}

int main(int argc, char *argv[])
{
    double exec_time;
    double **pts;
    node root;
    long last_id;
    /*
    exec_time = -omp_get_wtime();

    pts = get_points(argc, argv);

    root.id = 0;
    root.lnode = NULL;
    root.rnode = NULL;
    build_tree(&root, pts, n_points, n_dims, -1);

    exec_time += omp_get_wtime();

    //fprintf(stderr, "%.1lf\n", exec_time);
    printf("%.1lf\n", exec_time);

    dump_tree(&root);
    */
}
