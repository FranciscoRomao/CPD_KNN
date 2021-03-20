#include <omp.h>
#include <stdio.h>
#include "gen_points.h"
#include "geometry.h"

void flag(int n)
{
    printf("Flag %d\n", n);
    fflush(stdout);
}

node build_tree(node *newNode, double **pts, int npoints, int ndims)
{
    printf("Building Tree...\n");

    int *limits;          //pontos a e b
    double **projections; //projections from all points on the line ab (including a and b)
    double radius;
    double *center;
    long idx2project[1] = -1; //repensar
    double *distances2a;

    limits = furthest_apart(npoints, ndims, pts, 0, npoints - 1);

    projections = project_pts2line(ndims, limits[0], limits[1], pts, idx2project, npoints);

    distances2a = calc_distances_to_left_limit(limits[0], projections, npoints, ndims);

    center = getCenter(projections);

    newNode->center = center;

    radius = distance(ndims, center, limits[0]);

    newNode->center = center;

    set_splitter(radius, distances2a);

    //....

    printf("Done\n");

    return root;
}

void dump_tree(node *root)
{
    printf("Tree pro lixo\n");
}

int main(int argc, char *argv[])
{
    double exec_time;
    double **pts;
    node root;

    exec_time = -omp_get_wtime();

    pts = get_points(argc, argv);

    root = build_tree(&root, double **pts, int npoints, int ndims);

    exec_time += omp_get_wtime();

    //fprintf(stderr, "%.1lf\n", exec_time);
    printf("%.1lf\n", exec_time);

    dump_tree(&root);
}
