#include <omp.h>
#include <stdio.h>
#include "gen_points.h"
#include "geometry.h"

void flag(int n)
{
    printf("Flag %d\n", n);
    fflush(stdout);
}

node build_tree(double **pts, int npoints, int ndims)
{
    printf("Building Tree...\n");

    node root;

    double **limits; //pontos a e b
    double **projections; //projections from all points on the line ab (including a and b)
    double radius;

    limits = furthest_apart(npoints, ndims, pts, 0, npoints-1);

    projections = project_on_ab(limits, pts, npoints, ndims);

    radius = calc_Radius(projections, npoints, ndims);
    
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

    root = build_tree();

    exec_time += omp_get_wtime();

    //fprintf(stderr, "%.1lf\n", exec_time);
    printf("%.1lf\n", exec_time);

    dump_tree(&root);
}
