#include <omp.h>
#include <stdio.h>
#include "gen_points.h"
#include "distance.h"

void flag(int n)
{
    printf("Flag %d\n", n);
    fflush(stdout);
}

node build_tree(double **pts, int npoints, int dimensions)
{
    printf("Building Tree...\n");

    node root;

    double **limits; //pontos a e b

    limits = furthest_apart(npoints, dimensions, pts, 0, npoints-1);

    
    
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
