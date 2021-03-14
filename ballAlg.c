#include <omp.h>
#include <stdio.h>
#include "gen_points.h"

typedef struct
{
    int id;
    int lid;
    int rid;
    int radius;
    double *center;
} node;

void flag(int n)
{
    printf("Flag %d\n", n);
    fflush(stdout);
}

node build_tree()
{
    node root;

    printf("Building Tree...\n");
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
