#include <omp.h>
#include <stdio.h>
#include "gen_points.h"
#include "geometry.h"

void flag(int n)
{
    printf("Flag %d\n", n);
    fflush(stdout);
}

node build_tree(node *newNode, double **pts, int npoints, int n_dims, int last_id)
{
    printf("Building Tree...\n");

    newNode = (*node)malloc(sizeof(node));

    int *limits;          //pontos a e b
    double **projections; //projections from all points on the line ab (including a and b)
    double radius;
    double *center;
    long idx2project[1] = -1; //repensar
    double *distances2a;
    int aux_int;
    double aux_dbl1;
    double aux_dbl2;
    node *lnode;
    node *rnode;

    limits = furthest_apart(npoints, n_dims, pts, 0, npoints - 1);

    projections = project_pts2line(n_dims, limits[0], limits[1], pts, idx2project, npoints);

    distances2a = calc_distances_to_left_limit(limits[0], projections, npoints, n_dims);

    quick_sort(pts, projections, distances2a);

    center_idx = getMedian(projections, npoints, n_dims, newNode->center);

    aux_dbl1 = distance(n_dims, newNode->center, distances2a[0]);
    aux_dbl2 = distance(n_dims, newNode->center, distances2a[npoints - 1]);

    if (aux_dbl1 > aux_dbl2)
    {
        radius = aux_dbl1;
    }
    else
    {
        radius = aux_dbl2;
    }

    if (center_idx == 1)
    {
        newNode->lnode
    }

    if ()

        build_tree()

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
    node *root;
    long last_id;

    exec_time = -omp_get_wtime();

    pts = get_points(argc, argv);

    root = build_tree(&root, double **pts, int npoints, int ndims, -1);

    exec_time += omp_get_wtime();

    //fprintf(stderr, "%.1lf\n", exec_time);
    printf("%.1lf\n", exec_time);

    dump_tree(&root);
}
