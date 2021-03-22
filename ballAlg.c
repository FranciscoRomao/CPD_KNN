#include <omp.h>
#include <stdio.h>
#include "gen_points.h"
#include "geometry.h"

void flag(int n)
{
    printf("Flag %d\n", n);
    fflush(stdout);
}

build_tree(
    count++ if (size_array == 1)
        root->left = -1 right..raio = 0 center = own_point return;
    root->median = ....root->dasd..alocas(root->left) metes count++ alocas(root->right) metes count++ build_tree(root->left, )
                       build_tree(roof->right, ) return )

void build_tree(node *newNode, double **pts, int npoints, int n_dims, long start, long end)
{

    int *limits;          //pontos a e b
    double **projections; //projections from all points on the line ab (including a and b)
    double radius;
    double *center;
    long idx2project[1] = -1; //repensar
    double *distances2a;
    int aux_int;
    double aux1;
    double aux2;
    node lnode;
    node rnode;

    limits = furthest_apart(npoints, n_dims, pts, 0, npoints - 1);

    projections = project_pts2line(n_dims, limits[0], limits[1], pts, idx2project, npoints);

    distances2a = calc_distances_to_left_limit(limits[0], projections, npoints, n_dims);

    quick_sort(pts, projections, distances2a);

    center_idx = getMedian(projections, npoints, n_dims, newNode->center);

    newNode->radius = 0;

    aux = furthest_point(n_dims, np, pts, newNode->center); //Talvez temos de criar uma nova função que
                                                                             //recebe as coordenadas de um ponto em vez de indice

    free(distances2a);

    for(int i=0; i<npoints; i++)
    {
        free(projections[i]);
    }    
    free(projections);


    if (aux1 > aux2)
    {
        radius = aux1;
    }
    else
    {
        radius = aux2;
    }

    if(center_idx == 0) //Se já não existir pontos à esquerda
    {
        newNode->lid = -1;
        newNode->lnode = NULL;
    }
    else //Se ainda houver pontos à esquerda
    {
        node lnode;
        lnode->id = (newNode->id) * 2 + 1;
        build_tree(lnode, pts, center_idx, n_dims);//center_idx happens to be the number of points in the set
    }

    if(npoints - center_idx == 1) //Significa que já não existe pontos à direita
    {
        newNode->rid = -1;
        newNode->rnode = NULL;
    }
    else
    {
        node rnode;
        rnode->id = (newNode->id) * 2 + 2;
        build_tree(rnode, pts+center_idx, npoints-center_idx, n_dims);
    }

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

    exec_time = -omp_get_wtime();

    pts = get_points(argc, argv);

    root.id = 0;
    root.lnode = NULL;
    root.rnode = NULL;
    build_tree(&root, pts, npoints, ndims, -1);

    exec_time += omp_get_wtime();

    //fprintf(stderr, "%.1lf\n", exec_time);
    printf("%.1lf\n", exec_time);

    dump_tree(&root);
}
