#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "gen_points.h"
#include "geometry.h"

#define OMP_NUM_THREADS 4

void build_tree(node* tree, long node_idx, double **pts, double* projections, long n_points, int n_dims)
{   
    if (n_points == 1)
    {
        tree[node_idx].R = -1;
        tree[node_idx].L = -1;
        tree[node_idx].radius = 0;
        tree[node_idx].center = pts[0];
        return;
    }

    //double radius;
    //double *center;
    
    long center_idx;
    long fapart_idx;
    long idx_fp[2] = {0, 0}; //pontos a e b

    long lnode_id = node_idx + 1;
    tree[node_idx].L=lnode_id;

    //Calculate furthest points and make projections on their line
    recursive_furthest_apart(n_dims, n_points, pts, idx_fp);
    project_pts2line(n_dims, projections, pts[idx_fp[0]], pts[idx_fp[1]], pts, n_points);
    //Sort points and calculate Median
    //distances2a = calc_distances_to_left_limit(pts[idx_fp[0]], projections, n_points, n_dims);
    quick_sort(pts, projections, 0, n_points - 1);
    tree[node_idx].center = (double *)malloc(n_dims * sizeof(double));
    center_idx = getMedian(pts, n_points, n_dims, tree[node_idx].center);
    
    //Calculate the radius
    fapart_idx = furthest_point_from_coords(n_dims, n_points, pts, tree[node_idx].center);
    tree[node_idx].radius = distance(n_dims, pts[fapart_idx], tree[node_idx].center);

    build_tree(tree,lnode_id, pts, projections, center_idx, n_dims); //center_idx happens to be the number of points in the set
    long rnode_id = node_idx + 2 *  center_idx;
    tree[node_idx].R = rnode_id;
    build_tree(tree,rnode_id, pts + center_idx, projections+center_idx, n_points - center_idx, n_dims);

    return;
}

void print_Node(long node_id ,node* tree, int n_dims)
{   
    node foo=tree[node_id];
    //node_id left_child_id right_child_id radius center_coordinates
    if (foo.L != -1)
    {
        printf("%ld %ld %ld %.6lf", node_id, foo.L, foo.R, foo.radius);

        for (int i = 0; i < n_dims; i++)
            printf(" %.6lf", foo.center[i]);

        printf(" \n");

        print_Node(foo.L,tree, n_dims);

        print_Node(foo.R,tree, n_dims);
    }
    else
    {
        printf("%ld -1 -1 0.000000", node_id);

        for (int i = 0; i < n_dims; i++)
            printf(" %.6lf", foo.center[i]);

        printf(" \n");
    }
    return;
}

void dump_tree(node* tree, int n_dims, long n_points,long n_nodes)
{
    printf("%d %ld\n", n_dims, n_nodes);
    print_Node(0,tree, n_dims);
}

void destroy_tree(long n_nodes, node* tree){
    for(long i=0;i<n_nodes;i++)
    {   
        if(tree[i].L!=-1)
            free(tree[i].center);
    }
    free(tree);
}

int main(int argc, char *argv[])
{
    double exec_time;
    double **pts;
    node* tree;
    //long last_id;

    exec_time = -omp_get_wtime();
    pts = get_points(argc, argv);
    double* pts_first_position=pts[0];
    
    int n_dims = atoi(argv[1]);
    long n_points = atoi(argv[2]);
    long n_nodes = 2 * n_points - 1;

    double* projections=(double*)malloc(n_points*sizeof(double));
    tree=(node*)malloc(n_nodes*sizeof(node));
    build_tree(tree, 0, pts, projections, n_points, n_dims);
  
    exec_time += omp_get_wtime();

    //dump_tree(tree, n_dims, n_points,n_nodes);
    destroy_tree(n_nodes,tree);
    free(projections);
    free(pts_first_position);
    free(pts);
    fprintf(stderr, "%.1lf\n", exec_time); 
}
