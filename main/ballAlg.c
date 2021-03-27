#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "gen_points.h"
#include "geometry.h"

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
    double *first_coord = (double *)malloc(n_points*sizeof(double));
    
    long center_idx;
    long fapart_idx;
    long idx_fp[2] = {0, 0}; //pontos a e b

    node *lnode = (node *)malloc(sizeof(node));
    if (lnode == NULL)
    {
        printf("Error allocating memory");
        exit(1);
    }

    lnode->id = newNode->id + 1;
    newNode->lnode = lnode;

    node *rnode = (node *)malloc(sizeof(node));
    if (rnode == NULL)
    {
        printf("Error allocating memory");
        exit(1);
    }

    
    //Calculate furthest points and make projections on their line
    recursive_furthest_apart(n_dims, n_points, pts, idx_fp);
    projections = project_pts2line(n_dims, pts[idx_fp[0]], pts[idx_fp[1]], pts, n_points);
    //Sort points and calculate Median
    //distances2a = calc_distances_to_left_limit(pts[idx_fp[0]], projections, n_points, n_dims);

    if(n_dims==1)
        quick_sort(pts, projections, (double *)pts, 0, n_points - 1);
    else
        get_first_coord(n_points, pts, first_coord);
        quick_sort(pts, projections, first_coord, 0, n_points - 1);

    newNode->center = (double *)malloc(n_dims * sizeof(double));
    center_idx = getMedian(projections, n_points, n_dims, newNode->center);
    //Calculate the radius
    fapart_idx = furthest_point_from_coords(n_dims, n_points, pts, newNode->center);
    newNode->radius = distance(n_dims, pts[fapart_idx], newNode->center);

    free(first_coord);
    for (int i = 0; i < n_points; i++)
    {
        free(projections[i]);
    }
    free(projections);

    build_tree(lnode, pts, center_idx, n_dims); //center_idx happens to be the number of points in the set
    rnode->id = newNode->id + 2 *  center_idx;
    newNode->rnode = rnode;
    build_tree(rnode, pts + center_idx, n_points - center_idx, n_dims);
 

    return;
}

void print_Node(node foo, int n_dims)
{
    //node_id left_child_id right_child_id radius center_coordinates
    if (foo.lnode != NULL)
    {
        printf("%ld %ld %ld %.6lf", foo.id, foo.lnode->id, foo.rnode->id, foo.radius);

        for (int i = 0; i < n_dims; i++)
            printf(" %.6lf", foo.center[i]);

        printf("\n");

        print_Node(*(foo.lnode), n_dims);

        print_Node(*(foo.rnode), n_dims);
    }
    else
    {
        printf("%ld -1 -1 0.000000", foo.id);

        for (int i = 0; i < n_dims; i++)
            printf(" %.6lf", foo.center[i]);

        printf("\n");
    }
    return;
}

void dump_tree(node tree_root, int n_dims, long n_points)
{
    long n_nodes = 2 * n_points - 1;
    printf("%d %ld\n", n_dims, n_nodes);
    print_Node(tree_root, n_dims);
}

int main(int argc, char *argv[])
{
    double exec_time;
    double **pts;
    node tree_root;
    tree_root.id = 0;
    long last_id;

    exec_time = -omp_get_wtime();
    pts = get_points(argc, argv);
    int n_dims = atoi(argv[1]);
    
    //printf("DIMS:%d\n", n_dims);
    long n_points = atoi(argv[2]);
    /*printf("NPOINTS:%ld\n", n_points);
    for (long i = 0; i < n_points; i++)
    {
        printf("POINT %ld\n", i);
        for (int j = 0; j < n_dims; j++)
        {
            printf("%f ", pts[i][j]);
        }
        printf("\n");
    }
    */
    build_tree(&tree_root, pts, n_points, n_dims);
    

    exec_time += omp_get_wtime();
    printf("%.1lf\n", exec_time);

    dump_tree(tree_root, n_dims, n_points);
}
