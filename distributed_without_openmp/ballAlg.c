#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "gen_points.h"
#include "geometry.h"
#include "unsorted_median.h"

void build_tree(node* tree, long node_idx, double **pts, double* projections, long n_points, int n_dims)
{
    if (n_points == 1) //if the node is a leaf
    {
        tree[node_idx].R = -1;
        tree[node_idx].L = -1;
        tree[node_idx].radius = 0;
        tree[node_idx].center = pts[0];
        return;
    }

    long center_idx; //indice of the center of the pts array where the split for the childs is made
    long lnode_id = node_idx + 1; //indice of the left child

    tree[node_idx].L=lnode_id;

    tree[node_idx].center = (double *)malloc(n_dims * sizeof(double));

    if(n_points == 2) //only two points in the set -> easier/lesser operations
    {
        double* aux;
        if (pts[0][0] > pts[1][0]) //sort ascending if not already
        {
            aux = pts[0];
            pts[0] = pts[1];
            pts[1] = aux;
        }
        //median of two points is its average
        for(int i=0; i<n_dims; i++)
        {
            tree[node_idx].center[i] = (pts[0][i] + pts[1][i]) / 2;
        }
        center_idx = 1;
        //both points are equidistant to the center
        tree[node_idx].radius = distance(n_dims, pts[0], tree[node_idx].center);
        //build leafs
        build_tree(tree, lnode_id, pts, NULL, 1, n_dims); //center_idx happens to be the number of points in the set
        long rnode_id = node_idx + 2 * center_idx;
        tree[node_idx].R = rnode_id;
        build_tree(tree, rnode_id, pts+1, NULL, 1, n_dims);

        return;
    }
    else // > 2 points in the set
    {
        long median_idx; //indice of the median
        double median; //value of the median
        long fapart_idx; // indice of the point furthest apart from the center
        long idx_fp[2] = {0, 0}; // indices of the points furthest apart
        double radius_candidate[2] = {0, 0 }; //possible radius
        
        //compute furthest apart points in the current set    
        /** MPI_plan - master machine sends the coordinate of a point to all machines
         *           - each machine computes the furthest point from and returns to master
         *           - master says which of this points is the most further (pt a) and sends it to all machines
         *           - each machine computes the furthest point from pt a and rets to master
         *           - master who syas which is the furthest from a, being pt b
        */  
        recursive_furthest_apart(n_dims, n_points, pts, idx_fp);
        //pseudo-projection of all points (enough to know relative positions) 
        /** MPI_plan - furthest pts are sent to all machines and the project the points 
         * 
        */
        project_pts2line(n_dims, projections, pts[idx_fp[0]], pts[idx_fp[1]], pts, n_points);

        if(n_points % 2 == 0) //even n_pts -> median is the avergae of 2 central values
        {
            long median_left_idx; //indice of the immediatly smaller value than the median
            long median_right_idx; //indice of the immediatly bigger value than the median
            double median_left; //immediatly smaller value than the median
            double median_right; //immediatly bigger value than the median

            //compute median as the avergae of the two central pts
            median_right = getKsmallest(projections, n_points/2 , n_points);
            median_right_idx = find_idx_from_value(projections, n_points, median_right);
            median_left_idx = getLowerNeighborIdx(projections, n_points, median_right);
            median_left = projections[median_left_idx];
            median = (median_left + median_right) / 2;
            //compute and set center of the node
            double* center_l=(double*)malloc(n_dims*sizeof(double));
            double* center_r=(double*)malloc(n_dims*sizeof(double));
            orthogonal_projection(n_dims, pts[median_left_idx], pts[idx_fp[0]], pts[idx_fp[1]], center_l);
            orthogonal_projection(n_dims, pts[median_right_idx], pts[idx_fp[0]], pts[idx_fp[1]], center_r);
            for(int i=0; i<n_dims; i++)
            {
                tree[node_idx].center[i] = (center_l[i] + center_r[i]) / 2;
            }
            free(center_l);
            free(center_r);           
            //place pts which projection is smaller than the median to left half of the array and greater to right half 
            compare_with_median(projections, pts, median, n_points);
            center_idx = (n_points / 2);

        }
        else //odd n_pts -> median is the central value
        {   
            //compute median
            median = getKsmallest(projections, n_points/2, n_points);
            median_idx = find_idx_from_value(projections, n_points, median);
            //compute and set center of the node
            orthogonal_projection(n_dims, pts[median_idx], pts[idx_fp[0]], pts[idx_fp[1]], tree[node_idx].center);
            compare_with_median(projections, pts, median, n_points);
            center_idx = (n_points - 1) / 2;
        }
/*
        //compute the furthest point of the ball center and thus the radius
        radius_candidate[0] = distance(n_dims, pts[idx_fp[0]], tree[node_idx].center);
        radius_candidate[1] = distance(n_dims, pts[idx_fp[1]], tree[node_idx].center);

        if (radius_candidate[1]>radius_candidate[0]) 
        {
            tree[node_idx].radius = radius_candidate[1];
        }
        else 
        {
            tree[node_idx].radius = radius_candidate[0];
        }            
*/

	fapart_idx = furthest_point_from_coords(n_dims, n_points, pts, tree[node_idx].center);
    	tree[node_idx].radius = distance(n_dims, pts[fapart_idx], tree[node_idx].center);	
        // build child
        build_tree(tree, lnode_id, pts, projections, center_idx, n_dims); //center_idx happens to be the number of points in the set
        long rnode_id = node_idx + 2 * center_idx;
        tree[node_idx].R = rnode_id;
        build_tree(tree, rnode_id, pts + center_idx, projections + center_idx, n_points - center_idx, n_dims);

        return;
    } 
}

/**
 * Prints recursively all the nodes
 * @param node_id : id of the node
 * @param tree : array where all the tree nodes are stored 
 * @param n_dims : # of dimensions
 */
void print_Node(long node_id ,node* tree, int n_dims)
{   
    node foo = tree[node_id];
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

/**
 * Prints the tree to the standard output
 * @param n_dims: # of dimensions
 * @param n_points : # of points
 * @param n_nodes : # of nodes
 */
void dump_tree(node* tree, int n_dims, long n_points,long n_nodes)
{
    printf("%d %ld\n", n_dims, n_nodes);
    print_Node(0,tree, n_dims);
}

/**
 * Free tree memory
 * @param tree : array where all the tree nodes are stored 
 * @param n_nodes : # of nodes
 */
void destroy_tree(long n_nodes, node* tree)
{
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
    long n_local_points;

    //____________START_TIME_BENCHMARK_____________ 
    //----------N??o sei o que fazer com isto
    exec_time = -omp_get_wtime(); 
    MPI_Init(&argc, &argv);

    pts = get_points_mpi(argc, argv, MPI_COMM_WORLD,&n_local_points);
    
        double* pts_first_position = pts[0]; //position of the first element of pts (pts is sorted so it would be lost and hard to free)
        int n_dims = atoi(argv[1]); //number of dimensions 
        long n_points = atoi(argv[2]); //number of points in the set
        double* projections = (double*)malloc(n_local_points*sizeof(double)); //array to store pseudo-projections the locate the point in the line 
        long n_nodes = 2 * n_points - 1; //number of nodes in the tree
        long n_local_nodes = 2*n_points-1;
        tree = (node*)malloc(n_nodes*sizeof(node));

        // construct THE TREE

        //Odd_ranks={1, 3, 5}, Even_ranks={0, 2, 4}
        /*
        MPI_Comm_group(MPI_COMM_WORLD, Old_group)
        MPI_Group_incl(Old_group, 3, Odd_ranks, &Odd_group)
        MPI_Group_incl(Old_group, 3, Even_ranks, &Even_group)
        MPI_Comm_create(MPI_COMM_WORLD, Odd_group, Odd_Comm)
        MPI_Comm_create(MPI_COMM_WORLD, Even_group, Even_Comm) 

        build_tree(tree, 0, pts, projections, n_points, n_dims, Odd_Comm); 10
        build_tree(tree, 0, pts, projections, n_points, n_dims, Even_Comm); 10
        */
    

    MPI_Finalize();
    //____________END_TIME_BENCHMARK_____________
    /*
    exec_time += omp_get_wtime();

    dump_tree(tree, n_dims, n_points,n_nodes);
    destroy_tree(n_nodes,tree);
    */
    free(projections);
    free(pts_first_position);
    free(pts);
    //fprintf(stderr, "%.1lf\n", exec_time);
}
