#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "gen_points.h"
#include "geometry.h"
#include "unsorted_median.h"
#include <string.h>

#define PRINT_NODE

int node_idx = 0;

/**
 * Prints recursively all the nodes
 * @param node_id : id of the node
 * @param tree : array where all the tree nodes are stored 
 * @param n_dims : # of dimensions
 */
void print_Node(node* foo, int n_dims)
{   
    if (foo->L != -1)
    {
        printf("%ld %ld %ld %.6lf", foo->id, foo->L, foo->R, foo->radius);

        for (int i = 0; i < n_dims; i++)
            printf(" %.6lf", foo->center[i]);

        printf(" \n");
    }
    else
    {
        printf("%ld -1 -1 0.000000", foo->id);

        for (int i = 0; i < n_dims; i++)
            printf(" %.6lf", foo->center[i]);

        printf(" \n");
    }
    fflush(stdout);
    return;
}

/**
 * Prints node vector
 * @param vector : Vector of pointers to nodes
 * @param n_dims : # of dimensions
 * @param n_nodes: Total number of nodes 
 */
void print_Nodes(node** vector, int n_dims, long n_nodes)
{
    for(int i=0; i<n_nodes; i++)
    {
        if(vector[i] != NULL)
        {
            print_Node(vector[i], n_dims);

            if(vector[i]->L != -1)
                free(vector[i]->center);
            
            free(vector[i]);
        }
    }
    return;
}

void build_tree(long node_id, double **pts, double* projections, long n_points, int n_dims, MPI_Comm comm, int rank, long start_npoints, int threads_available, node** node_dump)
{   
    long lnode_id;
    long rnode_id;
    long center_idx; //indice of the center of the pts array where the split for the childs is made
    long int fapart_idx = 0;

    node* foo = (node*)malloc(sizeof(node));

    foo->id=node_id;

    if(n_points == 1) //if the node is a leaf
    {
        foo->R = -1;
        foo->L = -1;
        foo->radius = 0;
        foo->center = pts[0];
        #ifdef PRINT_NODE
            if(rank==0)
                node_dump[node_id] = foo;
        #endif
        return;
    }

    lnode_id = node_id + 1; //indice of the left child
    foo->L=lnode_id;

    foo->center = (double *)malloc(n_dims * sizeof(double));

    if(n_points == 2) //only two points in the set -> easier/lesser operations
    {
        center_idx = 1;
        rnode_id = node_id + 2 * center_idx;
        foo->R = rnode_id;
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
            foo->center[i] = (pts[0][i] + pts[1][i]) / 2;
        }

        //both points are equidistant to the center
        foo->radius = distance(n_dims, pts[0], foo->center);
        
        //build leafs
        #ifdef PRINT_NODE
            if(rank==0)
                node_dump[node_id] = foo;
        #endif

        build_tree(lnode_id, pts, NULL, 1, n_dims,comm,rank,start_npoints,threads_available, node_dump); //center_idx happens to be the number of points in the set
        build_tree(rnode_id,pts+1, NULL, 1, n_dims,comm,rank,start_npoints,threads_available, node_dump);
        return;
    }
    else // > 2 points in the set
    {
        long median_idx; //indice of the median
        double median; //value of the median
        long idx_fp[2] = {0, 0}; // indices of the points furthest apart
        
        //compute furthest apart points in the current set
        furthest_apart(n_dims, n_points, pts, idx_fp, threads_available);
        //pseudo-projection of all points (enough to know relative positions)
        project_pts2line(n_dims, projections, pts[idx_fp[0]], pts[idx_fp[1]], pts, n_points, threads_available);

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
                foo->center[i] = (center_l[i] + center_r[i]) / 2;
            
            free(center_l);
            free(center_r);

            //place pts which projection is smaller than the median to left half of the array and greater to right half 
            compare_with_median(projections, pts, median, n_points, n_dims);
            center_idx = (n_points / 2);
        }
        else //odd n_pts -> median is the central value
        {
            //compute median
            median = getKsmallest(projections, n_points/2, n_points);

            median_idx = find_idx_from_value(projections, n_points, median);

            //compute and set center of the node
            orthogonal_projection(n_dims, pts[median_idx], pts[idx_fp[0]], pts[idx_fp[1]], foo->center);

            compare_with_median(projections, pts, median, n_points,n_dims);

            center_idx = (n_points - 1) / 2;
        }

        // finds the most distant to center, the distance between thems is the radius
        fapart_idx = furthest_point_from_coords(n_dims, n_points, pts, foo->center, threads_available);
    	foo->radius = distance(n_dims, pts[fapart_idx], foo->center);	

        rnode_id = node_id + 2 * center_idx;
        foo->R = rnode_id;
        
        #ifdef PRINT_NODE
            if(rank==0)
                node_dump[node_id] = foo;
        #endif
    }

    int size_world;
    MPI_Comm_size(comm, &size_world);
    
    if(size_world>=2)
    {
        MPI_Comm above_comm;
        MPI_Comm below_comm;
        //MPI_Comm new_comm;

        int above_median=(((double)rank/size_world)>=0.5)?1:0;
        int new_rank;
        if(above_median)
        {
            MPI_Comm_split(comm, 0, 0, &below_comm);
            MPI_Comm_rank(below_comm, &new_rank);
            build_tree(lnode_id, pts, projections, center_idx, n_dims, below_comm, new_rank, start_npoints, threads_available, node_dump); //center_idx happens to be the number of points in the set
        }
        else
        {
            MPI_Comm_split(comm, 1, 0, &above_comm);
            MPI_Comm_rank(above_comm, &new_rank);
            build_tree(rnode_id, pts + center_idx, projections + center_idx, n_points - center_idx, n_dims, above_comm, 
                       new_rank, start_npoints, threads_available, node_dump);
        }
    }
    else
    {
        if(threads_available == -1) //number of threads available is always one
        {
            build_tree(lnode_id,pts, projections, center_idx, n_dims,comm,rank,start_npoints,threads_available, node_dump); //center_idx happens to be the number of points in the set
            build_tree(rnode_id,pts + center_idx, projections + center_idx, n_points - center_idx, n_dims,comm,rank,start_npoints,threads_available, node_dump);
        }
        else if (threads_available == 1) //number of threads available in this recursion is down to one
        {
            build_tree(lnode_id,pts, projections, center_idx, n_dims,comm,rank,start_npoints,threads_available, node_dump); //center_idx happens to be the number of points in the set
            build_tree(rnode_id,pts + center_idx, projections + center_idx, n_points - center_idx, n_dims,comm,rank,start_npoints,threads_available, node_dump);
        }
        else if (node_idx != 0 && threads_available > 1)//still threads left to put to work
        {   
                #pragma omp task
                {
                    threads_available = threads_available/2;
                    omp_set_nested(2);
                    omp_set_num_threads(threads_available);
                    build_tree(lnode_id,pts, projections, center_idx, n_dims,comm,rank,start_npoints,threads_available, node_dump); //center_idx happens to be the number of points in the set
                }
                #pragma omp task
                {
                    threads_available = threads_available/2;
                    omp_set_nested(2);
                    omp_set_num_threads(threads_available);
                    build_tree(rnode_id,pts + center_idx, projections + center_idx, n_points - center_idx, n_dims,comm,rank,start_npoints,threads_available, node_dump);
                }
        }
        else if (node_idx == 0 && threads_available != -1)//multiple threads available and root node -> start parallel in the building of the tree
        {   
            node_idx = 1;
            #pragma omp parallel
            {
                #pragma omp single
                {
                    threads_available = threads_available/2;
                    #pragma omp task
                    {
                        omp_set_nested(2);
                        omp_set_num_threads(threads_available);
                        build_tree(lnode_id,pts, projections, center_idx, n_dims,comm,rank,start_npoints,threads_available, node_dump); //center_idx happens to be the number of points in the set
                    }
                    #pragma omp task
                    {
                        omp_set_nested(2);
                        omp_set_num_threads(threads_available);
                        build_tree(rnode_id,pts + center_idx, projections + center_idx, n_points - center_idx, n_dims,comm,rank,start_npoints,threads_available, node_dump);
                    }
                }
            }
        }
    }
    return;
}

int main(int argc, char *argv[])
{   
    MPI_Init(&argc, &argv); //START MPI
    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank); /*DETERMINE RANK OF THIS PROCESSOR*/
    MPI_Comm_size(MPI_COMM_WORLD, &size); /*DETERMINE TOTAL NUMBER OF PROCESSORS*/
    
    double exec_time;
    double starttime, endtime;
    double **pts;
    int n_dims = atoi(argv[1]); //number of dimensions
    long n_points = atoi(argv[2]); //number of points in the set
    long n_nodes = 2 * n_points - 1; //number of nodes in the tree

    node** node_dump = (node**)malloc(n_nodes*sizeof(node*));

    for(int i=0; i<n_nodes;i++)
    {
        node_dump[i] = NULL;
    }

    //____________START_TIME_BENCHMARK_____________ 
    exec_time = -omp_get_wtime();
    starttime = MPI_Wtime();

    //generates dataset
    pts = get_points(argc, argv);

    double* pts_first_position = pts[0]; // pointer to pt first pos for freeing purposes

    double* projections = (double*)malloc(n_points*sizeof(double)); //array to store pseudo-projections the locate the point in the line 

    #ifdef PRINT_NODE
    if(rank==0)
        printf("%d %ld\n", n_dims, n_nodes);
    #endif
    
    MPI_Comm comm = MPI_COMM_WORLD;
    int threads_available = omp_get_max_threads();
    
    build_tree(0, pts, projections, n_points, n_dims, comm ,rank, n_points, threads_available, node_dump);

    //____________END_TIME_BENCHMARK_____________
    MPI_Barrier(MPI_COMM_WORLD);

    endtime = MPI_Wtime();
    exec_time += omp_get_wtime();

    print_Nodes(node_dump, n_dims, n_nodes);
    free(node_dump);

    if(rank==0)
    {
        fprintf(stderr, "%.1lf\n", exec_time);
        //fprintf(stderr, " %.1lf seconds\n", endtime-starttime);
    }
    
    free(projections);
    free(pts_first_position);
    free(pts);

    MPI_Finalize();
    return 0;
}
