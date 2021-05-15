#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "gen_points.h"
#include "geometry.h"
#include "unsorted_median.h"
#include <string.h>

MPI_Datatype mpi_node_type;

void flag (int x)
{   
    printf("%d\n", x);
    fflush(stdout);
}

void print_proc_name(char * name){
    printf("%s\n",name);
    fflush(stdout);
}

void my_print(char* name, int x){
    if(strcmp(name,"lab13p3")==0){
        printf("%d\n",x);
        fflush(stdout);
    }
}

/**
 * Prints recursively all the nodes
 * @param node_id : id of the node
 * @param tree : array where all the tree nodes are stored 
 * @param n_dims : # of dimensions
 */
void print_Node(node foo, int n_dims)
{   
    //node_id left_child_id right_child_id radius center_coordinates
    if (foo.L != -1)
    {
        printf("%ld %ld %ld %.6lf", foo.id, foo.L, foo.R, foo.radius);

        for (int i = 0; i < n_dims; i++)
            printf(" %.6lf", foo.center[i]);

        printf(" \n");
    }
    else
    {
        printf("%ld -1 -1 0.000000", foo.id);

        for (int i = 0; i < n_dims; i++)
            printf(" %.6lf", foo.center[i]);

        printf(" \n");
    }
    return;
}


void build_tree(node* tree, long node_id, long local_start_id, double **pts, double* projections, long n_points, int n_dims, MPI_Comm comm, int rank, long start_npoints, long start_ndims)
{   
    long lnode_id;
    long rnode_id;
    long center_idx; //indice of the center of the pts array where the split for the childs is made
    MPI_Comm above_comm;
    MPI_Comm below_comm;
    char my_name[MPI_MAX_PROCESSOR_NAME];
    int my_len;
    MPI_Get_processor_name(my_name, &my_len);
    //print_proc_name(my_name);
    if (rank!=0 && n_points == 1) //if the node is a leaf
        return;
    if(rank!=0 && n_points==2)
        return;

    if(rank==0){
        long node_idx=node_id;
        if (n_points == 1) //if the node is a leaf
        {
            tree[node_idx].R = -1;
            tree[node_idx].L = -1;
            tree[node_idx].radius = 0;
            tree[node_idx].center = pts[0];
            tree[node_idx].id=node_id;
            print_Node(tree[node_idx],n_dims);
            return;
        }

        tree[node_idx].id=node_id;
        lnode_id = node_id + 1; //indice of the left child

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
            print_Node(tree[node_idx],n_dims);
            build_tree(tree, lnode_id,local_start_id, pts, NULL, 1, n_dims,comm,rank,start_npoints,start_ndims); //center_idx happens to be the number of points in the set
            rnode_id = node_id + 2 * center_idx;
            tree[node_idx].R = rnode_id;
            build_tree(tree, rnode_id,local_start_id,pts+1, NULL, 1, n_dims,comm,rank,start_npoints,start_ndims);

            return;
        }
        else // > 2 points in the set
        {
            long median_idx; //indice of the median
            double median; //value of the median
            long idx_fp[2] = {0, 0}; // indices of the points furthest apart
            double radius_candidate[2] = {0, 0}; //possible radius
            
            //compute furthest apart points in the current set
            recursive_furthest_apart(n_dims, n_points, pts, idx_fp); 
            //pseudo-projection of all points (enough to know relative positions) 
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
                //printf("foo1\n");
                fflush(stdout);
                median_idx = find_idx_from_value(projections, n_points, median);
                //printf("foo2\n");
                fflush(stdout);
                //compute and set center of the node
                orthogonal_projection(n_dims, pts[median_idx], pts[idx_fp[0]], pts[idx_fp[1]], tree[node_idx].center);
                //printf("foo3\n");
                fflush(stdout);
                compare_with_median(projections, pts, median, n_points);
                //printf("foo4\n");
                fflush(stdout);
                center_idx = (n_points - 1) / 2;
            }

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
            rnode_id = node_id + 2 * center_idx;
            tree[node_idx].R = rnode_id;
            print_Node(tree[node_idx],n_dims);
        }
    }

    int size_world;
    MPI_Comm_size(MPI_COMM_WORLD, &size_world);
    if(size_world>=2)
    {   
        MPI_Bcast(pts,start_npoints*start_npoints,MPI_DOUBLE,0,comm);
        MPI_Bcast(projections,start_npoints,MPI_DOUBLE,0,comm);
        MPI_Bcast(tree,2*start_npoints-1,mpi_node_type,0,comm);
        MPI_Bcast(&rnode_id,1,MPI_LONG,0,comm);
        MPI_Bcast(&lnode_id,1,MPI_LONG,0,comm);
        MPI_Bcast(&center_idx,1,MPI_LONG,0,comm);
        MPI_Bcast(&n_points,1,MPI_LONG,0,comm);
        MPI_Comm above_comm;
        MPI_Comm below_comm;
        int above_median=(((double)rank/size_world)>0.5)?1:0;
        int new_rank;  
        if(above_median)
        {   
            MPI_Comm_split(comm, 0, 0, &above_comm);
            MPI_Comm_rank(comm, &new_rank);
            build_tree(tree, lnode_id,local_start_id, pts, projections, center_idx, n_dims,above_comm,new_rank,start_npoints,start_ndims); //center_idx happens to be the number of points in the set
        }
        else
        {   
            MPI_Comm_split(comm, 1, 0, &below_comm);
            MPI_Comm_rank(comm, &new_rank);
            build_tree(tree, rnode_id, local_start_id,pts + center_idx, projections + center_idx, n_points - center_idx, n_dims,below_comm,new_rank,start_npoints,start_ndims);
        }
    }
    else
    {   
        build_tree(tree, lnode_id, local_start_id,pts, projections, center_idx, n_dims,comm,rank,start_npoints,start_ndims); //center_idx happens to be the number of points in the set
        build_tree(tree, rnode_id, local_start_id,pts + center_idx, projections + center_idx, n_points - center_idx, n_dims,comm,rank,start_npoints,start_ndims);
    }
        
    return;
     
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

void init_tree(node* tree,long n_nodes)
{
    for(long i=0;i<n_nodes;i++){
        tree[i].id=-1;
    }
}

int main(int argc, char *argv[])
{   
    MPI_Init(&argc, &argv); /*START MPI */
    int rank, size, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); /*DETERMINE RANK OF THIS PROCESSOR*/
    MPI_Comm_size(MPI_COMM_WORLD, &size); /*DETERMINE TOTAL NUMBER OF PROCESSORS*/
    MPI_Get_processor_name(processor_name, &name_len);

    double exec_time;
    double **pts;
    node* tree;

    //____________START_TIME_BENCHMARK_____________ 
    exec_time = -omp_get_wtime();
    //generates dataset
    pts = get_points(argc, argv);
    double* pts_first_position = pts[0];


    int n_dims = atoi(argv[1]); //number of dimensions
    long n_points = atoi(argv[2]); //number of points in the set

    int blocklengths[5] = {1,n_dims,1,1,1};
    MPI_Datatype types[5] = {MPI_DOUBLE, MPI_DOUBLE, MPI_LONG,MPI_LONG,MPI_LONG};
    MPI_Aint     offsets[5];
    offsets[0] = offsetof(node, radius);
    offsets[1] = offsetof(node, center);
    offsets[2] = offsetof(node, L);;
    offsets[3] = offsetof(node, R);
    offsets[4] = offsetof(node, id);
    MPI_Type_create_struct(5, blocklengths, offsets, types, &mpi_node_type);
    MPI_Type_commit(&mpi_node_type);


    double* projections = (double*)malloc(n_points*sizeof(double)); //array to store pseudo-projections the locate the point in the line 
    long n_nodes = 2 * n_points - 1; //number of nodes in the tree
    tree = (node*)malloc(n_nodes*sizeof(node));
    init_tree(tree,n_nodes);
    
    // construct THE TREE
    long n_local_nodes = 2*(n_points/size+size);
    long local_start_id= rank*n_local_nodes;

    if(rank==0)
        printf("%d %ld\n", n_dims, n_nodes);
    
    build_tree(tree,0,local_start_id, pts, projections, n_points, n_dims,MPI_COMM_WORLD,rank,n_points,n_dims);
    
    MPI_Barrier(MPI_COMM_WORLD);    

    MPI_Finalize();
    //____________END_TIME_BENCHMARK_____________

    exec_time += omp_get_wtime();
    if(rank==0)
        fprintf(stderr, "%.1lf\n", exec_time);
    
    free(projections);
    free(pts_first_position);
    free(pts);
}
