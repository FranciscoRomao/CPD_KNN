#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "gen_points.h"
#include "geometry.h"
#include "unsorted_median.h"
#include <string.h>

void flag (int x, char* name)
{   
    printf("%d %s\n", x, name);
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

void print_pjs(double* pjs, long n_pts){
    for(long i=0;i<n_pts;i++){
        printf("%lf ",pjs[i]);
    }
}

void print_pts(double** pts, long n_pts, int n_dims){
    for(long i=0;i<n_pts;i++){
        for(int j=0;j<n_dims;j++)
            printf("%lf ",pts[i][j]);
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


void build_tree(long node_id, double **pts, double* projections, long n_points, int n_dims, MPI_Comm comm, int rank, long start_npoints)
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
    if (rank!=0 && (n_points == 1 || n_points==2)) //if the node is a leaf
    {
        return;
    }
    node foo;
    if(rank==0){
        //flag(76,my_name);
        foo.id=node_id;

        if (n_points == 1) //if the node is a leaf
        {
            foo.R = -1;
            foo.L = -1;
            foo.radius = 0;
            //flag(77,my_name);
            foo.center = pts[0];
            //flag(78,my_name);
            print_Node(foo,n_dims);
            //flag(79,my_name);
            return;
        }

        lnode_id = node_id + 1; //indice of the left child

        foo.L=lnode_id;

        foo.center = (double *)malloc(n_dims * sizeof(double));

        //flag(94,my_name);

        if(n_points == 2) //only two points in the set -> easier/lesser operations
        {   
            center_idx = 1;
            rnode_id = node_id + 2 * center_idx;
            foo.R = rnode_id;
            /*printf("After comm:\n");
            printf("r:%ld l:%ld c:%ld n:%ld\n",rnode_id,lnode_id,center_idx,n_points);
            print_pts(pts,start_npoints,n_dims);printf("\n");
            print_pjs(projections,start_npoints);printf("\n");*/
            //flag(100,my_name);
            double* aux;
            //printf("node_id:%ld pts[0][0]:%lf pts[1][0]:%lf\n\n",node_id,pts[0][0],pts[1][0]);
            if (pts[0][0] > pts[1][0]) //sort ascending if not already
            {   
                aux = pts[0];
                pts[0] = pts[1];
                pts[1] = aux;
            }
            //flag(107,my_name);
            //median of two points is its average
            for(int i=0; i<n_dims; i++)
            {
                foo.center[i] = (pts[0][i] + pts[1][i]) / 2;
            }
            //flag(113,my_name);
            //both points are equidistant to the center
            foo.radius = distance(n_dims, pts[0], foo.center);
            //printf("radius:%lf\n",foo.radius);
            //flag(115,my_name);
            //build leafs
            print_Node(foo,n_dims);
            build_tree(lnode_id, pts, NULL, 1, n_dims,comm,rank,start_npoints); //center_idx happens to be the number of points in the set
            build_tree(rnode_id,pts+1, NULL, 1, n_dims,comm,rank,start_npoints);
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
            //flag(136,my_name);
    
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
                //flag(156,my_name);
                for(int i=0; i<n_dims; i++)
                {
                    foo.center[i] = (center_l[i] + center_r[i]) / 2;
                }
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
                //flag(174,my_name);
                //compute and set center of the node
                orthogonal_projection(n_dims, pts[median_idx], pts[idx_fp[0]], pts[idx_fp[1]], foo.center);

                compare_with_median(projections, pts, median, n_points,n_dims);

                center_idx = (n_points - 1) / 2;
            }

            //compute the furthest point of the ball center and thus the radius
            radius_candidate[0] = distance(n_dims, pts[idx_fp[0]], foo.center);
            radius_candidate[1] = distance(n_dims, pts[idx_fp[1]], foo.center);

            if (radius_candidate[1]>radius_candidate[0]) 
            {
                foo.radius = radius_candidate[1];
            }
            else 
            {
                foo.radius = radius_candidate[0];
            }
            //flag(195,my_name);
            rnode_id = node_id + 2 * center_idx;
            foo.R = rnode_id;
            print_Node(foo,n_dims);
            //flag(201,my_name);
        }
    }

    int size_world;
    MPI_Comm_size(comm, &size_world);
    if(size_world>=2 && (n_points != 1 && n_points!=2))
    {   
        MPI_Bcast(pts[0],start_npoints*n_dims,MPI_DOUBLE,0,comm);
        MPI_Bcast(projections,start_npoints,MPI_DOUBLE,0,comm);
        MPI_Bcast(&rnode_id,1,MPI_LONG,0,comm);
        MPI_Bcast(&lnode_id,1,MPI_LONG,0,comm);
        MPI_Bcast(&center_idx,1,MPI_LONG,0,comm);
        MPI_Bcast(&n_points,1,MPI_LONG,0,comm);
        /*if(rank==0){
        printf("Before comm:\n");
        printf("r:%ld l:%ld c:%ld n:%ld\n",rnode_id,lnode_id,center_idx,n_points);
        print_pts(pts,n_points,n_dims);printf("\n");
        print_pjs(projections,n_points);printf("\n");}
        else{
            printf("After comm:\n");
            printf("r:%ld l:%ld c:%ld n:%ld\n",rnode_id,lnode_id,center_idx,n_points);
            print_pts(pts,n_points,n_dims);printf("\n");
            print_pjs(projections,n_points);printf("\n");
        }*/

        MPI_Comm above_comm;
        MPI_Comm below_comm;
        int above_median=(((double)rank/size_world)>=0.5)?1:0;
        //printf("above median:%d rank:%d\n",above_median,rank);
        int new_rank;
        if(above_median)
        {   
            MPI_Comm_split(comm, 0, 0, &below_comm);
            MPI_Comm_rank(below_comm, &new_rank);
            //printf("NEW_RANK: %d above\n", new_rank);
            build_tree(lnode_id, pts, projections, center_idx, n_dims,below_comm,new_rank,start_npoints); //center_idx happens to be the number of points in the set
        }
        else
        {  
            //printf("splitting below!\n");
            MPI_Comm_split(comm, 1, 0, &above_comm);
            MPI_Comm_rank(above_comm, &new_rank);
            //printf("NEW_RANK: %d below\n", new_rank);
            build_tree(rnode_id,pts + center_idx, projections + center_idx, n_points - center_idx, n_dims,above_comm,new_rank,start_npoints);
        }
    }
    else
    {   
        /*printf("Before comm:\n");
        printf("r:%ld l:%ld c:%ld n:%ld\n",rnode_id,lnode_id,center_idx,n_points);
        print_pts(pts,n_points,n_dims);printf("\n");
        print_pjs(projections,n_points);printf("\n");*/
        build_tree(lnode_id,pts, projections, center_idx, n_dims,comm,rank,start_npoints); //center_idx happens to be the number of points in the set
        build_tree(rnode_id,pts + center_idx, projections + center_idx, n_points - center_idx, n_dims,comm,rank,start_npoints);
    }
        
    return;
     
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
    int n_dims = atoi(argv[1]); //number of dimensions
    long n_points = atoi(argv[2]); //number of points in the set
    double* _p_arr;
    //____________START_TIME_BENCHMARK_____________ 
    exec_time = -omp_get_wtime();

    //generates dataset
    if(rank==0)
        pts = get_points(argc, argv);
    else{
        _p_arr = (double *)malloc(n_dims * n_points * sizeof(double));
        pts = (double **)malloc(n_points * sizeof(double *));
        for (long i = 0; i < n_points; i++)
        pts[i] = &_p_arr[i * n_dims];
    }
    double* pts_first_position = pts[0];


    double* projections = (double*)malloc(n_points*sizeof(double)); //array to store pseudo-projections the locate the point in the line 
    long n_nodes = 2 * n_points - 1; //number of nodes in the tree

    if(rank==0)
        printf("%d %ld\n", n_dims, n_nodes);
    
    MPI_Comm comm = MPI_COMM_WORLD;
    build_tree(0, pts, projections, n_points, n_dims, comm ,rank,n_points);

    //____________END_TIME_BENCHMARK_____________
    MPI_Barrier(MPI_COMM_WORLD);
    exec_time += omp_get_wtime();
    if(rank==0)
        fprintf(stderr, "%.1lf\n", exec_time);
    
    free(projections);
    free(pts_first_position);
    free(pts);
    MPI_Finalize();
    return 0;
}
