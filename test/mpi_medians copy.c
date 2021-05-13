#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <mpi.h>
#include "gen_points.h"
#include "unsorted_median.h"


int main(int argc, char *argv[])
{
    double exec_time;
    double **pts;
    long n_local_points;

    MPI_Init(&argc, &argv);

    pts = get_points_mpi(argc, argv, MPI_COMM_WORLD,&n_local_points);
    
    double* pts_first_position = pts[0]; //position of the first element of pts (pts is sorted so it would be lost and hard to free)
    int n_dims = atoi(argv[1]); //number of dimensions 
    long n_points = atoi(argv[2]); //number of points in the set
    double* projections = (double*)malloc(n_local_points*sizeof(double)); //array to store pseudo-projections the locate the point in the line 
    long n_nodes = 2 * n_points - 1; //number of nodes in the tree
    long n_local_nodes = 2*n_points-1;
    //tree = (node*)malloc(n_nodes*sizeof(node));

    //Construct THE TREE
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
	
	//median_right = getKsmallest(projections, n_points/2 , n_points);

    MPI_Finalize();

    free(projections);
    free(pts_first_position);
    free(pts);

}