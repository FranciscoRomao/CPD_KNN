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
    double* pts;
    long n_points;

    MPI_Init(&argc, &argv);

    pts = get_points_1dim_mpi(argc, argv, MPI_COMM_WORLD, &n_points);
    
    double median;

    //pseudo-projection of all points (enough to know relative positions)

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
	
	median = getKsmallest(pts, n_points/2 , n_points);

    MPI_Finalize();

    free(pts);

}