#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <mpi.h>
#include "gen_points.h"
#include "unsorted_median.h"

#define SAMPLE 3

#define DEBUG

int main(int argc, char *argv[])
{
    double* pts;
    double* pts_test;
    long n_points = atol(argv[2]);;
    long n_points_test = atol(argv[2]);
    int rank;
    int n_procs;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

    if(n_procs * SAMPLE > n_points)
    {
        if(!rank)
        {
            printf("Isto tem de ser implementado para correr numa máquina só, porque senão dá bugada\n"); fflush(stdout);
        }
        return -1;
    }
    
    pts = get_points_1dim_mpi(argc, argv, MPI_COMM_WORLD, &n_points);


    #ifdef DEBUG
        if(!rank)
        {
            pts_test = get_points_1dim(argc, argv);
            //printf("pontos gerados\n");
            //printArray(pts_test, n_points_test);
            printf("#-------------------------Mediana real: %lf\n", medianSort(pts_test, n_points_test));
            fflush(stdout);
        }
    #endif

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
	
	//median = getKsmallest(pts, n_points/2 , n_points);
    MPI_Barrier(MPI_COMM_WORLD);
    PSRS(pts, n_points);

    //printf("FIMMMMMM, rank %d\n", rank);
    //fflush(stdout);
    MPI_Finalize();

    //free(pts);

}