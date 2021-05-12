#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define RANGE 10
#define DEBUG 1

extern void print_point(double *, int);

double **create_array_pts(int n_dims, long np)
{
    double *_p_arr;
    double **p_arr;

    _p_arr = (double *)malloc(n_dims * np * sizeof(double));
    p_arr = (double **)malloc(np * sizeof(double *));
    if ((_p_arr == NULL) || (p_arr == NULL))
    {
        printf("Error allocating array of points, exiting.\n");
        exit(4);
    }

    for (long i = 0; i < np; i++)
        p_arr[i] = &_p_arr[i * n_dims];

    return p_arr;
}

double **get_points(int argc, char *argv[])
{
    double **pt_arr;
    unsigned seed;
    int n_dims;
    long np;
    long i;
    int j;

    if (argc != 4)
    {
        printf("Usage: %s <n_dims> <n_points> <seed>\n", argv[0]);
        exit(1);
    }

    n_dims = atoi(argv[1]);
    if (n_dims < 2)
    {
        printf("Illegal number of dimensions (%d), must be above 1.\n", n_dims);
        exit(2);
    }

    np = atol(argv[2]);
    if (np < 1)
    {
        printf("Illegal number of points (%ld), must be above 0.\n", np);
        exit(3);
    }

    seed = atoi(argv[3]);
    srandom(seed);

    pt_arr = (double **)create_array_pts(n_dims, np);

    for (i = 0; i < np; i++)
        for (j = 0; j < n_dims; j++)
            pt_arr[i][j] = RANGE * ((double)random()) / RAND_MAX;


    return pt_arr;
}

double **get_points_mpi(int argc, char *argv[], MPI_Comm comm, long* n_local_points)
{
    double **pt_arr;
    unsigned seed;
    int rank;
    int n_procs;
    int n_dims;
    long np;
    long i;
    int j;
    int k;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &n_procs);
    
    n_dims = atoi(argv[1]);
    if (n_dims < 2)
    {
        printf("Illegal number of dimensions (%d), must be above 1.\n", n_dims);
        exit(2);
    }

    np = atol(argv[2]);
    if (np < 1)
    {
        printf("Illegal number of points (%ld), must be above 0.\n", np);
        exit(3);
    }

    int full_split = (int) np/n_procs;
    int last_split = (int) full_split + np%n_procs;
    pt_arr = (double **)create_array_pts(n_dims, (long) last_split);

    if(rank==0)
    {
        if (argc != 4)
        {
            printf("Usage: %s <n_dims> <n_points> <seed>\n", argv[0]);
            exit(1);
        }

        seed = atoi(argv[3]);
        srandom(seed);
        
        for(k = 1; k<n_procs; k++)
        {
            for (i = 0; i<full_split; i++)
            {
                for (j = 0; j<n_dims; j++)
                {
                    pt_arr[i][j] = RANGE * ((double)random()) / RAND_MAX;
                }
            }
            MPI_Send(&(pt_arr[0][0]), full_split*n_dims, MPI_DOUBLE, k, 0, comm);
        }

        for (i = 0; i<last_split; i++)
        {
            for (j = 0; j<n_dims; j++)
            {
                pt_arr[i][j] = RANGE * ((double)random()) / RAND_MAX;
            }
        }
	np = last_split;
    }
    if(rank != 0)
    {   
        MPI_Recv(&(pt_arr[0][0]), full_split*n_dims, MPI_DOUBLE, 0, 0, comm, MPI_STATUS_IGNORE);
    	np = full_split;
    }

    *n_local_points=np;

    #ifdef DEBUG
    	printf("Rank: %d\n", rank);
    	fflush(stdout);
    	for (i = 0; i < np; i++)
        {
		printf("pt[%ld]: ", i);
 		fflush(stdout);
		for(j = 0; j < n_dims; j++)
		{
			printf("%lf ", pt_arr[i][j]);
			fflush(stdout);
		}
		printf("\n");
		fflush(stdout);
	}   
    #endif

    return pt_arr;
}
