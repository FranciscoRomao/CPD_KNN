#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include <mpi.h>
#include <string.h>
#include "unsorted_median.h"
#include "gen_points.h"

#define IGNORE_NUM -1010101010
#define SAMPLE 3

int rank;

/**
 * Compare function used inside quicksort
 * @param   a   First number to compare
 * @param   b   Second number to compare
 * @return  Returns 1 if a>b, -1 if a<b and 0 if a=b  
 */
int cmpfunc (const void * a, const void * b)
{
  if (*(double*)a > *(double*)b) return 1;
  else if (*(double*)a < *(double*)b) return -1;
  else return 0;
}

/**
 * Copies n_items from the source array to the destination array
 * @param   dest    Destination array
 * @param   sourc   Source array
 * @param   n_items Number of items to copy  
 */
void arraycpy(double* dest, double* sourc, int n_items)
{
    for(int i=0; i<n_items; i++)
    {
        dest[i] = sourc[i];
    }
}


/**
 * Compute the median of a vector using quicksort
 * @param   vector    Vector to find the median of
 * @param   n_items   Number of items of the vector
 * @return  Median of the vector  
 */
double sorted_median(double *vector, int n_items)
{
    double result;
    
    if(n_items==1)
    {
        result = (double)(vector[0]);
        return result;
    }

    if(n_items==2)
    {
        result = (double)(vector[0] + vector[1])/2;
        return result;
    }

    qsort(vector, n_items, sizeof(double), cmpfunc);

    if(n_items % 2 != 0)
    {
        return vector[n_items/2];
    }
    else
    {
        //printf("Sorted: %lf | %lf\n", vector[n_items/2], vector[n_items/2 - 1]);
        return 0.5 * ((double)vector[n_items/2] + vector[n_items/2 - 1]);
    }
}


/**
 * Builds two sets from the given vector, the left set has all the points lower than the median
 * and the right set has all the points higher than the median.
 * [Warning: The median is a prediction so the sets can be larger than half of the vector] 
 * @param   setL    Left set 
 * @param   setR    Rigth set
 * @param   counterL    Return by reference the number of elements on the left set
 * @param   counterR    Return by reference the number of elements on the rigth set
 * @param   vector  Vector to separate in two sets
 * @param   n_items     Number of elements on the vector
 * @param   median  Median to use to separate the vector
 * @return  Number of elements on the left set which is equal to the indice of the media
 */
double buildSet(double *setL, double *setR, int* counterL, int* counterR, double *vector, int n_items, double median)
{
    *counterL=0;
    *counterR=0;

    for(int i=0; i<n_items; i++)
    {
        //printArray(setL, *counterL);
        if(vector[i]<median)
        {
            (*counterL)++;
            setL[(*counterL)-1] = vector[i];
        }
        else
        {
            (*counterR)++;
            setR[(*counterR)-1] = vector[i];
        }
    }

    //printArray(setL, *counterL);
    //printf("Sets criados: setL %d items, setR %d items, idx da mediana: %d\n\n", *counterL, *counterR, (*counterL));
    return (*counterL);
}


/**
 * Compute the median of a vector using quicksort
 * @param   vector    Vector to find the median of
 * @param   n_items   Number of items of the vector
 * @return  Median of the vector
 */
double median(double *vector, int n_items)
{
    int full_splits = n_items/5;
    int semi_splits = n_items % 5;
    int i=0;
    double result;
    int n_medians = 0;
    double *medians = (double *)malloc(n_items/2 * sizeof(double));
    
    //printArray(vector, 6);

    for(i=0; i<full_splits; i++)
    {
        medians[i] = sorted_median(vector + 5*i, 5);
        // printf("Median do grupo %d: %lf\n", i, medians[i]);
    }

    if(semi_splits != 0)
    {
        medians[i] = sorted_median(vector + 5*full_splits, semi_splits);
        // printf("Median do semi grupo %d: %lf\n", i, medians[i]);
    }

    n_medians = full_splits + (semi_splits!=0 ? 1 : 0);

    if(n_medians <= 5)
    {
        result = sorted_median(medians, n_medians);
        //printf("Mediana das medianas: %lf\n", result);
        //result = buildSet(setL, setR, vector, n_items, result);
        free(medians);
        //printf("Mediana das medianas: %lf\n", result);
        return result;
    }

    result = median(medians, full_splits + (semi_splits!=0 ? 1 : 0));
    free(medians);

    //printf("Mediana das medianas: %lf\n", result);

    //result = buildSet(setL, setR, vector, n_items, result);

    return result;
}


/**
 * Finds the smallest element of a vector
 * @param   vector    Vector to find the median of
 * @param   n_items   Number of items of the vector
 * @return  Values of the smallest element
 */
double getSmallest(double *vector, int n_items)
{
    double min = DBL_MAX;

    for(int i=0; i<n_items; i++)
    {
        if(vector[i] < min)
            min = vector[i];
    }
    return min;
}

int vectorSum(int* vector, int n)
{
    int count=0;

    for(int i=0; i<n; i++)
    {
        count += vector[i];
    }
    return count;
}

double PSRS(double* vector, long n_items)
{
    int n_procs;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    
    double *vector_cpy = (double *)malloc(n_items * sizeof(double));
    double* sample = (double *)malloc(SAMPLE * sizeof(double));
    double* recvValues;
    int* aux_vector = (int*)malloc(n_procs*2*sizeof(int));;
    int result_idx;
    int n_total;
    double** partitions = (double**)create_array_pts(n_items, n_procs);
    double** partitions2recv = (double**)create_array_pts(n_items, n_procs);
    double_memset(partitions2recv[0], IGNORE_NUM, n_items*n_procs);

    double* finalVector = (double*)malloc(sizeof(double)*2*n_items); //_--------------------------Não sei bem que dimensão colocar aqui
    int elems_recvd=0;
    
    MPI_Request* requests = (MPI_Request*)malloc(sizeof(MPI_Request)*n_procs);
    double* partitions_limits = (double*)malloc(sizeof(double)*n_procs);
    int n_elem;
    int n_my_elem;
    int n_to_recv;
    int i;
    double aux_limit;
    int aux;
    double aux1_db;
    double aux2_db;

    //printf("aqui1, rank: %d\n", rank);
    //printArray(vector, n_items);
    
    if(!rank)
    {
        recvValues = (double*)malloc(n_procs*SAMPLE*sizeof(double));
    }

    //printf("aqui2, rank: %d\n", rank);
    //fflush(stdout);

    for(i=0; i<SAMPLE; i++) //Podemos experimentar retirar a sample apartir de i=1
    {
        //printf("aqui3, rank: %d\n", rank);
        //fflush(stdout);
        sample[i] = getKsmallest(vector, (n_items/SAMPLE)*i, n_items);
        //printf("sample %d = %lf, rank %d\n", i, sample[i], rank);
        //fflush(stdout);
    }

    MPI_Gather(&(sample[0]), SAMPLE, MPI_DOUBLE, &(recvValues[0]), SAMPLE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    //printf("aqui5, rank: %d\n", rank);
    //fflush(stdout);
    
    /*if(!rank)
    {
        printf("Samples recebidas: \n");
        printArray(recvValues, SAMPLE*n_procs);
    }*/

    if(!rank)
    {
        for(i=1; i<n_procs; i++)
        {
            partitions_limits[i-1] = getKsmallest(recvValues, SAMPLE*i-1, n_procs*SAMPLE); //cada processador enviou uma SAMPLE
        }
        //printf("Limites do rank 0:\n");
        //printArray(partitions_limits, n_procs-1);
    }

    MPI_Bcast(partitions_limits, n_procs-1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    //printf("Limites do rank %d\n", rank);
    //printArray(partitions_limits, n_procs-1);

    for(i=0; i<n_procs; i++)
    {
        if(i==0)
        {
            if(i==rank)
            {
                aux_limit = getKsmallest(vector, 0, n_items)-1;

                //printf("aux_limit -0: %lf\n", aux_limit);
                //printf("vector inicial: ");
                //printArray(vector, n_items);
                n_elem = getIntervalSet(partitions[i], vector, n_items, aux_limit, partitions_limits[i]);//O zero é arriscado depende dos valores que estamos a guardar
                //printf("Particao %d do rank %d\n-----------------\n", i, rank);
                //printArray(partitions[i], n_elem);
                //printf("-----------------\n", i, rank);
                //fflush(stdout);
                elems_recvd = appendArray(finalVector, 0, partitions[i], n_elem, IGNORE_NUM);
                continue;
            }
            n_elem = getIntervalSet(partitions[i], vector, n_items, partitions_limits[i-1], partitions_limits[i]);
            //printf("Particao %d do rank %d\n-----------------\n", i, rank);
            //printArray(partitions[i], n_elem);
            //printf("-----------------\n", i, rank);
            //fflush(stdout);
            MPI_Isend(&n_elem, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &(requests[i]));
            MPI_Isend(partitions[i], n_elem, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &(requests[i]));
        }
        else
        {
            if(i==n_procs-1)
            {
                if(i==rank)
                {
                    //printArray(vector, n_items);
                    aux_limit = getKsmallest(vector, n_items-1, n_items);
                    //printf("rank %d, Limit superior: %lf\n, inferior: %lf", rank, aux_limit, partitions_limits[i-1]);
                    //printArray(vector, n_items);
                    n_elem = getIntervalSet(partitions[i], vector, n_items, partitions_limits[i-1], aux_limit);//O zero é arriscado depende dos valores que estamos a guardar
                    //printf("---Particao %d do rank %d\n-----------------\n", i, rank);
                    //printArray(partitions[i], n_elem);
                    //printf("-----------------\n", i, rank);
                    //fflush(stdout);
                    elems_recvd = appendArray(finalVector, 0, partitions[i], n_elem, IGNORE_NUM);
                    continue;
                }
                aux_limit = getKsmallest(vector, n_items-1, n_items);
                //printf("rank %d, Limit superior: %lf\n, inferior: %lf", rank, aux_limit, partitions_limits[i-1]);
                //printArray(vector, n_items);
                
                n_elem = getIntervalSet(partitions[i], vector, n_items, partitions_limits[i-1], aux_limit);
                //printf("---Particao %d do rank %d\n-----------------\n", i, rank);
                //printArray(partitions[i], n_elem);
                //printf("-----------------\n", i, rank);
                //fflush(stdout);
                MPI_Isend(&n_elem, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &(requests[i]));
                MPI_Isend(partitions[i], n_elem, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &(requests[i]));
            }
            else
            {
                if(i==rank)
                {
                    n_elem = getIntervalSet(partitions[i], vector, n_items, partitions_limits[i-1], partitions_limits[i]);
                    //printf("Particao %d do rank %d\n-----------------\n", i, rank);
                    //printArray(partitions[i], n_elem);
                    //printf("-----------------\n", i, rank);
                    //fflush(stdout);
                    elems_recvd = appendArray(finalVector, 0, partitions[i], n_elem, IGNORE_NUM);
                    continue;
                }
                n_elem = getIntervalSet(partitions[i], vector, n_items, partitions_limits[i-1], partitions_limits[i]);
                //printf("Particao %d do rank %d\n-----------------\n", i, rank);
                //printArray(partitions[i], n_elem);
                //printf("-----------------\n", i, rank);
                //fflush(stdout);
                MPI_Isend(&n_elem, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &(requests[i]));
                MPI_Isend(partitions[i], n_elem, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &(requests[i]));
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    for(i=0; i<n_procs; i++)
    {
        if(i==rank)
            continue;

        MPI_Recv(&n_to_recv, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //printf("Rank %d recebeu %d do rank %d\n", rank, n_to_recv, i);
        MPI_Recv(partitions2recv[i], n_to_recv , MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //printf("--------Rank %d recebeu do rank %d--------\n", rank, i);
        //printArray(partitions2recv[i], n_to_recv);
        elems_recvd = appendArray(finalVector, elems_recvd, partitions2recv[i], n_to_recv, IGNORE_NUM);
    }

    //printf("Valores finais do rank %d\n", rank);
    //printArray(finalVector, elems_recvd);

    aux = elems_recvd;

    n_items = elems_recvd;

    //o nome partitions_limits não tem lógica aqui, mas o objectivo é reutilizar o vetor que já não seja necessário

    MPI_Gather(&aux, 1, MPI_INT, &(aux_vector[0]), 1, MPI_INT, 0, MPI_COMM_WORLD);
    

    /*if(!rank)
    {
        printf("Isto: ");
        printf("%d, %d, %d\n", aux_vector[0], aux_vector[1], aux_vector[2]);
        fflush(stdout);
    }*/

    int median_holder1;
    int median_holder2;
    int median_idx1;
    int median_idx2;

    if(!rank)
    {
        n_total = getVectorSum(aux_vector, n_procs);
        
        //printf("N total: %d\n", n_total);
        //fflush(stdout);

        if(n_total%2 != 0)
        {
            median_idx1 = find_K_idx(aux_vector, n_total/2, n_total, &median_holder1);

            //printf("median1: %d\n", median_idx1);
            //printf("holder1: %d\n", median_holder1);
            //fflush(stdout);

            for(i = 0; i<n_procs*2; i++)
                aux_vector[i] = -1;
            
            aux=0;

            for(i=0; i<n_procs; i++)
            {
                if(i == median_holder1)
                    aux_vector[aux] = median_idx1;

                MPI_Isend(&(aux_vector[aux]), 2, MPI_INT, i, 2, MPI_COMM_WORLD, &(requests[i]));
                aux++;
            }
        }
        else
        {
            median_idx1 = find_K_idx(aux_vector, n_total/2, n_total, &median_holder1);
            median_idx2 = find_K_idx(aux_vector, (n_total/2)-1, n_total, &median_holder2);

            //printf("median1: %d\n", median_idx1);
            //printf("median2: %d\n", median_idx2);
            //printf("holder1: %d\n", median_holder1);
            //printf("holder2: %d\n", median_holder1);
            //fflush(stdout);

            for(i = 0; i<n_procs*2; i++)
                aux_vector[i] = -1;

            aux=0;
            
            for(i = 0; i<n_procs; i++)
            {
                if(i == median_holder1)
                {
                    aux_vector[aux] = median_idx1;
                }
                else
                    aux_vector[aux] = -1;
                aux++;

                if(i == median_holder2)
                {
                    aux_vector[aux] = median_idx2;
                }
                else
                    aux_vector[aux] = -1;
                aux++;
            }

            //printf("%d %d %d %d %d %d\n", aux_vector[0], aux_vector[1], aux_vector[2], aux_vector[3], aux_vector[4], aux_vector[5]);
            //fflush(stdout);

            for(i=0; i<n_procs; i++)
                MPI_Isend(&(aux_vector[2*i]), 2, MPI_INT, i, 2, MPI_COMM_WORLD, &(requests[i]));
        }
    }
    
    MPI_Recv(&(aux_vector[0]), 2, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    //printf("rank %d recebeu: %d %d\n", rank, aux_vector[0], aux_vector[1]);
    //fflush(stdout);
    
    for(i=0; i<2; i++)
    {
        if(aux_vector[i] != -1)//aux_limit é reutilizado
        {
            aux_limit = getKsmallest(finalVector, aux_vector[i], n_items);
            //printf("Buscado valor: %lf", aux_limit);
            MPI_Isend(&aux_limit, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &(requests[i]));
        }
    }

    if(!rank)
    {
        if(n_total%2 == 0)
        {
            MPI_Recv(&aux1_db, 1, MPI_DOUBLE, median_holder1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&aux2_db, 1, MPI_DOUBLE, median_holder1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("---------------------------------Mediana final: %lf\n", (aux1_db + aux2_db)/2);
            fflush(stdout);
            return (aux1_db + aux2_db)/2;
        }
        else
        {
            MPI_Recv(&aux1_db, 1, MPI_DOUBLE, median_holder1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("---------------------------------Mediana final: %lf\n", aux1_db);
            fflush(stdout);
            return aux1_db;
        }
    }

    //printf("chegou ao fim||||||||||||||\n");
    //fflush(stdout);

    //Todos menos o rank 0 vão fazer return de 0. O rank 0 faz return da mediana
    return 0;
}


/**
 * Finds the index of the closer and lower element to result on the vector
 * @param   vector  Vector to find the element
 * @param   n_items Number of items of the array
 * @return  Index of the lower neightbor
 */
long getLowerNeighborIdx(double* vector, long n_items, double result)
{
    long neighbor_idx;
    double min_diff = -1;
    double diff;
    int flag = 0;

    for (int i = 0; i < n_items; i++)
    {
        if (vector[i] > result)
        {
            continue;
        } 
        else if ((vector[i] < result))
        {
            diff = result-vector[i];
            if ( diff < min_diff || min_diff == -1 )
            {
                min_diff = diff;
                neighbor_idx = i;
            }
        }
        else if (vector[i] == result)
        {
            if (flag == 0)
            {
                flag = 1;
            }
            else if (flag == 1)
            {
                neighbor_idx = i;
                break;
            }
        }
    }

    return neighbor_idx;
}


/**
 * Reorders the vectors projections and *pts. The left set has all the points with values of the projections
 * lower than the median and the right set has all the points higher than the median. The *pts array is swapped
 * acordingly to the projections array
 * @param   projections Reference array to reorder the arrays
 * @param   pts         Array to be swapped acordingly
 * @param   median      Median to use to separate the projections vector
 * @param   n_points    Number of elements on the array
 */
void compare_with_median(double* projections, double** pts, double median, long n_points)
{
    long left_pivot = 0;
    long right_pivot = n_points -1;
    long points_remaining;
    double aux;
    double* aux_pt;

    while (left_pivot < right_pivot) {
       if (projections[left_pivot] >= median) {
           if (projections[right_pivot] < median)
           {
                aux = projections[right_pivot];
                projections[right_pivot] = projections[left_pivot];
                projections[left_pivot] = aux;

                aux_pt = pts[right_pivot];
                pts[right_pivot] = pts[left_pivot];
                pts[left_pivot] = aux_pt;
            } else {
               right_pivot--;
           }
       } else {
           left_pivot++;
       }
    }

    points_remaining = n_points/2 - left_pivot;
    if (points_remaining == 0)
        return;

    for(int i=left_pivot; i<n_points; i++)
    {
        if(projections[i] == median)
        {
            aux = projections[left_pivot];
            projections[left_pivot] = projections[i];
            projections[i] = aux;

            aux_pt = pts[right_pivot];
            pts[right_pivot] = pts[left_pivot];
            pts[left_pivot] = aux_pt;
            
            left_pivot++;
            points_remaining--;
            if (points_remaining == 0)
                break;
        }
    }
}

/**
 * Searchs one value in an array and returns its idx
 * @param   projections    Array to search in
 * @param   n_points   number of points in the array 
 * @param   value   value to find
 * @return  Returns the index of value from the projection vector
 */
long find_idx_from_value(double *projections, long n_points, double value)
{
    for(int i=0; i<n_points; i++)
    {
        if(projections[i] == value)
            return i;
    }
    return -1;
}

double getKsmallest(double* vector, long k, long n_items)
{
    int L_items;
    int R_items;
    double result;
    int result_idx;
    double *setL = (double *)malloc(n_items * sizeof(double));
    double *setR = (double *)malloc(n_items * sizeof(double));
    double *vector_cpy = (double *)malloc(n_items * sizeof(double));
    
    arraycpy(vector_cpy, vector, n_items);
    
    //printArray(vector, n_items);

    result = median(vector_cpy, n_items);

    arraycpy(vector_cpy, vector, n_items);

    result_idx = buildSet(setL, setR, &L_items, &R_items, vector_cpy, n_items, result);

    //printf("Target idx: %d, result idx: %d\n", k, result_idx);

    while(k > result_idx || k < result_idx)
    {
        if(k < result_idx)
        {
            arraycpy(vector_cpy, setL, L_items);

            result = median(vector_cpy, L_items);
            result_idx = buildSet(setL, setR, &L_items, &R_items, vector_cpy, L_items, result);
            //printf("Target idx: %d, result idx: %d\n", k, result_idx);
        }

        if(k > result_idx)
        {
            arraycpy(vector_cpy, setR, R_items);
            {
                k = k - result_idx;
                result = median(vector_cpy, R_items);
                result_idx = buildSet(setL, setR, &L_items, &R_items, vector_cpy, R_items, result);
                //printf("Target idx: %d, result idx: %d\n", k, result_idx);
            }
        }
    }

    result = getSmallest(setR, R_items);
    free(vector_cpy);
    free(setL);
    free(setR);
    
    return result;
}

double* linerize(double** matrix, int nx, int ny, double* vector)//x é o numero de vetores na matrix e y é o numero de elementos dentro do vetor
{
    for(int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            vector[i*nx+j] = matrix[i][j];
        }
    }
}

double getIntervalSet(double *interval, double *vector, int n_items, double comparatorL, double comparatorR)
{
    int counter=0;

    for(int i=0; i<n_items; i++)
    {
        //printArray(setL, *counterL);
        if(vector[i]>comparatorL && vector[i]<=comparatorR)
        {
            counter++;
            interval[counter-1] = vector[i];
        }
    }
    //printArray(setL, *counterL);
    //printf("Sets criados: setL %d items, setR %d items, idx da mediana: %d\n\n", *counterL, *counterR, (*counterL));
    return counter;
}

int appendArray(double* dest, int begin_idx, double* source, int n_elem, double ignore)
{
    int j=0;

    for(int i=0; i<n_elem; i++)
    {
        /*if(source[i] == (double)ignore)
        {
            //printf("rank aqui ----- %d", rank);
            //if(rank == 1)
            //    printf("entrou\n");
            continue;
        }*/

        //dest[j+begin_idx] = source[i];
        dest[i+begin_idx] = source[i];
        //j++;
    }

    return begin_idx+n_elem;
}

void printArray(double* vector, int n_elem)
{
    for(int i=0; i<n_elem; i++)
    {
        printf("%lf ", vector[i]);
    }
    printf("\n");
    fflush(stdout);
}

void double_memset(double* vector, double to_set, int n_elem)
{
    for(int i=0; i<n_elem; i++)
    {
        vector[i] = to_set;
    }
}

//Retornar o indice do k dentro da partição com k
int find_K_idx(int *vector, int k, int n_elem, int* median_holder)
{
    int counter = 0;

    for(int i=0; i<n_elem; i++)
    {
        counter += vector[i];

        if(counter>=k)
        {
            *median_holder = i;
            return k - (counter - vector[i]);
        }
    }
    return -1;
}

int getVectorSum(int* vector, int n_elem)
{
    int counter = 0;

    for(int i=0; i<n_elem; i++)
    {
        counter += vector[i];
    }
    return counter;
}

double medianSort(double* vector, int n_items)
{
    qsort(vector, n_items, sizeof(double), cmpfunc);

    if(n_items % 2 != 0)
    {
        return vector[n_items/2];
    }
    else
    {
        return 0.5 * (vector[n_items/2] + vector[n_items/2 - 1]);
    }
}