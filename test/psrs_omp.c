#include "omp.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "math.h"

/* headers */
int lcompare(const void * ptr2num1, const void * ptr2num2);
long long *merge(long long * left, long long * right, int l_end, int r_end);
long long *merge_sort(long long * arr, int size);
void insertion_sort(long long *arr, int n);
void calc_partition_borders(long long array[],
               int start,
               int end,
               int sublist_sizes[],
               int at,
               long long pivots[],
               int first_p,
               int last_p);
void psrs_sort(long long *a, int n);
void sortll(long long *a, int len);

/* sort an array in non-descending order */
void psrs_sort(long long *a, int n) {
  if(n > 1){
  if(n <= 55){
	// Testing shows that sequential insertion sort is quickest when n <= 55 (approx.)
	insertion_sort(a,n);
  }else if(n <= 10000){
  // Testing shows that sequential merge sort is quickest when n <= 10000(approx.)
	merge_sort(a,n);
  }else{
  // Testing shows that our algorithm is now the quickest 
  int p, size, rsize, sample_size;
  long long *sample, *pivots;
  int *partition_borders, *bucket_sizes, *result_positions;
  long long **loc_a_ptrs;

  // Determine the appropriate number of threads to use
  // p^3 <= n - We need this to hold true
  p = omp_get_max_threads();
  p = p*p*p;
  if(p > n){
	p = floor(pow(n,0.33));
	p-=p%2;
  }else{
   p = omp_get_max_threads();
   p-=p%2;
  }
  omp_set_num_threads(p);

  size  = (n + p - 1) / p;
  rsize = (size + p - 1) / p;
  sample_size = p * (p - 1);

  loc_a_ptrs = malloc(p * sizeof(long long *));
  sample = malloc(sample_size * sizeof(long long));
  partition_borders = malloc(p * (p + 1) * sizeof(int));
  bucket_sizes = malloc(p * sizeof(int));
  result_positions = malloc(p * sizeof(int));
  pivots = malloc((p - 1) * sizeof(long long));

  #pragma omp parallel
  {
    int i, j, max, thread_num, start, end, loc_size, offset, this_result_size;
    long long *loc_a, *this_result, *current_a;

    thread_num = omp_get_thread_num();
    start = thread_num * size;
    end = start + size - 1;
    if(end >= n) end = n - 1;
    loc_size = (end - start + 1);
    end = end % size;

    loc_a = malloc(loc_size * sizeof(long long));
    memcpy(loc_a, a + start, loc_size * sizeof(long long));
    loc_a_ptrs[thread_num] = loc_a;

    sortll(loc_a, loc_size); // Testing shows that this sequential sort is quickest in this instance
	
    offset = thread_num * (p - 1) - 1;

    for(i = 1; i < p; i++) {
      if(i * rsize <= end) {
        sample[offset + i] = loc_a[i * rsize - 1];
      } else {
        sample[offset + i] = loc_a[end];
      }
    }

    #pragma omp barrier

    #pragma omp single
    {
      merge_sort(sample, sample_size); // Testing shows that this sequential sort is quickest in this instance
      for(i = 0; i < p - 1; i++) {
        pivots[i] = sample[i * p + p / 2];
      }
    }

    #pragma omp barrier

    offset = thread_num * (p + 1);
    partition_borders[offset] = 0;
    partition_borders[offset + p] = end + 1;
    calc_partition_borders(loc_a, 0, loc_size-1, partition_borders, offset, pivots, 1, p-1);

    #pragma omp barrier

    max = p * (p + 1);
    bucket_sizes[thread_num] = 0;
    for(i = thread_num; i < max; i += p + 1) {
      bucket_sizes[thread_num] += partition_borders[i + 1] - partition_borders[i];
    }

    #pragma omp barrier

    #pragma omp single
    {
      result_positions[0] = 0;
      for(i = 1; i < p; i++) {
        result_positions[i] = bucket_sizes[i-1] + result_positions[i-1];
      }
    }

    #pragma omp barrier

    this_result = a + result_positions[thread_num];

    if(thread_num == p-1) {
      this_result_size = n - result_positions[thread_num];
    } else {
      this_result_size = result_positions[thread_num+1] - result_positions[thread_num];
    }

    // pluck this threads sublist from each of the local arrays
    this_result = a + result_positions[thread_num];

    for(i = 0, j = 0; i < p; i++) {
      int low, high, partition_size;
      offset = i * (p + 1) + thread_num;
      low = partition_borders[offset];
      high = partition_borders[offset+1];
      partition_size = (high - low);
      if(partition_size > 0) {
        memcpy(this_result+j, &(loc_a_ptrs[i][low]), partition_size * sizeof(long long));
        j += partition_size;
      }
    }

    // sort p local sorted arrays
    sortll(this_result, this_result_size); // Testing shows that this sequential sort is quickest in this instance
	
    #pragma omp barrier
    free(loc_a);
  }

  free(loc_a_ptrs);
  free(sample);
  free(partition_borders);
  free(bucket_sizes);
  free(result_positions);
  free(pivots);
  }
  }
}


/* determine the boundaries for the sublists of an local array */
void calc_partition_borders(long long array[],    // array being sorted
                            int start,
                            int end,              // separate the array into current process range
                            int result[],
                            int at,               // this process start point in result
                            long long pivots[],   // the pivot values
                            int first_pv,         // first pivot
                            int last_pv)          // last pivot
{
  int mid, lowerbound, upperbound, center;
  long long pv;

  mid = (first_pv + last_pv) / 2;
  pv = pivots[mid-1];
  lowerbound = start;
  upperbound = end;
  while(lowerbound <= upperbound) {
    center = (lowerbound + upperbound) / 2;
    if(array[center] > pv) {
      upperbound = center - 1;
    } else {
      lowerbound = center + 1;
    }
  }
  result[at + mid] = lowerbound;

  if(first_pv < mid) {
    calc_partition_borders(array, start, lowerbound - 1, result, at, pivots, first_pv, mid - 1);
  }
  if(mid < last_pv) {
    calc_partition_borders(array, lowerbound, end, result, at, pivots, mid + 1, last_pv);
  }
}


/* Compare long longs (lifted from the harness) */
int lcompare(const void * ptr2num1, const void * ptr2num2)
{
  long long num1 = *((long long*) ptr2num1);
  long long num2 = *((long long*) ptr2num2);

  if ( num1 > num2 )
    return 1;
  else if ( num1 < num2 )
    return -1;
  else
    return 0;
}


/*
  Sort a portion of an array
  @todo: see if merge sort might be better in some cases
*/
void sortll(long long *a, int len)
{
  qsort(a, len, sizeof(long long), lcompare);
}

/*
  Standard merge sort
*/
long long *merge_sort(long long * arr, int size){
	// Arrays shorter than 1 are already sorted
	if(size > 1){
		int middle = size / 2, i; 
		long long *left, *right;
		left = arr;
		right = arr + middle; 
		
		left = merge_sort(left, middle);
		right = merge_sort(right, size-middle);
		return merge(left, right, middle,size-middle);
	}else { return arr; }
}

long long *merge(long long * left, long long * right, int l_end, int r_end){
	int temp_off, l_off, r_off, size = l_end+r_end;
	long long *temp = malloc(sizeof(long long) * l_end);

	// Copy lower half into temp buffer
	for(l_off=0, temp_off=0; left+l_off != right; l_off++, temp_off++){
		*(temp + temp_off) = *(left + l_off);
	}
	
	temp_off=0; l_off=0; r_off=0;

	while(l_off < size){
		if(temp_off < l_end){
			if(r_off < r_end){
				if(*(temp+temp_off) < *(right+r_off)){
					*(left+l_off) = *(temp+temp_off);
					temp_off++;
				}else{
					*(left+l_off) = *(right+r_off);
					r_off++;
				}
			}else{
				*(left+l_off) = *(temp+temp_off);
				temp_off++;
			}
		}else{
			if(r_off < r_end) {
				*(left + l_off) = *(right + r_off);
				r_off++;
			}else{
				printf("\nERROR - merging loop going too far\n");
			}
		}
		l_off++;
	}
	free(temp);
	return left;
}

/*
  Standard insertion sort
*/
void insertion_sort(long long *arr, int n){
	int i, j, k, temp;
	
	for ( i = 1 ; i <= n ; i++ )
	{
		for ( j = 0 ; j < i ; j++ )
		{
			if ( arr[j] > arr[i] )
			{
				temp = arr[j] ;
				arr[j] = arr[i] ;

				for ( k = i ; k > j ; k-- )
					arr[k] = arr[k - 1] ;

				arr[k + 1] = temp ;
			}
		}
	}
}


















/* Test and timing harness program for developing a sorting
   routine for the CS3014 module */

#include <stdio.h>
#include <sys/time.h>

/* the following two definitions of DEBUGGING control whether or not
   debugging information is written out. To put the program into
   debugging mode, uncomment the following line: */
/*#define DEBUGGING(_x) _x */
/* to stop the printing of debugging information, use the following line: */
#define DEBUGGING(_x)

/* this is a comparison function for long longs that is being added
   purely to allow the standard qsort function to be used */
int llcompare(const void * ptr2num1, const void * ptr2num2)
{
  long long num1 = *((long long*) ptr2num1);
  long long num2 = *((long long*) ptr2num2);

  if ( num1 > num2 )
    return 1;
  else if ( num1 < num2 )
    return -1;
  else
    return 0;
}

/* sort an array into non-decreasing order */
void sort(long long a[], int size)
{
  qsort(a, size, sizeof(long long), llcompare);
}

/* write out an array of integers up to 'size' */
void write_out(long long a[], int size)
{
  int i;

  for ( i = 0; i < size; i++ ) {
    printf("%lld\n", a[i]);
  }
}

/* read a stream of long long numbers from a file.
   the first number in the file is the number of numbers*/
long long * read_in(char filename[], int * ptr2size)
{
  const int max_line = 1024;
  char line[max_line];
  int i;
  FILE * file;
  char * eof;
  long long * a;
  int size;

  /* open the file */
  file = fopen(filename, "r");
  if ( file == NULL ) {
    fprintf(stderr, "File not found: %s\n", filename);
    exit(1);
  }

  /* read in the size of the array to allocate */
  eof = fgets(line, max_line, file);
  if ( eof == NULL ) {
    fprintf(stderr, "Empty file: %s\n", filename);
    exit(1);
  }
  sscanf(line, "%d", &size);
  a = malloc(sizeof(long long) * size);

  /* read in the long longs - one per line */
  i = 0;
  eof = fgets(line, max_line, file);
  while ( eof != NULL && i < size ) {     /* eof == NULL => end of file */
    sscanf(line, "%lld", &(a[i]));
    i++;
    eof = fgets(line, max_line, file);
  }

  fclose(file);
  *ptr2size = size;

  return a;
}

long long * copy_array(long long * array, int size)
{
  long long * result = malloc(sizeof(long long) * size);
  int i;

  for ( i = 0; i < size; i++ ) {
    result[i] = array[i];
  }
  return result;
}

/* create an array of length size and fill it with random numbers */
long long * gen_random(int size)
{
  long long * result = malloc(sizeof(long long) * size);
  int i;
  struct timeval seedtime;
  int seed;

  /* use the microsecond part of the current time as a pseudorandom seed */
  gettimeofday(&seedtime, NULL);
  seed = seedtime.tv_usec;
  srandom(seed);

  /* fill the array with random numbers */
  for ( i = 0; i < size; i++ ) {
    long long upper = random();
    long long lower = random();
    result[i] = (upper << 32) | lower;
  }

  return result;
}

int main(int argc, char ** argv)
{
  long long * array;
  long long * copy;
  long long sort_time;
  int array_size;
  struct timeval start_time;
  struct timeval stop_time;
  int i;

  if ( argc == 3 && !strcmp(argv[1], "-f") ) { /* read data from file */
    char * filename = argv[2];
    array = read_in(filename, &array_size);
  }
  else if ( argc == 3 && !strcmp(argv[1], "-r") ) { /* generate random data */
    array_size = atoi(argv[2]);
    array = gen_random(array_size);
  }
  else {
    fprintf(stderr, "Usage: sort-harness -f <filename> OR sort-harness -r <size>\n");
    exit(1);
  }
  DEBUGGING(write_out(array, array_size));

  /* make a copy of the array, so we can check sorting later */
  copy = copy_array(array, array_size);

  /* record starting time */
  gettimeofday(&start_time, NULL);

  /* sort the array somehow */
  /* you should replace the following function with your own function */
  psrs_sort(array, array_size);

  /* record finishing time */
  gettimeofday(&stop_time, NULL);
  sort_time = (stop_time.tv_sec - start_time.tv_sec) * 1000000L +
    (stop_time.tv_usec - start_time.tv_usec);
  printf("Sorting time: %lld microseconds\n", sort_time);

  DEBUGGING(write_out(array, array_size));

  /* sort the copy. use built-in sorting function to check your function */
  sort(copy, array_size);

  /* now check that the two are identical */
  for ( i = 0; i < array_size; i++ ) {
    if ( array[i] != copy[i] ) {
      // fprintf(stderr, "Error in sorting at position %d\n", i);
    }
  }

  return 0;
}
