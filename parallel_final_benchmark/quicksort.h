#ifndef _QUICKSORT_PARALLEL_H_
#define _QUICKSORT_PARALLEL_H_

void swap(double *pts, long a, long b);
long partition(double *pts, long low, long high);
void quick_sort_serial(double *pts, long low, long high);
void quick_sort(double *pts, long low, long high, int cutoff);
void quick_sort_parallel(double *pts, int n_points, int cutoff);


#endif