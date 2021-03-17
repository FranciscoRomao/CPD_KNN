#include <math.h>
#include <stdlib.h>
#include "geometry.h"


typedef struct
{
    int id;
    double *coord;
}point;


typedef struct
{
    int id;
    int lid;
    int rid;
    int radius;
    double *center;
} node;



int cmpfunc(double* a, double* b)
   return (a - b);
   

double calc_Radius(double **projections, int npoints)
{
    double ans = 0;
    double *distances2a;

    distances2a = calc_distances_to_left_limit(projections);

    aux_pt = furthest_point(dimensions, npoints, pts, median);

    void qsort(distances[0], npoints, sizeof(double), cmpfunc);

    if(npoints%2 != 0)
    {
        ans = distances[(npoints-1)/2];
        return ans;

    }
    else
    {
        ans = (distances[npoints/2-1] + distances[npoints/2])/2;
        return ans
    }
}

double distance(int n_dims, double* a, double* b)
{
    double distance=0;
    for(int i=0;i<n_dims;i++){
        distance+=pow(a[i]-b[i],2);
    }
    distance=sqrt(distance);
    return distance;
}

long furthest_point(int n_dims, long np, double **pts, long base_idx)
{
    double max_dist = -1,curr_dist=0;
    long  idx_newpt= 0;
    
    for (long i = 0; i < np; i++)
    {   
        if((curr_dist=distance_(n_dims, pts[base_idx], pts[i]))>max_dist){
            max_dist=curr_dist;
            idx_newpt=i;
        }
    }
    return idx_newpt;
}

/**
 * take two indexes idx_furthest_pts[]={a, b}
 * find the furthest point from a
 * compare the pt found with b
 *  if ==
 *      return both idx, they are the furthest apart
 *  if !=
 *      ifp[1] = a
 *      ifp[0] = newpoint
 *      call again
 * 
 */
long* recursive_furthest_apart(int n_dims, long np, double **pts, long* idx_fp){
    long idx_new_fp = 0;
    idx_new_fp = furthest_point(n_dims, np, **pts, idx_fp[0]);
    if(idx_fp_a != idx_fp[1])
    {
        idx_fp[1] = idx_fp[0];
        idx_fp[0] = idx_new_fp;
        idx_fp = recursive_furthest_apart(n_dims, long np, double **pts, idx_fp)
    }  
    return idx_fp;   
}

long* recursive_furthest_apart(int n_dims,long np,long base_index,double max_distance)
{   
    if(base_index==np-1){//only one point remaining
        return max_distance;
    }
    double new_max_distance=furthest_point(n_dims,np,pts,base_index);
    if(max_distance<new_max_distance) max_distance=new_max_distance;
    return recursive_furthest_apart(n_dims,np,base_index+1);
}