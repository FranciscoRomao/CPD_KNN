#include <math.h>
#include <stdlib.h>
#include "distance.h"

double distance(int dimensions, double* x_coords, double* y_coords){
    double distance=0;
    for(int i=0;i<dimensions;i++){
        distance+=pow(x_coords[i]-y_coords[i],2);
    }
    distance=sqrt(distance);
    return distance;
}

point* furthest_apart(int numb_points,int dimensions, point* pts, int initial_index, int start_index){
    point* furthest_pts=(point*)malloc(2*sizeof(point));
    double max_distance=-1;
    for(int i=0;i<numb_points;i++){
        for(int j=i+1;j<numb_points;j++){
            if(distance(dimensions,pts[i].coord,pts[j].coord)>max_distance){
                furthest_pts[0]=pts[i];
                furthest_pts[1]=pts[j];
            }
        }
    }
    return furthest_pts;
}