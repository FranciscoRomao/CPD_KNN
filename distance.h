#ifndef _DISTANCE_H_
#define _DISTANCE_H_

#include "point.h"

double distance(int dimensions, double* x_coords, double* y_coords);//compute the distance between two points
point* furthest_apart(int numb_points,int dimensions, point* pts, int initial_index, int start_index);//comput the points that are furthest apart from each other

#endif