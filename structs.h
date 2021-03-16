#ifndef _STRUCTS_H_
#define _STRUCTS_H_

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



#endif