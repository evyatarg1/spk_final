#include <stdio.h>
#include <stdlib.h>

#define EPSILON pow(10,-15)
#define MAX_LOOPS 100
#define MAX_ITER 300

typedef struct Cluster{
    int size;
    double* prev_centroid;
    double* curr_centroid;
} Cluster;

typedef struct Point{
    int cluster;
    double* coordinates;
} Point;

typedef struct Matrix{
    int rows;
    int cols;
    double** vertices;
} Matrix;

typedef struct Eigen{
    int n;
    double * eigen_vector;
    double eigen_value;
    int origin_index;
}Eigen;


void kmeans(double** observations, int* initial_centroid_indices, Cluster** clusters, int k, int n, int dim);
Matrix* main_logic(int k, char * goal, Point** point_arr, int n, int dim, int flag);
void free_matrix(Matrix* matrix);
