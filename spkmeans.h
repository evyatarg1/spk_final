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

int get_num_points(FILE *filename);
int get_dim(FILE *filename);
int is_int(char* p);
void input_to_points_struct(FILE *filename, Point** point_arr, int dim, int n);
void points_to_matrix(Point ** point_arr, Matrix * mat);
void to_weighted_adj_mat(Point** point_arr, Matrix* adj_matrix, int dim);
double calc_weight(Point *first, Point *sec, int dim);
void calc_diagonal_degree_mat(Matrix* diag_mat, Matrix * adj_matrix);
void to_sqrt_diag_mat(Matrix* diag_mat);
double calc_row_sum(Matrix * adj_matrix, int row);
void mult_diag_matrix(Matrix * first, Matrix * sec, int is_first_diag, Matrix * res);
void to_l_norm(Matrix *mat);
double off_diag_squares_sum(Matrix * mat);
void finding_pivot(Matrix* a_mat, int* pivot_arr);
void converting_a_and_v_mat(Matrix * a_mat, Matrix * v_mat);
void eigenvalues_into_arr(Matrix * mat, double * arr);
void mat_to_eigen_struct(Matrix * a_eigenvalues, Matrix * v_eigenvectors, Eigen** final_eigen);
int eigengap_heuristic(Eigen** eigen, int n);
int eigen_copm(const void * cp1, const void * cp2);
void sort_eigen_values(Eigen** final_eigen);
void print_mat(Matrix *mat);
void print_transpose_mat(Matrix * mat);
void print_arr(double* arr, int n);
void init_to_identity(Matrix* mat);
void init_mat(Matrix* mat);
void copy_mat_a_to_b(Matrix* mat_a, Matrix* mat_b);
void eigen_struct_to_matrix(Matrix * u_matrix, Eigen **final_eigen);

Matrix* normalize_matrix(Matrix* mat);
double* normalize_row(double* row, int elements_in_row);

void copy_centroid(double* from, double* to, int dim);
void kmeans(double** observations, int* initial_centroid_indices, Cluster** clusters, int k, int n, int dim);
void update_centroids(Cluster** clusters, Point** point_arr, int k, int dim, int n);
int check_difference_in_centroids(Cluster** clusters, int k, int dim);
void print_c(Cluster *cluster, int dim, int cluster_num, int k);
int assign_to_closest_cluster(int point_num, Cluster** clusters, Point** points, int first_insert, int dim, int k);
double find_distance(Point* point, Cluster* cluster, int dim);
void free_matrix(Matrix* matrix);
Matrix* main_logic(int k, char * goal, Point** point_arr, int n, int dim, int flag);
