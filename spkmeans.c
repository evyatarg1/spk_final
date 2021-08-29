#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#define EPSILON pow(10,-15)
#define MAX_LOOPS 100

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
}Eigen;

int get_num_points(FILE *filename);
int get_dim(FILE *filename);
int is_int(char* p);
void input_to_points_struct(FILE *filename, Point** point_arr, int dim);
void to_weighted_adj_mat(Point** point_arr, Matrix* adj_matrix, int dim);
double calc_weight(Point *first, Point *sec, int dim);
void calc_diagonal_degree_mat(Matrix* diag_mat, Matrix * adj_matrix);
double calc_row_sum(Matrix * adj_matrix, int row);
void mult_diag_matrix(Matrix * first, Matrix * sec, int is_first_diag, Matrix * res);
double off_diag_squares_sum(Matrix * mat);
void print_mat(Matrix *mat);
int eigengap_heuristic(double* eigenvalues, int n);
void converting_a_and_v_mat(Matrix * a_mat, Matrix * v_mat);
void init_to_identity(Matrix* mat);
void init_mat(Matrix* mat);
void copy_mat_a_to_b(Matrix* mat_a, Matrix* mat_b);
void to_l_norm(Matrix *mat);

Matrix* normalize_matrix(Matrix* mat);
double* normalize_row(double* row, int elements_in_row);


void kmeans_logic(Cluster** clusters, int k, int n, int dim, int max_iter, Matrix* mat);
int assign_to_closest_cluster(int point_num, Cluster** clusters, Point** points, int first_insert, int dim, int k);
double find_distance(Point* point, Cluster* cluster, int dim);
void print_c(Cluster *cluster, int dim);
int update_centroids(Cluster** clusters, Point** point_arr, int k, int dim, int n);
int check_difference_in_centroids(Cluster** clusters, int k, int dim);

int get_num_points(FILE *filename)
{
    int num_of_points;
    int ch;
    int last_ch;
    ch=0;
    last_ch = 0;
    num_of_points = 0;
    ch = fgetc(filename);

    while (ch != EOF){
        if (ch == '\n')
        {
            num_of_points++;
        }
        last_ch = ch;
        ch = fgetc(filename);
    }
    if(last_ch !='\n')
        num_of_points++;
    rewind(filename);
    return num_of_points;
}

int get_dim(FILE *filename){
    int dim;
    int ch;
    ch=0;
    dim = 0;
    ch = fgetc(filename);
    while (ch != '\n'){
        if (ch == ',') dim++;
        ch = fgetc(filename);
    }
    if (dim>0) dim++;
    rewind(filename);
    return dim;

}

int is_int(char* p){
    while (*p != '\0'){
        if (*p<'0' ||*p>'9') return 0;
        p++;
    }
    return 1;
}


void input_to_points_struct(FILE *filename, Point** point_arr, int dim)
{
    int cnt, i;
    char ch;
    double value;
    cnt = 0;
    i=0;

    while (fscanf(filename, "%lf%c", &value, &ch) == 2)
    {
        if (i==0)
        {
            point_arr[cnt] = malloc(sizeof(Point));
            assert(point_arr[cnt] != NULL);
            point_arr[cnt]->coordinates = calloc(dim+1, sizeof(double));
            assert(point_arr[cnt]->coordinates != NULL);

        }
        point_arr[cnt]->coordinates[i] = value;
        i++;
        if (i==dim)
        {
            i=0;
            cnt++;
        }
    }
}

void points_to_matrix(Point ** point_arr, Matrix * mat)
{
    /*Both of the objects must be from the same dim*/
    int i, j;
    for(i = 0; i < mat->rows; i++)
    {
        for (j = 0; j < mat->cols; j++)
        {
            mat->vertices[i][j] = point_arr[i]->coordinates[j];
        }
    }
}

void to_weighted_adj_mat(Point** point_arr, Matrix* adj_matrix, int dim)
{
    int i,j;
    for(i=0; i<adj_matrix->rows;i++)
    {
        for(j=0;j<=i;j++)
        {
            if(i!=j)
            {
                adj_matrix->vertices[i][j] = (double)calc_weight(point_arr[i], point_arr[j], dim);
                adj_matrix->vertices[j][i] = adj_matrix->vertices[i][j];
            }
        }
        adj_matrix->vertices[i][i] = 0;
    }
}


double calc_weight(Point *first, Point *sec, int dim)
{
    int i;
    double distance, res;
    res = 0.0000L;
    for (i=0; i<dim; i++)
    {
        distance = (first->coordinates[i] - sec->coordinates[i]);
        distance = distance * distance;
        res+=distance;
    }
    /*printf("(a1-b1)**2 + ...+ (a.n-b.n)**2= %f\n",res);*/
    res = sqrt(res);
    /*printf("sqrt= %f\n",res);*/
    /*printf("exp= %f\n",exp(res));*/
    res = 1/(exp(res/2));
    /*printf("final= %f\n",res);*/
    return res;
}

void calc_diagonal_degree_mat(Matrix* diag_mat, Matrix * adj_matrix)
{
    /*calc D*/
    int i;
    double row_sum;
    for(i=0; i<diag_mat->rows;i++)
    {
        row_sum = calc_row_sum(adj_matrix, i);
        diag_mat->vertices[i][i] = row_sum;
    }
}

void to_sqrt_diag_mat(Matrix* diag_mat)
{
    int i;
    for(i=0; i<diag_mat->rows;i++)
    {
        if(diag_mat->vertices[i][i]>0)
        {
            diag_mat->vertices[i][i] = 1/sqrt(diag_mat->vertices[i][i]);
        }
    }
}

double calc_row_sum(Matrix * adj_matrix, int row)
{
    int i;
    double res;
    res = 0;
    for(i=0; i<adj_matrix->cols;i++)
    {
        res+= adj_matrix->vertices[row][i];
    }
    return res;
}


void mult_diag_matrix(Matrix * first, Matrix * sec, int is_first_diag, Matrix * res)
{
    int i, j;

    /*One of the matrixex MUST be diagonal! */

    for (i = 0; i < first->rows; i++)
    {
        for (j = 0; j < first->rows; j++)
        {
            if(is_first_diag)
            {
                res->vertices[i][j] += (first->vertices[i][i])*(sec->vertices[i][j]);
            }
            else /*sec mat must be diagonal*/
            {
                res->vertices[i][j] += (first->vertices[i][j])*(sec->vertices[j][j]);

            }
        }
    }
    return;
}


void to_l_norm(Matrix *mat)
{
    int i,j;
    for (i = 0; i < mat->rows; i++)
    {
        for (j = 0; j < mat->rows; j++)
        {
            if(i==j)
            {
                mat->vertices[i][j] = 1-mat->vertices[i][j];

            }
            else
            {
                mat->vertices[i][j] = (-1)*(mat->vertices[i][j]);
            }
        }
    }
    return;
}

double off_diag_squares_sum(Matrix * mat)
{
    int i,j, res;
    res=0;
    for(i=0; i<mat->rows;i++)
    {
        for (j = 0; j < mat->rows; j++)
        {
            if(i!=j)
            {
                res += (mat->vertices[i][j]) * (mat->vertices[i][j]);
            }
        }

    }
    return res;
}

void converting_a_and_v_mat(Matrix * a_mat, Matrix * v_mat)
{
    int i, j, r, cur_row, cur_col, teta_sign, index;
    double t, c, s, teta, max_val;
    double a_ri, a_rj, a_ii, a_jj, a_ij;

    double * row_i_arr;
    double * row_j_arr;
    row_i_arr = (double *) calloc(a_mat->rows, sizeof(double));
    assert(row_i_arr != NULL);
    row_j_arr = (double *) calloc(a_mat->cols, sizeof(double));
    assert(row_j_arr != NULL);

    i=0;
    j=0;
    max_val=0;

    /*finding i,j*/
    for (cur_row = 0; cur_row < a_mat->rows; cur_row++)
    {
        for (cur_col = 0; cur_col < a_mat->cols; cur_col++)
        {
            if(cur_row!=cur_col && (fabs(a_mat->vertices[cur_row][cur_col])>max_val))
            {
                i = cur_row;
                j = cur_col;
                max_val = fabs(a_mat->vertices[cur_row][cur_col]);
            }
        }
    }
    /*
    printf("FINAL: i=%d and j=%d \n",i,j);
    printf("FINAL max val is: %f\n", max_val);
    */

    /* coping j-row and i-col*/
    for(index = 0; index < a_mat->rows; index++)
    {
        row_i_arr[index] = a_mat->vertices[index][i];
        row_j_arr[index] = a_mat->vertices[index][j];
    }

    /*
    for(index = 0; index < a_mat->rows; index++)
    {
        printf("%f ,", row_i_arr[index]);
    }

    printf("\n row j is: \n");
    for(index = 0; index < a_mat->rows; index++)
    {
        printf("%f ,", row_j_arr[index]);
    }
    printf("\n");
    printf("a_mat->vertices[j][j]=%f\n",a_mat->vertices[j][j]);
    printf("a_mat->vertices[i][i]=%f\n", a_mat->vertices[i][i]);
    printf("max_val=%f\n",max_val);
    */

    /*calculating teta, c,t and s*/
    teta = (a_mat->vertices[j][j] - a_mat->vertices[i][i]) / (2*a_mat->vertices[i][j]);

    /*
    printf("teta = (a_mat->vertices[j][j] - a_mat->vertices[i][i]) / (2*max_val) is %f\n", teta);
    */

    if(teta >=0)
    {
        teta_sign = 1;
    }
    else
        teta_sign = -1;

    t = teta_sign / (fabs(teta) + sqrt(teta*teta + 1));
    c = 1 / sqrt(t*t + 1);
    s = t*c;
    /*
      printf("t is: %f\n", t);
      printf("c is: %f\n", c);
      printf("s is: %f\n", s);


  */
    /* Calc eigenvalues*/
    for(r = 0; r<a_mat->rows; r++)
    {
        if(r!=i && r!=j)
        {
            a_ri = c*row_i_arr[r] - s*row_j_arr[r];
            a_rj = c*row_j_arr[r] + s*row_i_arr[r];
            /*
            printf("a_%d%d and a_%d%d are %f\n",i,r, r,i, a_ri);
            printf("a_%d%d and a_%d%d are %f\n",j,r, r,j, a_rj);
            */

            a_mat->vertices[r][i] = a_ri;
            a_mat->vertices[i][r] = a_ri;
            a_mat->vertices[r][j] = a_rj;
            a_mat->vertices[j][r] = a_rj;
        }
    }
    a_ii = (c*c*a_mat->vertices[i][i]) + (s*s*a_mat->vertices[j][j]) - (2*s*c*a_mat->vertices[i][j]);
    a_jj = (s*s*a_mat->vertices[i][i]) + (c*c*a_mat->vertices[j][j]) + (2*s*c*a_mat->vertices[i][j]);
    a_ij = 0;

    a_mat->vertices[i][i] = a_ii;
    a_mat->vertices[j][j] = a_jj;
    a_mat->vertices[i][j] = a_ij;
    a_mat->vertices[j][i] = a_ij;

    /* coping j-col and i-col*/
    for(index = 0; index < a_mat->rows; index++)
    {
        row_i_arr[index] = v_mat->vertices[index][i];
        row_j_arr[index] = v_mat->vertices[index][j];
    }

    /* Calc eigenvalues */
    for (r = 0; r < v_mat->rows; r++)
    {
        v_mat->vertices[r][i] = (c*row_i_arr[r]) - (s*row_j_arr[r]);
        v_mat->vertices[r][j] = (s*row_i_arr[r]) + (c*row_j_arr[r]);
        /*
        printf("v_mat->vertices[index][i]=%f\n", v_mat->vertices[index][i]);
        printf("v_mat->vertices[index][j]=%f\n", v_mat->vertices[index][j]);
        */
    }

    free(row_i_arr);
    free(row_j_arr);
    return;
}

void eigenvalues_into_arr(Matrix * mat, double * arr)
{
    int i;
    for (i = 0; i < mat->rows; i++)
    {
        arr[i] = mat->vertices[i][i];
    }
}

void mat_to_eigen_struct(Matrix * a_eigenvalues, Matrix * v_eigenvectors, Eigen** final_eigen)
{
    int i,j;
    for(i = 0; i<a_eigenvalues->rows; i++)
    {
        final_eigen[i] = malloc(sizeof(Eigen));
        assert(final_eigen[i] != NULL);
        final_eigen[i]->eigen_vector = calloc(a_eigenvalues->rows, sizeof(double));
        assert(final_eigen[i]->eigen_vector != NULL);
    }
    for (i = 0; i < a_eigenvalues->rows; i++)
    {
        for (j = 0; j < a_eigenvalues->rows; j++)
        {
            final_eigen[i]->eigen_vector[j] = v_eigenvectors->vertices[i][j];
        }
        final_eigen[i]->eigen_value =  a_eigenvalues->vertices[i][i];
        final_eigen[i]->n = a_eigenvalues->rows;
    }
}

int eigengap_heuristic(double* eigenvaluse, int n)
{
    int i, max_eigengap, max_eigengap_index, cur_gap;
    max_eigengap = 0;
    max_eigengap_index = 0;
    for (i = 0; i < (int)floor(n/2); i++)
    {
        cur_gap = (eigenvaluse[i] - eigenvaluse[i+1]);
        if(fabs(cur_gap)>max_eigengap)
        {
            max_eigengap = fabs(cur_gap);
            max_eigengap_index = i;
        }
    }
    return max_eigengap_index;
}

void sort_eigen_values(Eigen** final_eigen)
{
    printf("%d...\n", final_eigen[0]->n);
    return;
}

void print_mat(Matrix * mat)
{
    int i, j;
    for (i = 0; i < mat->rows; i++)
    {
        for (j = 0; j < mat->cols; j++)
        {
            /*converting -0.0000 to 0.0000 - all the values that negative and greaten then -0.00005 should be zero */
            if(mat->vertices[i][j] > -0.00005 && mat->vertices[i][j]<0)
            {
                printf("0.0000");
            }
            else
            {
                printf("%.4f", mat->vertices[i][j]);
            }

            if(j!=mat->rows-1)
            {
                printf(",");
            }
        }
        if(i!=mat->rows-1)
        {
            printf("\n");
        }

    }
    /*printf("there are %d cells that are not 0\n", count);*/

}

void init_to_identity(Matrix* mat)
{
    int i;
    mat->vertices = calloc(mat->rows, sizeof(double*));
    for(i=0; i<mat->rows;i++)
    {
        mat->vertices[i] = calloc(mat->rows, sizeof(double));
        assert(mat->vertices[i] != NULL);
        mat->vertices[i][i] = 1;
    }
}

void init_mat(Matrix* mat)
{
    int i;
    mat->vertices = calloc(mat->rows, sizeof(double*));
    for(i=0; i<mat->rows;i++)
    {
        mat->vertices[i] = calloc(mat->cols, sizeof(double));
        assert(mat->vertices[i] != NULL);
    }
}

void copy_mat_a_to_b(Matrix* mat_a, Matrix* mat_b)
{
    int i, j;
    for(i=0; i<mat_a->rows;i++)
    {
        for(j=0; j<mat_b->rows; j++)
        {
            mat_b->vertices[i][j] = mat_a->vertices[i][j];
        }
    }
}

void eigen_struct_to_matrix(Matrix * u_mat, Eigen ** final_eigen, int k)
{
    printf("%d..\n", final_eigen[0]->n);
    printf("%d..\n", u_mat->rows);
    printf("%d...\n", k);
    return;
}



/**step 5 functions**/

Matrix* normalize_matrix(Matrix* mat){
    int i;
    Matrix* normalized_matrix;
    double* normalized_row;

    normalized_matrix = (Matrix*) calloc(1, sizeof(Matrix));
    assert(normalized_row!=NULL);
    normalized_matrix->rows = mat->rows;
    normalized_matrix->cols = mat->cols;
    normalized_matrix->vertices = (double **) calloc(mat->rows, sizeof (double *));
    assert(normalized_matrix->vertices != NULL);
    for (i=0; i<mat->rows; i++){
        normalized_row = normalize_row(mat->vertices[i], mat->cols);
        normalized_matrix->vertices[i] = normalized_row;
    }
    return normalized_matrix;
}

double* normalize_row(double* row, int elements_in_row){
    double sum_of_squares, denominator;
    double* normalized_row;
    int i;

    sum_of_squares=0;
    for (i=0; i<elements_in_row; i++){
        sum_of_squares+= pow(row[i], 2);
    }
    denominator = pow(sum_of_squares, 0.5);

    normalized_row = (double *) calloc(elements_in_row, sizeof (double ));
    assert(normalized_row!=NULL);
    for (i=0; i<elements_in_row; i++){
        normalized_row[i] = row[i]/denominator;
    }
    return normalized_row;
}

/**step 6 functions**/

void copy_centroid(double* from, double* to, int dim){
    int i;
    for (i=0;i<dim; i++){
        to[i] = from[i];
    }
}

void kmeans_logic(Cluster** clusters, int k, int n, int dim, int max_iter, Matrix* mat){
    int i, count, q;
    int differ;
    Point** point_arr;
    differ =1;
    i=0;
    point_arr = (Point**) calloc(n, sizeof(Point*));
    assert(point_arr != NULL);

    for (i=0; i< n; i++){
        point_arr[i] = (Point*) calloc(1, sizeof (Point));
        assert(point_arr[i]!=NULL);
        point_arr[i]->coordinates = mat->vertices[i];
        if (i<k){
            clusters[i] = calloc(1, sizeof (Point));
            assert(clusters[i]!=NULL);
            clusters[i]->curr_centroid = (double *) calloc(dim+1, sizeof (double ));
            assert(clusters[i]->curr_centroid!=NULL);
            copy_centroid(mat->vertices[i], clusters[i]->curr_centroid, dim);
            clusters[i]->prev_centroid = (double *) calloc(dim+1, sizeof (double ));
            assert(clusters[i]->prev_centroid!=NULL);
            copy_centroid(mat->vertices[i], clusters[i]->prev_centroid, dim);
            clusters[i]->size++;
            point_arr[i]->cluster = i;
        }
    }
    for (i=k; i<n; i++){
        assign_to_closest_cluster(i, clusters, point_arr, 1, dim, k);
    }

    count = 1;

    while((count < max_iter) && differ)
    {
        for(q=0; q<n; q++)
        {
            assign_to_closest_cluster(q, clusters, point_arr, 0, dim, k);
        }

        update_centroids(clusters, point_arr, k, dim, n);
        count++;
        differ = check_difference_in_centroids(clusters, k, dim);

    }
    for (i=0; i<n; i++){
        free(point_arr[i]->coordinates);
    }
    free(point_arr);

}

/*This function gets a cluster 2D array, a point array, and k+dim. Its return a new Cluster 2D object
 where the curr_centroid array is updated  */
int update_centroids(Cluster** clusters, Point** point_arr, int k, int dim, int n)
{
    int i;
    int j;
    int t;
    int first_loop;
    first_loop =1;

    /*for each cluster*/
    for(i=0; i<k; i++)
    {
        /*for each point*/
        for(j=0; j<n; j++)
        {
            /* if this point in the cluster*/
            if(point_arr[j]->cluster == i)
            {
                /* each coordinate */
                for(t=0; t<dim; t++)
                {
                    if(first_loop)
                    {
                        clusters[i]->prev_centroid[t] = clusters[i]->curr_centroid[t];
                        clusters[i]->curr_centroid[t] = 0.0000L;
                    }
                    clusters[i]->curr_centroid[t] += point_arr[j]->coordinates[t];

                }
                first_loop = 0;
            }
        }
        /* after we sum the coordinated, we will divied each one of them in the sum*/
        if(clusters[i]->size > 1)
        {
            for(t=0; t<dim; t++)
            {
                clusters[i]->curr_centroid[t] = (clusters[i]->curr_centroid[t]) / (clusters[i]->size);

            }
        }
        first_loop = 1;
    }
    return 1;
}

/*This function gets a cluster (and k+dim). Its returns:
1 - if there is a any difference between the "cerrent" and "prev" of any cluster
0 - otherwise   */
int check_difference_in_centroids(Cluster** clusters, int k, int dim)
{
    int i, j;
    /*for each cluster*/
    for(i=0; i<k; i++)
    {
        /*for each coordinate in the centroid */
        for(j=0; j<dim; j++)
        {
            if(clusters[i]->curr_centroid[j]!=clusters[i]->prev_centroid[j])
            {
                return 1; /* we found a difference */
            }
        }
    }

    return 0; /* all of the centroides are the same*/
}

int assign_to_closest_cluster(int point_num, Cluster** clusters, Point** points, int first_insert, int dim, int k)
{
    int i, index, prev_index;
    double min_distance, curr_distance;
    if(first_insert)
    {
        min_distance = find_distance(points[point_num], clusters[0], dim);
        index = 0;
    }
    else
    {
        index = points[point_num]->cluster;
        min_distance = find_distance(points[point_num], clusters[points[point_num]->cluster], dim);
    }
    /*for all k clusters in*/
    for (i=0; i<k; i++)
    {
        curr_distance = find_distance(points[point_num], clusters[i], dim);
        if (curr_distance<min_distance)
        {
            min_distance = curr_distance;
            index = i; /* change point index to this cluster*/
        }
    }
    /* need to add the point to relevant cluster only */
    if (first_insert)
    {
        points[point_num]->cluster = index;
        clusters[index]->size+=1;
    }

        /* need to add the point to relevant cluster + delete it from previos one*/
    else
    {
        /*remove from old cluster*/
        prev_index = points[point_num]->cluster;
        clusters[prev_index]->size--;
        /*add to new cluster*/
        points[point_num]->cluster = index;
        clusters[index]->size++;
    }
    return 0;
}


double find_distance(Point* point, Cluster* cluster, int dim)
{
    int i;
    double distance, res;
    res = 0.0000L;
    for (i=0; i<dim; i++){
        distance = (point->coordinates[i] - cluster->curr_centroid[i]);
        distance = distance * distance;
        res+=distance;
    }
    return res;
}

void print_c(Cluster *cluster, int dim)
{
    int i;
    for (i=0; i<dim; ++i)
    {
        printf("%.4f",cluster->curr_centroid[i]);
        if(i!=(dim-1))
        {
            printf(",");
        }
    }
}


int kmeans(Matrix* mat, int k){

    int max_iter, n, dim, i;
    Cluster** clusters;

    max_iter = 300;

    n = mat->rows;
    dim = mat->cols;
    clusters = (Cluster**) calloc(k, sizeof(Cluster*));
    assert(clusters != NULL);

    kmeans_logic(clusters, k, n, dim, max_iter, mat);

    for (i=0; i<k; i++)
    {
        print_c(clusters[i], dim);
        if (i!=k)
            printf("\n");
    }

    for (i=0; i<k; i++){
        free(clusters[i]->curr_centroid);
        free(clusters[i]->prev_centroid);
    }

    free(clusters);


    return 1;
}


int main_logic(int k, char * goal, char * f_name, int flag)
{
    /*  ---- Declaration ----  */
    int  n, dim, i, count;
    FILE *filename;
    Matrix* adj_matrix;
    Matrix* diag_degree_mat;
    Matrix* sqrt_diag_degree_mat;
    Matrix* l_norm_mat;
    Matrix* temp_mat;

    Matrix* a_eigenvalues;
    Matrix* v_eigenvectors;
    Point** point_arr;
    Eigen** final_eigen;
    Matrix* u_matrix;
    double* arr_eigenvalues;
    double prev_off_a, cur_off_a;

    /*  ---- Basic Validation ----  */
    filename = fopen(f_name, "r");
    if(filename==NULL){
        printf("Invalid Input!\n");
        exit(1);
    }

    /* ---- Calculation of input dimentions (n, dim) ----*/
    n = get_num_points(filename);
    dim = get_dim(filename);
    if (n == 0 || dim == 0 || k>=n){
        printf("Invalid Input!\n");
        fclose(filename);
        exit(1);
    }

    /*  ---- Input to struct ----  */
    point_arr = (Point**) calloc(n, sizeof(Point*));
    assert(point_arr != NULL);
    input_to_points_struct(filename, point_arr, dim);

    /* ######## Calculate all relevant matrixes ########*/


    /*---- Weighted Adjacency Matrix (W) ----*/
    adj_matrix = (Matrix*) calloc(n, sizeof(double*));
    assert(adj_matrix != NULL);
    adj_matrix->rows = n;
    adj_matrix->cols = n;
    init_mat(adj_matrix);
    to_weighted_adj_mat(point_arr, adj_matrix, dim);

    /*---- Diagonal Degree Matrix (D) ----*/
    diag_degree_mat = (Matrix*) calloc(n, sizeof(double*));
    assert(diag_degree_mat != NULL);
    diag_degree_mat->rows = n;
    diag_degree_mat->cols = n;
    init_mat(diag_degree_mat);
    calc_diagonal_degree_mat(diag_degree_mat, adj_matrix);

    /*----  Making D into D^(-1/2) ----*/
    sqrt_diag_degree_mat = (Matrix*) calloc(n, sizeof(double*));
    assert(sqrt_diag_degree_mat != NULL);
    sqrt_diag_degree_mat->rows = n;
    sqrt_diag_degree_mat->cols = n;
    init_mat(sqrt_diag_degree_mat);
    copy_mat_a_to_b(diag_degree_mat, sqrt_diag_degree_mat); /*coping D */
    to_sqrt_diag_mat(sqrt_diag_degree_mat); /*changing the matrix to D^(-1/2)*/

    /*---- L Norm Matrix (lnorm)----*/
    l_norm_mat = (Matrix*) calloc(n, sizeof(double*));
    assert(l_norm_mat != NULL);
    l_norm_mat->rows = n;
    l_norm_mat->cols = n;
    init_mat(l_norm_mat);

    /* init additional temp matrix (for calculation purposes) */
    temp_mat = (Matrix*) calloc(n, sizeof(double*));
    assert(temp_mat != NULL);
    temp_mat->rows = n;
    temp_mat->cols = n;
    init_mat(temp_mat);

    /*---- Calc lnorm matrix ----*/
    mult_diag_matrix(sqrt_diag_degree_mat, adj_matrix, 1, temp_mat); /*   D^(-0.5) x (W)  */
    mult_diag_matrix(temp_mat, sqrt_diag_degree_mat, 0, l_norm_mat); /*  (D^(-0.5) * W)  x  (D^(-0.5)) */
    to_l_norm(l_norm_mat);     /*   I - (A*B)*A        */




    if(strcmp(goal,"spk")==0 || strcmp(goal, "jacobi")==0)
    {
        /*---- Calc A matrix (eiganvalues) and V matrix (eiganvectors) ----*/
        /* init eigenvalues */
        a_eigenvalues = (Matrix*) calloc(n, sizeof(double*));
        assert(a_eigenvalues != NULL);
        a_eigenvalues->rows = n;
        a_eigenvalues->cols = n;
        init_mat(a_eigenvalues);
        if(strcmp(goal, "jacobi")==0)
        {
            /*takes from input*/
            points_to_matrix(point_arr, a_eigenvalues);
        }
        else
        {
            /*takes from l norm*/
            copy_mat_a_to_b(l_norm_mat, a_eigenvalues);

        }

        /* init eigenvectors */
        v_eigenvectors = (Matrix*) calloc(n, sizeof(double*));
        assert(v_eigenvectors != NULL);
        v_eigenvectors->rows = n;
        v_eigenvectors->cols = n;
        init_to_identity(v_eigenvectors);

        /*----  Calc the eigenvectors and eigenvalues  ----*/
        prev_off_a = 0;
        cur_off_a = off_diag_squares_sum(l_norm_mat);
        count = 0;
        while (count<MAX_LOOPS && (prev_off_a-cur_off_a) <= EPSILON) /*until Convergence or Max loops*/
        {
            prev_off_a = cur_off_a;
            count++;
            converting_a_and_v_mat(a_eigenvalues, v_eigenvectors);
            cur_off_a = off_diag_squares_sum(a_eigenvalues);
        }


        arr_eigenvalues = (double*) calloc(n, sizeof(double));
        assert(arr_eigenvalues != NULL);
        eigenvalues_into_arr(a_eigenvalues, arr_eigenvalues);

        /* creating Eigen struct and coping the values of the relevant matrices to it*/
        final_eigen = (Eigen**) calloc(n, sizeof(Eigen*));
        assert(final_eigen != NULL);


        /*############  need to complete from here #############*/
        mat_to_eigen_struct(a_eigenvalues, v_eigenvectors, final_eigen);
        sort_eigen_values(final_eigen); /*sorting the eigenvalues for creating U anf for k_heuristic*/
        /*############  up to here #############*/


        if(k==0 && strcmp(goal, "spk")==0)
        {
            /*############  need to change the input - from double* to Eigen** #############*/
            k = eigengap_heuristic(arr_eigenvalues, n);
        }

        u_matrix = (Matrix*) calloc(n, sizeof(double*));
        assert(u_matrix != NULL);
        u_matrix->rows = n;
        u_matrix->cols = k;
        init_mat(u_matrix);

        /*############  need to complete from here #############*/
        eigen_struct_to_matrix(u_matrix, final_eigen, k);
        /*############  up to here #############*/


        if(strcmp(goal, "spk")==0)
        {
            if(flag) /* 1 is python, 0 is C*/
            {
                /*python - T goes to ex2*/
            }
            else
            {
                /*C - T goes to ex1*/
            }
            printf("In construction :)\n");
        }
        if(strcmp(goal, "jacobi")==0)
        {
            /* print the eigen values*/
            for (i = 0; i < n; i++)
            {
                if(arr_eigenvalues[i]<0 && arr_eigenvalues[i] > -0.00005)
                {
                    printf("0.0000");
                }
                else
                {
                    printf("%.4f", arr_eigenvalues[i]);
                }
                if(i+1!=n)
                {
                    printf(", ");
                }
            }
            printf("\n");
            print_mat(v_eigenvectors);
        }
    }
    else if(strcmp(goal,"wam")==0)
    {
        print_mat(adj_matrix);
    }
    else if(strcmp(goal,"ddg")==0)
    {
        print_mat(diag_degree_mat);
    }
    else if(strcmp(goal,"lnorm")==0)
    {
        print_mat(l_norm_mat);
    }
    else
    {
        printf("Invalid Input!\n");
        exit(1);
    }

    /* ######### FREE ######### */
    for(i=0; i<n; i++)
    {
        free(point_arr[i]->coordinates);
        free(point_arr[i]);

        free(adj_matrix->vertices[i]);
        free(diag_degree_mat->vertices[i]);
        free(sqrt_diag_degree_mat->vertices[i]);
        free(l_norm_mat->vertices[i]);
        free(temp_mat->vertices[i]);
        if(strcmp(goal, "spk")==0 || strcmp(goal, "spk")==0)
        {
            free(a_eigenvalues->vertices[i]);
            free(v_eigenvectors->vertices[i]);
            free(final_eigen[i]->eigen_vector);
            free(final_eigen[i]);
        }
    }

    free(point_arr);
    free(adj_matrix);
    free(diag_degree_mat);
    free(sqrt_diag_degree_mat);
    free(l_norm_mat);
    free(temp_mat);
    if(strcmp(goal, "spk")== 0 || strcmp(goal, "spk")==0)
    {
        free(a_eigenvalues);
        free(v_eigenvectors);
        free(final_eigen);
        free(arr_eigenvalues);
    }
    fclose(filename);
    return 1;
}

int main(int argc, char** argv)
{
    /*CMD: prog_name, k, goal, file.mane*/

    int k;
    char* f_name;
    char* goal;

    if (argc != 4 || (!is_int(argv[1])))
    {
        printf("Invalid Input!\n");
        exit(1);
    }

    k = atoi(argv[1]);
    goal = argv[2];
    f_name = argv[3];

    main_logic(k, goal, f_name, 0);
    return 1;
}
