#include "spkmeans.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

/*This function gets a file object, and returns the number of points (non-empty rows) in it*/
int get_num_points(FILE *filename)
{
    int last_ch;
    int num_of_points;
    int ch;
    last_ch = 0;
    ch=0;
    num_of_points = 0;
    ch = fgetc(filename);

    while (ch != EOF)
    {
        if (ch == '\n' && last_ch!='\n')
        {
            num_of_points++;
        }
        last_ch = ch;
        ch = fgetc(filename);
    }
    if(last_ch!='\n')
    {
        /*we need to count the last point*/
        num_of_points++;
    }
    rewind(filename);
    return num_of_points;
}

/*This function gets a file obgect, and returns the dim of its points (asuuming each point has the same dim)*/
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

/*This function gets a char, and checks if its int (retruns 1) else 0*/
int is_int(char* p){
    while (*p != '\0'){
        if (*p<'0' ||*p>'9') return 0;
        p++;
    }
    return 1;
}

/*This function get a file object that contain a list of points, their dim and the amount of points (n),and insert the points into a Point** object*/
void input_to_points_struct(FILE *filename, Point** point_arr, int dim, int n)
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
            if(point_arr[cnt]==NULL)
                printf("An Error Has Occured");
            assert(point_arr[cnt] != NULL);
            point_arr[cnt]->coordinates = calloc(dim+1, sizeof(double));
            if(point_arr[cnt]->coordinates==NULL)
                printf("An Error Has Occured");
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
    /*in case that there is no '\n' in the end of the file,
    and therefore - we will need to init it manualy*/
    point_arr[n-1]->coordinates[dim-1] = value;
}

/*This function extracting the points from Point object, and insert them to a Matrix object (Both of the objects must be from the same dim)*/
void points_to_matrix(Point ** point_arr, Matrix * mat)
{
    int i, j;
    for(i = 0; i < mat->rows; i++)
    {
        for (j = 0; j < mat->cols; j++)
        {
            mat->vertices[i][j] = point_arr[i]->coordinates[j];
        }
    }
}

/*This function calculating the Weighted Adjacency Matrix (gets a point object, and insert the results into a matrix object)*/
void to_weighted_adj_mat(Point** point_arr, Matrix* adj_matrix, int dim)
{
    int i,j;
    for(i=0; i<adj_matrix->rows;i++)
    {
        for(j=0;j<=i;j++)
        {
            if(i!=j)
            {
                adj_matrix->vertices[i][j] = (double) calc_weight(point_arr[i], point_arr[j], dim);
                adj_matrix->vertices[j][i] = adj_matrix->vertices[i][j];
            }
        }
        adj_matrix->vertices[i][i] = 0;
    }
}


/*This function gets 2 poionts, and calculate the weight of their "Edge" */
double calc_weight(Point *first, Point *sec, int dim)
{
    int i;
    double distance, res;
    res = 0.00000000L;
    distance = 0.00000000L;
    for (i=0; i<dim; i++)
    {
        distance = (first->coordinates[i] - sec->coordinates[i]);
        distance = distance * distance;
        res+=distance;
    }
    res = sqrt(res);
    res = 1/(exp(res/2));
    return res;
}

/*This dunction calculate the Diagonal Degree Matrix*/
void calc_diagonal_degree_mat(Matrix* diag_mat, Matrix * adj_matrix)
{
    int i;
    double row_sum;
    for(i=0; i<diag_mat->rows;i++)
    {
        row_sum = calc_row_sum(adj_matrix, i);
        diag_mat->vertices[i][i] = row_sum;
    }
}

/*On given matrix (D), this function will calculate D^(-1/2)*/
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

/*This function gets a amtrix and a row index, and returns the sum of the cells in the given row in the matrix*/
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

/*This functions gets 2 matrixes, and mult them, when one of the matrixes must be diagonal (is_first_diag =1 if "first" is diagonal, is_first_diag=0 if "sec" is diagonal). the results will be insert into "res" matirx*/
void mult_diag_matrix(Matrix * first, Matrix * sec, int is_first_diag, Matrix * res)
{
    int i, j;
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

/*This function convert a given matrix [(D^−1/2)*(W)*(D^-1/2)]  into l norm.
It sub I (identity mat) with a the given mat */
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

/*This function gets a matrix, and calculate the sum of squares of all off-diagonal elements of the mat */
double off_diag_squares_sum(Matrix * mat)
{
    int i,j;
    double res;
    res=0;
    for(i=0; i<mat->rows;i++)
    {
        for (j = 0; j < mat->cols; j++)
        {
            if(i!=j)
            {
                res += (mat->vertices[i][j]) * (mat->vertices[i][j]);
            }
        }
    }
    return res;
}

/*This functions gets an eigenvalues mat, and search for the off-diagonal element
in it with the largest absolute value. Its insert the 2 indices to the given array*/
void finding_pivot(Matrix* a_mat, int* pivot_arr)
{
    int cur_row, cur_col;
    double max_val;
    pivot_arr[0] = 0; /*i*/
    pivot_arr[1] = 0; /*j*/
    max_val=0;

    /*finding i,j*/
    for (cur_row = 0; cur_row < a_mat->rows; cur_row++)
    {
        for (cur_col = 0; cur_col < a_mat->cols; cur_col++)
        {
            if(cur_row<cur_col)
            {
                if(fabs(a_mat->vertices[cur_row][cur_col])==max_val)
                {
                    if (cur_row<pivot_arr[0] || (cur_row==pivot_arr[0] && cur_col<pivot_arr[1]))
                    {
                        pivot_arr[0] = cur_row;
                        pivot_arr[1] = cur_col;
                        max_val = fabs(a_mat->vertices[pivot_arr[0]][pivot_arr[1]]);
                    }
                }
                else if(fabs(a_mat->vertices[cur_row][cur_col])>max_val)
                {
                    pivot_arr[0] = cur_row;
                    pivot_arr[1] = cur_col;
                    max_val = fabs(a_mat->vertices[pivot_arr[0]][pivot_arr[1]]);
                }
            }
        }
    }
}

/*This function gets an eigenvalues matrix (a_mat), and its eigenvectors matrix (v_mat), and calculates one step of the jacobi algorithm simultaneously. The results will be stored in a_mat and v_mat, respectively*/
void converting_a_and_v_mat(Matrix * a_mat, Matrix * v_mat)
{
    int i, j, r, teta_sign, index;
    double t, c, s, teta, a_ri, a_rj, a_ii, a_jj;

    int* pivot_arr;
    double * row_i_arr;
    double * row_j_arr;
    row_i_arr = (double *) calloc(a_mat->rows, sizeof(double));
    if(row_i_arr==NULL)
        printf("An Error Has Occured");
    assert(row_i_arr != NULL);
    row_j_arr = (double *) calloc(a_mat->cols, sizeof(double));
    if(row_j_arr==NULL)
        printf("An Error Has Occured");
    assert(row_j_arr != NULL);
    pivot_arr = (int *) calloc(2, sizeof(int));
    if(pivot_arr==NULL)
        printf("An Error Has Occured");
    assert(pivot_arr != NULL);

    finding_pivot(a_mat, pivot_arr); /*finding relevant i and j*/
    i=pivot_arr[0];
    j=pivot_arr[1];

    /* coping j-row and i-col*/
    for(index = 0; index < a_mat->rows; index++)
    {
        row_i_arr[index] = a_mat->vertices[index][i];
        row_j_arr[index] = a_mat->vertices[index][j];
    }

    /*calculating teta, c,t and s*/
    teta = (a_mat->vertices[j][j] - a_mat->vertices[i][i]) / (2*a_mat->vertices[i][j]);
    if(teta >=0)
    {
        teta_sign = 1;
    }
    else
        teta_sign = -1;

    t = teta_sign / (fabs(teta) + sqrt(teta*teta + 1));
    c = 1 / sqrt(t*t + 1);
    s = t*c;
    /* Calc eigenvalues*/
    for(r = 0; r<a_mat->rows; r++)
    {
        if(r!=i && r!=j)
        {
            a_ri = c*row_i_arr[r] - s*row_j_arr[r];
            a_rj = c*row_j_arr[r] + s*row_i_arr[r];
            a_mat->vertices[r][i] = a_ri;
            a_mat->vertices[i][r] = a_ri;
            a_mat->vertices[r][j] = a_rj;
            a_mat->vertices[j][r] = a_rj;
        }
    }
    a_ii = (c*c*a_mat->vertices[i][i]) + (s*s*a_mat->vertices[j][j]) - (2*s*c*a_mat->vertices[i][j]);
    a_jj = (s*s*a_mat->vertices[i][i]) + (c*c*a_mat->vertices[j][j]) + (2*s*c*a_mat->vertices[i][j]);
    a_mat->vertices[i][i] = a_ii;
    a_mat->vertices[j][j] = a_jj;
    a_mat->vertices[i][j] = 0;
    a_mat->vertices[j][i] = 0;

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
    }
    free(pivot_arr);
    free(row_i_arr);
    free(row_j_arr);
    return;
}

/*This function insert eigenvalues (in a matrix) into an array)*/
void eigenvalues_into_arr(Matrix * mat, double * arr)
{
    int i;
    for (i = 0; i < mat->rows; i++)
    {
        arr[i] = mat->vertices[i][i];
    }
}

/*This function gets an eigen valuematrix and eigenvector matrix, and matches each vector to its value in a given Eigen struct */
void mat_to_eigen_struct(Matrix * a_eigenvalues, Matrix * v_eigenvectors, Eigen** final_eigen)
{
    int i,j;
    for(i = 0; i<a_eigenvalues->rows; i++)
    {
        final_eigen[i] = malloc(sizeof(Eigen));
        if(final_eigen[i]==NULL)
                printf("An Error Has Occured");
        assert(final_eigen[i] != NULL);
        final_eigen[i]->eigen_vector = calloc(a_eigenvalues->rows, sizeof(double));
        if(final_eigen[i]->eigen_vector==NULL)
                printf("An Error Has Occured");
        assert(final_eigen[i]->eigen_vector != NULL);
    }
    for (i = 0; i < a_eigenvalues->rows; i++)
    {
        for (j = 0; j < a_eigenvalues->rows; j++)
        {
            final_eigen[i]->eigen_vector[j] = v_eigenvectors->vertices[j][i];
        }
        final_eigen[i]->eigen_value =  a_eigenvalues->vertices[i][i];
        final_eigen[i]->n = a_eigenvalues->rows;
        final_eigen[i]->origin_index = i;
    }
}

/*This function calculate the eigengap Heuristic that will determine the number of clusters k*/
int eigengap_heuristic(Eigen** eigen, int n)
{
    int i, max_eigengap_index;
    double cur_gap, max_eigengap;
    max_eigengap = 0;
    max_eigengap_index = 0;
    for (i = 0; i < (int)floor(n/2); i++)
    {
        cur_gap = (eigen[i]->eigen_value - eigen[i+1]->eigen_value);
        if(fabs(cur_gap)>max_eigengap)
        {
            max_eigengap = fabs(cur_gap);
            max_eigengap_index = i;
        }
    }
    /*max+1 due to the fact that the first eigen is 0, but in the constructions its 1*/
        return (max_eigengap_index+1);
}

/*This function is the comparator of the qsort function that compares between 2 eigenvalues.
This ia a stable comparator, therefor in case of 2 identical eigenvalues, the comprator will compare the original indices of the vestors*/
int eigen_copm(const void * cp1, const void * cp2)
{
    const Eigen* first = *(const Eigen**) cp1;
    const Eigen* sec = *(const Eigen**) cp2;
    if(first->eigen_value - sec->eigen_value < 0)
    {
        return -1;
    }
    if(first->eigen_value - sec->eigen_value > 0)
    {
        return 1;
    }
    return ((first->origin_index) - (sec->origin_index)); /*identical values, use index*/
}

/*This function gets an Eigen struct and sort it by its eigen values (stable sort).*/
void sort_eigen_values(Eigen** final_eigen)
{
    qsort(final_eigen, final_eigen[0]->n, sizeof(Eigen*), eigen_copm);
    return;
}


/*This funcion gets a matrix and print it s.t the matrix will be separated by a comma, such that each row is in a line of its own*/
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
}

/*This funcion gets a matrix and print its transpose matrix s.t the matrix will be separated by a comma, such that each row is in a line of its own*/
void print_transpose_mat(Matrix * mat)
{
    int i, j;
    for (i = 0; i < mat->rows; i++)
    {
        for (j = 0; j < mat->cols; j++)
        {
            /*converting -0.0000 to 0.0000 - all the values that negative and greaten then -0.00005 should be zero */
            if(mat->vertices[j][i] > -0.00005 && mat->vertices[j][i]<0)
            {
                printf("0.0000");
            }
            else
            {
                printf("%.4f", mat->vertices[j][i]);
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
}

/*This funcion gets an array and its length, and print it s.t each value will be separated by a comma*/
void print_arr(double* arr, int n)
{
    int i;
    for (i = 0; i < n; i++)
    {
        if(arr[i]<0 && arr[i] > -0.00005)
        {
            printf("0.0000");
        }
        else
        {
            printf("%.4f", arr[i]);
        }
        if(i+1!=n)
        {
            printf(",");
        }
    }
}

/*This function gets a matrix, and its initialize it to the identity*/
void init_to_identity(Matrix* mat)
{
    int i;
    mat->vertices = calloc(mat->rows, sizeof(double*));
    for(i=0; i<mat->rows;i++)
    {
        mat->vertices[i] = calloc(mat->rows, sizeof(double));
        if(mat->vertices[i]==NULL)
                printf("An Error Has Occured");
        assert(mat->vertices[i] != NULL);
        mat->vertices[i][i] = 1;
    }
}

/*This function gets a matrix, and its initialize it to an empty matrix*/
void init_mat(Matrix* mat)
{
    int i;
    mat->vertices = calloc(mat->rows, sizeof(double*));
    for(i=0; i<mat->rows;i++)
    {
        mat->vertices[i] = calloc(mat->cols, sizeof(double));
        if(mat->vertices[i]==NULL)
            printf("An Error Has Occured");
        assert(mat->vertices[i] != NULL);
    }
}

/*This function gets to matrixes, and copy the values of the first matrix to the second one*/
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

/*This function gets a matrix and an eigen struct, and copy the vectors in the eigen struc to the matrix as columns*/
void eigen_struct_to_matrix(Matrix * u_mat, Eigen ** final_eigen)
{
    int i,j;
    for (i = 0; i < u_mat->rows; i++)
    {
        for (j = 0; j < u_mat->cols; j++)
        {
            u_mat->vertices[i][j] = final_eigen[j]->eigen_vector[i];
        }
    }
    return;
}

/*This function gets a matrix, and normalizes it (for step 5)*/
Matrix* normalize_matrix(Matrix* mat){
    int i;
    Matrix* normalized_matrix;
    double* normalized_row;

    normalized_matrix = (Matrix*) calloc(1, sizeof(Matrix));
    if(normalized_matrix==NULL)
        printf("An Error Has Occured");
    assert(normalized_matrix!=NULL);
    normalized_matrix->rows = mat->rows;
    normalized_matrix->cols = mat->cols;

    normalized_matrix->vertices = (double **) calloc(mat->rows, sizeof (double *));
    if(normalized_matrix->vertices==NULL)
        printf("An Error Has Occured");
    assert(normalized_matrix->vertices != NULL);
    for (i=0; i<mat->rows; i++){
        normalized_row = normalize_row(mat->vertices[i], mat->cols);
        normalized_matrix->vertices[i] = normalized_row;
    }
    return normalized_matrix;
}

/*This function gets a double array that represents a row, and the num of elements in the row, and normalizes the row*/
double* normalize_row(double* row, int elements_in_row){
    double sum_of_squares, denominator;
    double* normalized_row;
    int i;

    sum_of_squares=0;
    for (i=0; i<elements_in_row; i++){
        sum_of_squares+= pow(row[i], 2);
    }
    if (sum_of_squares==0){
        denominator = 1;
    }
    else{
        denominator = pow(sum_of_squares, 0.5);
    }

    normalized_row = (double *) calloc(elements_in_row, sizeof (double ));
    if(normalized_row==NULL)
        printf("An Error Has Occured");
    assert(normalized_row!=NULL);
    for (i=0; i<elements_in_row; i++){
        normalized_row[i] = row[i]/denominator;
    }
    return normalized_row;
}

/*This function gets 2 double arrays and copy the first ones values to the second matrix (for step 6)*/
void copy_centroid(double* from, double* to, int dim){
    int i;
    for (i=0;i<dim; i++){
        to[i] = from[i];
    }
}

/* This is ta main function that calculates the kmeans algorithm. */
void kmeans(double** observations, int* initial_centroid_indices,
            Cluster** clusters, int k, int n, int dim){

    int i, j, t, count, q, centroid_index, w;
    int differ;
    struct Point** point_arr;

    differ =1;
    i=0;
    point_arr = (struct Point**) calloc(n+1, sizeof(struct Point*));
    if(point_arr==NULL)
        printf("An Error Has Occured");
    assert(point_arr != NULL);

    /* initialize+adapt Point array (as a struct)*/
    for(j=0; j<n; j++)
    {
        point_arr[j] = malloc(sizeof(struct Point));
        if(point_arr[j]==NULL)
            printf("An Error Has Occured");
        assert(point_arr[j] != NULL);
        point_arr[j]->coordinates = calloc(dim+1, sizeof(double));
        if(point_arr[j]->coordinates==NULL)
            printf("An Error Has Occured");
        assert(point_arr[j]->coordinates != NULL);
        /*adapting from python to c*/
        for (w = 0; w< dim; w++)
        {
            point_arr[j]->coordinates[w] = observations[j][w];
        }
        
        point_arr[j]->cluster = -1;
    }


    /* initialize+adapt Clusters array (as a struct)*/
    for(t=0; t<k; t++)
    {
        clusters[t] = malloc(sizeof(Cluster));
        if(clusters[t]==NULL)
            printf("An Error Has Occured");
        assert(clusters[t] != NULL);
        clusters[t]->curr_centroid = (double*) malloc((dim+1)*sizeof(double));
        if(clusters[t]->curr_centroid==NULL)
            printf("An Error Has Occured");
        assert(clusters[t]->curr_centroid != NULL);
        clusters[t]->prev_centroid = (double*) malloc((dim+1)*sizeof(double));
        if(clusters[t]->prev_centroid==NULL)
            printf("An Error Has Occured");
        assert(clusters[t]->prev_centroid != NULL);
        /*adapting from python to c:
        updating the clusters of the observations that were already classified as centroids*/
        point_arr[initial_centroid_indices[t]]->cluster = t;
        centroid_index = initial_centroid_indices[t];
        /*point_arr[centroid_index]->cluster = centroid_index;*/
        clusters[t]->size = 1;
        for(w=0; w<dim; w++)
        {
            clusters[t]->curr_centroid[w] = observations[centroid_index][w];
            clusters[t]->prev_centroid[w] = observations[centroid_index][w];
        }
    }
    for (i=0; i<n; i++){
        /*assign the observations that were not classfied yet to a centroid*/
        if(point_arr[i]->cluster==-1)
        {
            assign_to_closest_cluster(i, clusters, point_arr, 1, dim, k);
        }
    }
    update_centroids(clusters, point_arr, k, dim, n);
    count = 1;

    while((count < MAX_ITER) && differ)
    {
        for(q=0; q<n; q++)
        {
            assign_to_closest_cluster(q, clusters, point_arr, 0, dim, k);
        }

        update_centroids(clusters, point_arr, k, dim, n);
        count++;
        differ = check_difference_in_centroids(clusters, k, dim);
    }
    for (i=0;i<k;i++){
        print_c(clusters[i], dim, i, k);
    }
    for (i=0; i<n; i++){
        free(point_arr[i]->coordinates);
        free(point_arr[i]);
    }
    free(point_arr);
    for (i=0;i<k;i++){
        free(clusters[i]->curr_centroid);
        free(clusters[i]->prev_centroid);
        free(clusters[i]);
    }
}


/*This function gets a 2d cluster array, a point array, and k+dim. It returns a new Cluster 2D object
 where the curr_centroid array is updated  */
void update_centroids(Cluster** clusters, Point** point_arr, int k, int dim, int n)
{
    int i, j, t;
    /*update current centroid to previous one in each cluster*/
    for(i=0; i<k; i++)
    {
        /* each coordinate */
        for(t=0; t<dim; t++)
        {
            clusters[i]->prev_centroid[t] = clusters[i]->curr_centroid[t];
            clusters[i]->curr_centroid[t] = 0.0000L;
        }
    }

    i=0;
    /*for each point*/
    for(j=0; j<n; j++)
    {
        /* Save the point cluster*/
        i = point_arr[j]->cluster;
        /* each coordinate */
        for(t=0; t<dim; t++)
        {
            clusters[i]->curr_centroid[t] += point_arr[j]->coordinates[t];
        }

    }
    /* after we sum the coordinated, we divied each one of them by the sum*/
    for(i=0; i<k; i++)
    {
        if(clusters[i]->size > 1)
        {
            for(t=0; t<dim; t++)
            {
                clusters[i]->curr_centroid[t] = (clusters[i]->curr_centroid[t]) / (clusters[i]->size);
            }
        }
    }
}

/*This function gets a cluster (and k+dim). It returns:
1 - if there is a any difference between the "current" and "prev" of any cluster
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
            if(clusters[i]->curr_centroid[j]!=clusters[i]->prev_centroid[j]){
                return 1; /* we found a difference */
            }
        }
    }

    return 0; /* all of the centroides are the same*/
}

/*This function gets a point number, the points array, and the clusters array, and marches the point the the closest cluster*/
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

        /* need to add the point to relevant cluster + delete it from previous one*/
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

/*This function calculate the distance between a given point and cluster*/
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

/*This function gets a cluster centroid and prints is*/
void print_c(Cluster *cluster, int dim, int cluster_num, int k)
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
    if (cluster_num<(k-1)) {
        printf("\n");
    }
}

/*This function gets a matrix and free its memory*/
void free_matrix(Matrix* matrix){
    int i;

    for (i=0;i<matrix->rows;i++){
        free(matrix->vertices[i]);
    }
    free(matrix->vertices);
    free(matrix);
}


/*This is the main function of the program. Its gets the cmd arguments (k, goal, and the points that in the file path as a Point struct). As weel as it gets the n (num of points), dim and flag (if we came from python or C).
This function allocate and calculates the relevant matrixes, using the kmeans algorithem (if gaol=SPK), print the results and free the allocated memmory */
Matrix* main_logic(int k, char * goal, Point** point_arr, int n, int dim, int flag)
{
    /*  ---- Declaration ----  */
    int  i, count, *initial_centroid_indices;
    Matrix* adj_matrix;
    Matrix* diag_degree_mat;
    Matrix* sqrt_diag_degree_mat;
    Matrix* l_norm_mat;
    Matrix* temp_mat;

    Matrix* a_eigenvalues;
    Matrix* v_eigenvectors;
    Eigen** final_eigen;
    Matrix* u_matrix;
    Matrix* t_matrix;
    double* arr_eigenvalues;
    double prev_off_a, cur_off_a;
    struct Cluster** clusters;

    /* ######## Calculate all relevant matrixes ########*/
    /*---- Weighted Adjacency Matrix (W) ----*/
    adj_matrix = (Matrix*) calloc(n, sizeof(double*));
    if(adj_matrix==NULL)
            printf("An Error Has Occured");
    assert(adj_matrix != NULL);
    adj_matrix->rows = n;
    adj_matrix->cols = n;
    init_mat(adj_matrix);
    to_weighted_adj_mat(point_arr, adj_matrix, dim);

    if(strcmp(goal,"wam")==0)
    {
        print_mat(adj_matrix);
        free_matrix(adj_matrix);
        return NULL;
    }

    /*---- Diagonal Degree Matrix (D) ----*/
    diag_degree_mat = (Matrix*) calloc(n, sizeof(double*));
    if(diag_degree_mat==NULL)
            printf("An Error Has Occured");
    assert(diag_degree_mat != NULL);
    diag_degree_mat->rows = n;
    diag_degree_mat->cols = n;
    init_mat(diag_degree_mat);
    calc_diagonal_degree_mat(diag_degree_mat, adj_matrix);

    if(strcmp(goal,"ddg")==0)
    {
        print_mat(diag_degree_mat);
        free_matrix(adj_matrix);
        free_matrix(diag_degree_mat);
        return NULL;
    }

    /*----  Making D into D^(-1/2) ----*/
    sqrt_diag_degree_mat = (Matrix*) calloc(n, sizeof(double*));
    if(sqrt_diag_degree_mat==NULL)
            printf("An Error Has Occured");
    assert(sqrt_diag_degree_mat != NULL);
    sqrt_diag_degree_mat->rows = n;
    sqrt_diag_degree_mat->cols = n;
    init_mat(sqrt_diag_degree_mat);
    copy_mat_a_to_b(diag_degree_mat, sqrt_diag_degree_mat); /*coping D */
    to_sqrt_diag_mat(sqrt_diag_degree_mat); /*changing the matrix to D^(-1/2)*/

    /*---- L Norm Matrix (lnorm)----*/
    l_norm_mat = (Matrix*) calloc(n, sizeof(double*));
    if(l_norm_mat==NULL)
            printf("An Error Has Occured");
    assert(l_norm_mat != NULL);
    l_norm_mat->rows = n;
    l_norm_mat->cols = n;
    init_mat(l_norm_mat);

    /* init additional temp matrix (for calculation purposes) */
    temp_mat = (Matrix*) calloc(n, sizeof(double*));
    if(temp_mat==NULL)
            printf("An Error Has Occured");
    assert(temp_mat != NULL);
    temp_mat->rows = n;
    temp_mat->cols = n;
    init_mat(temp_mat);

    /*---- Calc lnorm matrix ----*/
    mult_diag_matrix(sqrt_diag_degree_mat, adj_matrix, 1, temp_mat); /*   D^(-0.5) x (W)  */
    mult_diag_matrix(temp_mat, sqrt_diag_degree_mat, 0, l_norm_mat); /*  (D^(-0.5) * W)  x  (D^(-0.5)) */
    to_l_norm(l_norm_mat);     /*   I - (A*B)*A        */

    if(strcmp(goal,"lnorm")==0)
    {
        print_mat(l_norm_mat);
        free_matrix(adj_matrix);
        free_matrix(diag_degree_mat);
        free_matrix(sqrt_diag_degree_mat);
        free_matrix(temp_mat);
        free_matrix(l_norm_mat);
        return NULL;
    }

    /*---- Calc A matrix (eiganvalues) and V matrix (eiganvectors) ----*/
    /* init eigenvalues */
    a_eigenvalues = (Matrix*) calloc(n, sizeof(double*));
    if(a_eigenvalues==NULL)
            printf("An Error Has Occured");
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
    if(v_eigenvectors==NULL)
            printf("An Error Has Occured");
    assert(v_eigenvectors != NULL);
    v_eigenvectors->rows = n;
    v_eigenvectors->cols = n;
    init_to_identity(v_eigenvectors);

    /*----  Calc the eigenvectors and eigenvalues  ----*/
    count = 0;
    cur_off_a = off_diag_squares_sum(a_eigenvalues);
    do{
        count++;
        converting_a_and_v_mat(a_eigenvalues, v_eigenvectors);
        prev_off_a = cur_off_a;
        cur_off_a=off_diag_squares_sum(a_eigenvalues);
        /*printf("prev is: %.16f\n", prev_off_a);
        printf("cur is: %.16f\n", cur_off_a);
        printf("prev - cur is: %.16f\n", prev_off_a-cur_off_a);*/
    } while (count<MAX_LOOPS && (prev_off_a-cur_off_a) > EPSILON);

    arr_eigenvalues = (double*) calloc(n, sizeof(double));
    if(arr_eigenvalues==NULL)
            printf("An Error Has Occured");
    assert(arr_eigenvalues != NULL);
    eigenvalues_into_arr(a_eigenvalues, arr_eigenvalues);

    if(strcmp(goal, "spk")==0)
    {
        /* creating Eigen struct and coping the values of the relevant matrices to it*/
        final_eigen = (Eigen**) calloc(n, sizeof(Eigen*));
        if(final_eigen==NULL)
            printf("An Error Has Occured");
        assert(final_eigen != NULL);
        mat_to_eigen_struct(a_eigenvalues, v_eigenvectors, final_eigen);
        sort_eigen_values(final_eigen); /*sorting the eigenvalues for creating U anf for k_heuristic*/
        if(k==0)
        {
            k = eigengap_heuristic(final_eigen, n);
        }

        u_matrix = (Matrix*) calloc(n, sizeof(double*));
        if(u_matrix==NULL)
            printf("An Error Has Occured");
        assert(u_matrix != NULL);
        u_matrix->rows = n;
        u_matrix->cols = k;

        init_mat(u_matrix);
        eigen_struct_to_matrix(u_matrix, final_eigen);
        t_matrix = normalize_matrix(u_matrix);

        if(flag) /* 1 is python, 0 is C*/
        {
            free_matrix(adj_matrix);
            free_matrix(diag_degree_mat);
            free_matrix(sqrt_diag_degree_mat);
            free_matrix(temp_mat);
            free_matrix(l_norm_mat);
            free_matrix(a_eigenvalues);
            free_matrix(v_eigenvectors);
            free(arr_eigenvalues);
            for(i=0;i<n;i++){
                free(final_eigen[i]->eigen_vector);
                free(final_eigen[i]);
            }
            free(final_eigen);
            free_matrix(u_matrix);

            return t_matrix;
        }
        else
        {
            initial_centroid_indices = (int*) calloc(k+1, sizeof (int));
            for (i=0;i<k;i++){
                initial_centroid_indices[i] = i;
            }
            clusters = (struct Cluster**) calloc(k+1, sizeof(struct Cluster*));
            if(clusters==NULL)
                printf("An Error Has Occured");
            assert(clusters != NULL);
            kmeans(t_matrix->vertices, initial_centroid_indices, clusters, k, n, k);

            free(initial_centroid_indices);
            free(clusters);
            for(i=0;i<n;i++){
                free(final_eigen[i]->eigen_vector);
                free(final_eigen[i]);
            }
            free(final_eigen);
            free_matrix(u_matrix);
            free_matrix(t_matrix);

        }
    }

    if(strcmp(goal, "jacobi")==0)
    {
        print_arr(arr_eigenvalues, n); /* print the eigen values*/
        printf("\n");
        print_transpose_mat(v_eigenvectors);
    }

    free_matrix(adj_matrix);
    free_matrix(diag_degree_mat);
    free_matrix(sqrt_diag_degree_mat);
    free_matrix(temp_mat);
    free_matrix(l_norm_mat);
    free_matrix(a_eigenvalues);
    free_matrix(v_eigenvectors);
    free(arr_eigenvalues);
    return NULL;
}

int main(int argc, char** argv)
{
    /*CMD: prog_name, k, goal, file.mane*/
    FILE *data;
    int k, n, dim, i;
    char *goal, *filename;
    Point **point_arr;

    /*Input validation*/
    if (argc != 4 || (!is_int(argv[1])))
    {
        printf("Invalid Input!");
        exit(1);
    }

    k = atoi(argv[1]);
    goal = argv[2];
    filename = argv[3];

    data = fopen(filename, "r");
    if(data==NULL){
        printf("Invalid Input!");
        exit(1);
    }

    /* ---- Calculation of input dimentions (n, dim) ----*/
    n = get_num_points(data);
    dim = get_dim(data);
    if (n == 0 || dim == 0 || k>=n){
        printf("Invalid Input!");
        fclose(data);
        exit(1);
    }
    /*  ---- Input to struct ----  */
    point_arr = (Point**) calloc(n+1, sizeof(Point*));
    if(point_arr==NULL)
                printf("An Error Has Occured");
    assert(point_arr != NULL);
    input_to_points_struct(data, point_arr, dim, n);

    if (strcmp(goal, "wam") != 0 && strcmp(goal, "ddg")!=0 && strcmp(goal, "lnorm")!=0 &&
            strcmp(goal, "jacobi")!=0 && strcmp(goal, "spk")!=0){
        printf("Invalid Input!\n");
        exit(1);
    }
    main_logic(k, goal, point_arr, n, dim, 0);

    for(i=0;i<n;i++){
        free(point_arr[i]->coordinates);
        free(point_arr[i]);
    }
    free(point_arr);

    fclose(data);
    return 1;
}