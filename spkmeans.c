#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

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
    int columns;
    double** vertices;
} Matrix;


int get_num_points(FILE *filename);
int get_dim(FILE *filename);
int is_int(char* p);
void input_to_points_struct(FILE *filename, Point** point_arr, int dim);
void to_weighted_adj_mat(Point** point_arr, Matrix* adj_matrix, int dim);
double calc_weight(Point *first, Point *sec, int dim);
void calc_diagonal_degree_mat(Matrix* diag_mat, Matrix * adj_matrix);
double calc_row_sum(Matrix * adj_matrix, int row);
void mult_diag_matrix(Matrix * first, Matrix * sec, int is_first_diag, Matrix * res);
void print_mat(Matrix *mat);
void to_l_norm(Matrix *mat);
Matrix* normalize_matrix(Matrix* mat);
double* normalize_row(double* row, int elements_in_row);


void kmeans_logic(Cluster** clusters, int k, int n, int dim, int max_iter, Matrix* mat);
int assign_to_closest_cluster(int point_num, Cluster** clusters, Point** points, int first_insert, int dim, int k);
double find_distance(Point* point, Cluster* cluster, int dim);
void print_c(Cluster *cluster, int dim);
int update_centroids(Cluster** clusters, Point** point_arr, int k, int dim, int n);
int check_difference_in_centroids(Cluster** clusters, int k, int dim);
void print_matrix(Matrix *mat);


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
    printf("int get num3\n");
    rewind(filename);
    printf("Num of points %d\n", num_of_points);

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
    printf("Dim is %d\n", dim);

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

void to_weighted_adj_mat(Point** point_arr, Matrix* adj_matrix, int dim)
{
    int i,j;
    adj_matrix->vertices = calloc(adj_matrix->rows, sizeof(double*));
    
    for(i=0; i<adj_matrix->rows;i++)
    {
        adj_matrix->vertices[i] = calloc(adj_matrix->rows, sizeof(double));
        assert(adj_matrix->vertices[i] != NULL);
    }
    
    for(i=0; i<adj_matrix->rows;i++)
    {
        for(j=0;j<=i;j++)
        {
            if(i!=j)
            {
                /*if(j==0 && i==1)
                {
                    printf("point 1: %f %f %f\n", point_arr[i]->coordinates[0], point_arr[i]->coordinates[1], point_arr[i]->coordinates[2]);
                    printf("point 2: %f %f %f\n", point_arr[j]->coordinates[0], point_arr[j]->coordinates[1], point_arr[j]->coordinates[2]);

                }*/
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
    /*calc D^(-0.5)*/
    int i;
    double row_sum;
    diag_mat->vertices = calloc(diag_mat->rows, sizeof(double*));
    for(i=0; i<diag_mat->rows;i++)
    {
        diag_mat->vertices[i] = calloc(diag_mat->rows, sizeof(double));
        assert(diag_mat->vertices[i] != NULL);
        row_sum = calc_row_sum(adj_matrix, i);
        if(row_sum>0)
        {
            diag_mat->vertices[i][i] = 1 / sqrt(row_sum);
        }
    }
}

double calc_row_sum(Matrix * adj_matrix, int row)
{
    int i;
    double res;
    for(i=0; i<adj_matrix->rows;i++)
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
                res->vertices[i][j] += (first->vertices[i][i]) * (sec->vertices[i][j]);
            }
            else /*sec mat must be diagonal*/
            {
                res->vertices[i][j] += (first->vertices[i][j]) * (sec->vertices[j][j]);

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
                mat->vertices[i][j] = 1 - mat->vertices[i][j];

            }
            else
            {
                mat->vertices[i][j] = (-1) * (mat->vertices[i][j]);
            }
        }
    }
    return;  
}

void print_mat(Matrix *mat)
{
    int i, j, count;
    double a;
    count = 0;
    for (i = 0; i < mat->rows; i++)
    {
        for (j = 0; j < mat->columns; j++)
        {
            a = mat->vertices[i][j];
            printf("%f, ", mat->vertices[i][j]);
            if(mat->vertices[i][j] > 0.00001 || mat->vertices[i][j] < -0.0001)
            {
                count++;
            }
        }
         printf("\n \n");
    }
    printf("there are %d cells that are not 0\n", count);
    
}

/**step 5 functions**/

Matrix* normalize_matrix(Matrix* mat){
    int i;
    Matrix* normalized_matrix;
    double* normalized_row;

    normalized_matrix = (Matrix*) calloc(1, sizeof(Matrix));
    assert(normalized_row!=NULL);
    normalized_matrix->rows = mat->rows;
    normalized_matrix->columns = mat->columns;
    normalized_matrix->vertices = (double **) calloc(mat->rows, sizeof (double *));
    assert(normalized_matrix->vertices != NULL);
    for (i=0; i<mat->rows; i++){
        normalized_row = normalize_row(mat->vertices[i], mat->columns);
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
    int cnt, i, count, q;
    int differ;
    double value;
    char ch;
    Point** point_arr;
    differ =1;
    cnt = 0;
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
    dim = mat->columns;
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

int tester_main(int argc, char** argv){
    Matrix* mat, *normalized_matrix;
    int i;

    mat = (Matrix*) calloc(1, sizeof (Matrix));
    mat->rows=6;
    mat->columns=3;
    mat->vertices = (double **) calloc(mat->rows, sizeof (double *));
    for (i=0; i<mat->rows; i++){
        mat->vertices[i] = (double *) calloc(mat->columns, sizeof (double ));
    }
    mat->vertices[0][0] = 8.1402;
    mat->vertices[0][1] =-5.8022;
    mat->vertices[0][2] =-7.2376;
    mat->vertices[1][0] =10.1626;
    mat->vertices[1][1] =-7.4824;
    mat->vertices[1][2] =-6.5774;
    mat->vertices[2][0] =9.3153;
    mat->vertices[2][1] =-5.4974;
    mat->vertices[2][2] =-6.7025;
    mat->vertices[3][0] =9.7950;
    mat->vertices[3][1] =-5.2550;
    mat->vertices[3][2] =-9.1757;
    mat->vertices[4][0] =7.9095;
    mat->vertices[4][1] =-6.6488;
    mat->vertices[4][2] =-7.6088;
    mat->vertices[5][0] =10.3309;
    mat->vertices[5][1] =-5.3494;
    mat->vertices[5][2] =-5.9993;

    kmeans(mat, 3);
    normalized_matrix = normalize_matrix(mat);
    print_matrix(normalized_matrix);
    return 1;
}

int main(int argc, char** argv)
{
    /*CMD: prog_name, k, goal, file.mane*/
    /*  ---- Declaration ----  */
    FILE *filename;
    int k, n, dim, i;
    char* goal;
    Matrix* adj_matrix;
    Matrix* diag_degree_mat;
    Matrix* l_norm_mat;
    Matrix* temp_mat;
    Point** point_arr;    
    
    /*  ---- Validation ----  */
    if (argc != 4 || (!is_int(argv[1]))) 
    {
        printf("Invalid Input!\n");
        exit(1);
    }

    k = atoi(argv[1]);
    goal = argv[2];
    filename = fopen(argv[3], "r");
    if(filename==NULL)
    {
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


    /* ---- Checking the input goal ----*/
    printf("goal is %s\n", goal);
    if(k==0)
    {
        /*use the heuristic 1.3*/
    }
    /*if(goal == "spk")
    {
        printf("In construction\n");
    }
    else if(goal =="wam")
    {   
        printf("In construction\n");
    }
    else if(goal =="ddg")
    {
        printf("In construction\n");
    }
    else if(goal =="lnorm")
    {
        printf("In construction\n");
    }
    else if(goal =="jacobi")
    {
        printf("In construction\n");
    }
    else
    {
        printf("Invalid Input!\n");
        exit(1);  
    }*/


    /* ######## Calculate all relevant matrixes ########*/

    /*---- Weighted Adjacency Matrix----*/
    adj_matrix = (Matrix*) calloc(n, sizeof(double*));
    assert(adj_matrix != NULL);
    adj_matrix->rows = n;
    to_weighted_adj_mat(point_arr, adj_matrix, dim);

    /*---- Diagonal Degree Matrix ----*/
    diag_degree_mat = (Matrix*) calloc(n, sizeof(double*));
    assert(diag_degree_mat != NULL);
    diag_degree_mat->rows = n;
    calc_diagonal_degree_mat(diag_degree_mat, adj_matrix);    

    /*---- L Norm Matrix ----*/
    l_norm_mat = (Matrix*) calloc(n, sizeof(double*));
    assert(l_norm_mat != NULL);
    l_norm_mat->rows = n;
    l_norm_mat->vertices = calloc(n, sizeof(double*));
    for(i=0; i<n;i++)
    {
        l_norm_mat->vertices[i] = calloc(n, sizeof(double));
        assert(l_norm_mat->vertices[i] != NULL);
    }

    temp_mat = (Matrix*) calloc(n, sizeof(double*));
    assert(temp_mat != NULL);
    temp_mat->rows = n;
    temp_mat->vertices = calloc(n, sizeof(double*));
    for(i=0; i<n;i++)
    {
        temp_mat->vertices[i] = calloc(n, sizeof(double));
        assert(temp_mat->vertices[i] != NULL);
    }

    printf("The weight adj mat is:\n");
    print_mat(adj_matrix);
    
    printf("A = Diag:\n");
    print_mat(diag_degree_mat);

    printf("B = adj_matrix:\n");
    print_mat(adj_matrix);

    /*   D^(-0.5) x (W)  */
    mult_diag_matrix(diag_degree_mat, adj_matrix, 1, temp_mat); 

    printf("A*B \n");
    print_mat(temp_mat);

    /*  (D^(-0.5) * W)  x  (D^(-0.5)) */
    mult_diag_matrix(temp_mat, diag_degree_mat, 0, l_norm_mat);
    printf("(A*B)*A:\n");
    print_mat(l_norm_mat);

    to_l_norm(l_norm_mat);
    printf("I - (A*B)*A:\n");

    print_mat(l_norm_mat);


    printf("in main3\n");

    /* ######### FREE ######### */
    for(i=0; i<n; i++)
    {
        free(point_arr[i]->coordinates);
        free(point_arr[i]);

        free(adj_matrix->vertices[i]);
        free(diag_degree_mat->vertices[i]);
        free(l_norm_mat->vertices[i]);
        free(temp_mat->vertices[i]);
    }
    free(point_arr);
    free(adj_matrix);
    free(diag_degree_mat);
    free(l_norm_mat);
    free(temp_mat);
    fclose(filename);
    return 1;
}

