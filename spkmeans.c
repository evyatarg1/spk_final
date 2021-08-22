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
    int size;
    double** vertexs;
} Matrix;


int get_num_points(FILE *filename);
int get_dim(FILE *filename);
int is_int(char* p);
void input_to_points_struct(FILE *filename, Point** point_arr, int dim);
void to_weighted_adj_mat(Point** point_arr, Matrix* adj_matrix, int n, int dim);
double calc_weight(Point *first, Point *sec, int dim);
void calc_diagonal_degree_mat(Matrix* diag_mat, Matrix * adj_matrix, int n);
double calc_row_sum(Matrix * adj_matrix, int row, int n);
void mult_diag_matrix_first_diag(Matrix * first, Matrix * sec, Matrix * res, int n);
void mult_diag_matrix_sec_diag(Matrix * first_diag, Matrix * sec, Matrix * res, int n);
void print_mat(Matrix *mat, int n);
void to_l_norm(Matrix *mat, int n);

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

void to_weighted_adj_mat(Point** point_arr, Matrix* adj_matrix, int n, int dim)
{
    int i,j;
    adj_matrix->vertexs = calloc(n, sizeof(double*));
    
    for(i=0; i<n;i++)
    {
        adj_matrix->vertexs[i] = calloc(n, sizeof(double));
        assert(adj_matrix->vertexs[i] != NULL);
    }
    
    for(i=0; i<n;i++)
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
                adj_matrix->vertexs[i][j] = (double)calc_weight(point_arr[i], point_arr[j], dim);
                adj_matrix->vertexs[j][i] = adj_matrix->vertexs[i][j];
            }
        }
        adj_matrix->vertexs[i][i] = 0;
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

void calc_diagonal_degree_mat(Matrix* diag_mat, Matrix * adj_matrix, int n)
{
    /*calc D^(-0.5)*/
    int i;
    double row_sum;
    diag_mat->vertexs = calloc(n, sizeof(double*));
    for(i=0; i<n;i++)
    {
        diag_mat->vertexs[i] = calloc(n, sizeof(double));
        assert(diag_mat->vertexs[i] != NULL);
        row_sum = calc_row_sum(adj_matrix, i, n);
        if(row_sum>0)
        {
            diag_mat->vertexs[i][i] = 1/sqrt(row_sum);
        }
    }
}

double calc_row_sum(Matrix * adj_matrix, int row, int n)
{
    int i;
    double res;
    for(i=0; i<n;i++)
    {
        res+= adj_matrix->vertexs[row][i];
    }
    return res;
}


void mult_diag_matrix_first_diag(Matrix * first_diag, Matrix * sec, Matrix * res, int n)
{
    int i, j;
    
    /*The first mat must be diagonal! */
    
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            res->vertexs[i][j] += (first_diag->vertexs[i][i])*(sec->vertexs[i][j]); 
        }
    }
    return;
}

void mult_diag_matrix_sec_diag(Matrix * first_diag, Matrix * sec, Matrix * res, int n)
{
    int i, j;
    
    /*The sec mat must be diagonal! */
    
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            res->vertexs[i][j] += (first_diag->vertexs[i][j])*(sec->vertexs[j][j]); 
        }
    }
    return;
}

void to_l_norm(Matrix *mat, int n)
{
    int i,j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if(i==j)
            {
                mat->vertexs[i][j] = 1-mat->vertexs[i][j];

            }
            else
            {
                mat->vertexs[i][j] = (-1)*(mat->vertexs[i][j]);
            }
        }
    }
    return;  
}

void print_mat(Matrix *mat, int n)
{
    int i, j, count;
    count = 0;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("%f, ", mat->vertexs[i][j]);
            if(mat->vertexs[i][j]>0.00001 || mat->vertexs[i][j]<-0.0001)
            {
                count++;
            }
        }
         printf("\n \n");
    }
    printf("there are %d cells that are not 0\n", count);
    
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
    printf("bla\n");
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
    adj_matrix->size = n;
    to_weighted_adj_mat(point_arr, adj_matrix, n, dim);

    /*---- Diagonal Degree Matrix ----*/
    diag_degree_mat = (Matrix*) calloc(n, sizeof(double*));
    assert(diag_degree_mat != NULL);
    diag_degree_mat->size = n;
    calc_diagonal_degree_mat(diag_degree_mat, adj_matrix, n);    

    /*---- L Norm Matrix ----*/
    l_norm_mat = (Matrix*) calloc(n, sizeof(double*));
    assert(l_norm_mat != NULL);
    l_norm_mat->size = n;
    l_norm_mat->vertexs = calloc(n, sizeof(double*));
    for(i=0; i<n;i++)
    {
        l_norm_mat->vertexs[i] = calloc(n, sizeof(double));
        assert(l_norm_mat->vertexs[i] != NULL);
    }

    temp_mat = (Matrix*) calloc(n, sizeof(double*));
    assert(temp_mat != NULL);
    temp_mat->size = n;
    temp_mat->vertexs = calloc(n, sizeof(double*));
    for(i=0; i<n;i++)
    {
        temp_mat->vertexs[i] = calloc(n, sizeof(double));
        assert(temp_mat->vertexs[i] != NULL);
    }

    printf("The weight adj mat is:\n");
    print_mat(adj_matrix, n);
    
    printf("A = Diag:\n");
    print_mat(diag_degree_mat, n);

    printf("B = adj_matrix:\n");
    print_mat(adj_matrix, n);

    /*   D^(-0.5) x (W)  */
    mult_diag_matrix_first_diag(diag_degree_mat, adj_matrix, temp_mat, n); 

    printf("A*B \n");
    print_mat(temp_mat, n);

    /*  (D^(-0.5) * W)  x  (D^(-0.5)) */
    mult_diag_matrix_sec_diag(temp_mat, diag_degree_mat, l_norm_mat, n);
    printf("(A*B)*A:\n");
    print_mat(l_norm_mat, n);

    to_l_norm(l_norm_mat, n);
    printf("I - (A*B)*A:\n");

    print_mat(l_norm_mat, n);


    printf("in main3\n");

    /* ######### FREE ######### */
    for(i=0; i<n; i++)
    {
        free(point_arr[i]->coordinates);
        free(point_arr[i]);

        free(adj_matrix->vertexs[i]);
        free(diag_degree_mat->vertexs[i]);
        free(l_norm_mat->vertexs[i]);
        free(temp_mat->vertexs[i]);
    }
    free(point_arr);
    free(adj_matrix);
    free(diag_degree_mat);
    free(l_norm_mat);
    free(temp_mat);
    fclose(filename);
    return 1;
}

