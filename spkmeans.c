#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#define EPSILON pow(10,-15)

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
int eigengap_heuristic(int* eigenvalues, int n);
void converting_a_and_v_mat(Matrix * a_mat, Matrix * v_mat);
void init_to_identity(Matrix* mat);
void init_mat(Matrix* mat);
void copy_mat_a_to_b(Matrix* mat_a, Matrix* mat_b);
void to_l_norm(Matrix *mat);

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
            diag_mat->vertices[i][i] = row_sum;
        }
    }
}

void to_sqrt_diag_mat(Matrix* diag_mat)
{   
    int i;
    for(i=0; i<diag_mat->rows;i++)
    {
        diag_mat->vertices[i][i] = 1/sqrt(diag_mat->vertices[i][i]);
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
    teta = (a_mat->vertices[j][j] - a_mat->vertices[i][i]) / (2*max_val);
    
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

    free(row_i_arr);
    free(row_j_arr);

    a_ij = v_mat->rows;
    return;
    
}

int eigengap_heuristic(int* eigenvaluse, int n)
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

void print_mat(Matrix *mat)
{
    int i, j, count;
    count = 0;
    for (i = 0; i < mat->rows; i++)
    {
        for (j = 0; j < mat->rows; j++)
        {
            printf("%f", mat->vertices[i][j]);
            if(j!=mat->rows-1)
            {
                printf(", ");
            }
            if(mat->vertices[i][j]>0.00001 || mat->vertices[i][j]<-0.0001)
            {
                count++;
            }
        }
         printf("\n \n");
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
        mat->vertices[i] = calloc(mat->rows, sizeof(double));
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
    double prev_off_a, cur_off_a;

    /*  ---- Validation ----  */
    filename = fopen(f_name, "r");
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
    if(k==0)
    {
        /*use the heuristic 1.3*/
    }

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

    sqrt_diag_degree_mat = (Matrix*) calloc(n, sizeof(double*));
    assert(sqrt_diag_degree_mat != NULL);
    sqrt_diag_degree_mat->rows = n;    

    init_mat(sqrt_diag_degree_mat);
    copy_mat_a_to_b(diag_degree_mat, sqrt_diag_degree_mat);
    to_sqrt_diag_mat(sqrt_diag_degree_mat);

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
    /*   D^(-0.5) x (W)  */
    mult_diag_matrix(sqrt_diag_degree_mat, adj_matrix, 1, temp_mat);

    /*  (D^(-0.5) * W)  x  (D^(-0.5)) */
    mult_diag_matrix(temp_mat, sqrt_diag_degree_mat, 0, l_norm_mat);
    /*   I - (A*B)*A        */
    to_l_norm(l_norm_mat);


    /*  Calc a (eiganvalues) and v (eiganvectors) */
    a_eigenvalues = (Matrix*) calloc(n, sizeof(double*));
    assert(a_eigenvalues != NULL);
    a_eigenvalues->rows = n;
    a_eigenvalues->cols = n;
    init_mat(a_eigenvalues);
    copy_mat_a_to_b(l_norm_mat, a_eigenvalues);

    v_eigenvectors = (Matrix*) calloc(n, sizeof(double*));
    assert(v_eigenvectors != NULL);
    v_eigenvectors->rows = n;
    v_eigenvectors->cols = n;
    init_to_identity(v_eigenvectors);

    prev_off_a = 0;
    cur_off_a = off_diag_squares_sum(l_norm_mat);
    count = 0;
    while (count<50 && (prev_off_a-cur_off_a) <= EPSILON)
    {
        prev_off_a = cur_off_a;
        count++;
        converting_a_and_v_mat(a_eigenvalues, v_eigenvectors);
        cur_off_a = off_diag_squares_sum(a_eigenvalues);
        /*
        printf("At loop %d, eigenvalues mat is:\n", count);
        print_mat(a_eigenvalues);
        */
    }

    printf("Eigenvalues mat is:\n");
    print_mat(a_eigenvalues);
    printf("Was in while %d times\n",count);
    

    if(strcmp(goal,"spk")==0)
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

    else if(strcmp(goal,"wam")==0)
    {   
        printf("The weight adj mat is:\n");
        print_mat(adj_matrix);
    }
    else if(strcmp(goal,"ddg")==0)
    {
        printf("diag:\n");
        print_mat(diag_degree_mat);
        printf("diag^(-0.5):\n");
        print_mat(sqrt_diag_degree_mat);
    }
    else if(strcmp(goal,"lnorm")==0)
    {
        printf("l_norm:\n");
        print_mat(l_norm_mat);
    }
    else if(strcmp(goal,"jacobi")==0)
    {
        printf("In construction\n");
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
        free(a_eigenvalues->vertices[i]);
        free(v_eigenvectors->vertices[i]);
    }
    free(point_arr);
    free(adj_matrix);
    free(diag_degree_mat);
    free(sqrt_diag_degree_mat);
    free(l_norm_mat);
    free(temp_mat);
    free(a_eigenvalues);
    free(v_eigenvectors);
    fclose(filename);

    return 1;
}

int main(int argc, char** argv)
{
    /*CMD: prog_name, k, goal, file.mane*/

    int k;
    char* goal;
    char* f_name;
   
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

