#include "gsl_matrix_methods.hpp"

using namespace std;

//@author bjd2385
//https://gist.github.com/bjd2385/7f4685e703f7437e513608f41c65bbd7
gsl_matrix *invert_a_matrix(gsl_matrix *matrix) {

    int size = matrix->size1; //assume equal rows and cols
    gsl_matrix* matrix_copy = gsl_matrix_calloc(size, size);

    for (int i = 0; i < size; i++) {
	for (int j = 0; j < size; j++) {
	    gsl_matrix_set(matrix_copy, i, j, gsl_matrix_get(matrix, i, j));
	}
    }
    gsl_permutation *p = gsl_permutation_alloc(size);
    int s;

    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(matrix_copy, p, &s);

    // Compute the  inverse of the LU decomposition
    gsl_matrix *inv = gsl_matrix_alloc(size, size);
    gsl_linalg_LU_invert(matrix_copy, p, inv);

    gsl_permutation_free(p);

    return inv;
}

void print_mat_contents(gsl_matrix *matrix)
{
    size_t i, j;
    double element;
    int size1 = matrix->size1;
    int size2 = matrix->size2;

    for (i = 0; i < size1; ++i) {
        for (j = 0; j < size2; ++j) {
            element = gsl_matrix_get(matrix, i, j);
            printf("%f ", element);
        }
        printf("\n");
    }
}


//test if matrix inversion is correct
void gsl_matrix_inversion_test()
{
    gsl_matrix *m = gsl_matrix_calloc(4,4);
    double tdarr[4][4] = {{1,2,3,4},{1,3,5,6},{9,3,4,9},{5,6,1,1}};
    for (int i = 0; i < 4; i++) {
	for (int j = 0; j < 4; j++) {
	    gsl_matrix_set(m, i, j, tdarr[i][j]);
	}
    }
    gsl_matrix *minv = invert_a_matrix(m);
    print_mat_contents(m);
    std::cout<<"================inv================="<<endl;
    print_mat_contents(minv);
}

//copy contents of vertex.Coord to given gsl matrix
void copy_vertex_to_gsl_matrix(vector<Vertex>& vertex, gsl_matrix *matrix) {
    for (int i = 0; i < vertex.size(); i++) {
        vector<double> * coord = &(vertex[i].Coord);
        for (int j = 0; j < 3; j++) {
            gsl_matrix_set(matrix, i, j, (*coord)[j]);
        }
    }
}

/* @deprecated warning: NOT WOKRING AS INTENDED  
 *
 * This method serves to solve for the best fit / LS fit plane for
 * the giving (n, 3) vector in cartesian coordinate system.
 * The solution is by solving (LS-solution) the equation:
 * 
 * Av = B
 * 
 * where A = [X, Y, 1], B = [Z], and solution v is giving by:
 * 
 * v = [a, b, c] where the plane is giving by ax + by + c = z 
 * 
 * The LSS can be calculated with Moore-Penrose inverse:
 * 
 * v = (A^T A)^-1 A^T B
 * 
 * Input: (n, 3) vector for spline points
 * Output: vector<double>(3) as [a, b, c]
 */
vector<double> leastSquarePlane(vector<vector<double>>& splinePoints) {

    //get length of spline points
    int length = splinePoints.size();

    //alloc matrix A and B
    gsl_matrix *matA = gsl_matrix_alloc(length, 3);
    gsl_matrix *matB = gsl_matrix_alloc(length, 1);

    //iterate over length and assign value to A and B
    for (int i = 0; i < length; i++) {
        //set matA to [X, Y, 1]
        gsl_matrix_set(matA, i, 0, splinePoints[i][0]);
        gsl_matrix_set(matA, i, 1, splinePoints[i][1]);
        gsl_matrix_set(matA, i, 2, 1.0);
        //set matB to [Z]
        gsl_matrix_set(matB, i, 0, splinePoints[i][2]);
    }

    //calculate matA transpose (matAT)
    gsl_matrix *matAT = gsl_matrix_alloc(3, length);
    gsl_matrix_transpose_memcpy(matAT, matA);

    //calculate (A^T * A)^(-1) -> (matATAinv)
    gsl_matrix *matATA = gsl_matrix_alloc(3, 3);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, matAT, matA,
                        0.0, matATA);
    gsl_matrix *matATAinv = invert_a_matrix(matATA);

    //calculate ATAinv * AT * B -> [a, b, c]
    gsl_matrix *matATB = gsl_matrix_alloc(3, 1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, matAT, matB,
                        0.0, matATB);
    gsl_matrix *planeParams = gsl_matrix_alloc(3, 1);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, matATAinv, matATB,
                        0.0, planeParams);

    //extract value to vector<double>
    vector<double> planeParamsVec(3, 0.0);
    for (int i = 0; i < 3; i++) {
        planeParamsVec[i] = gsl_matrix_get(planeParams, i, 0);
    }

    return planeParamsVec;    
}