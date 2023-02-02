#pragma once

#include <math.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <omp.h>
#include <iomanip>
#include <ctime>       /* clock_t, clock, CLOCKS_PER_SEC */
#include "mesh.hpp"
#include "vector_math.hpp" 
#include "parameters.hpp"

#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

gsl_matrix *invert_a_matrix(gsl_matrix *matrix);
void print_mat_contents(gsl_matrix *matrix);
void gsl_matrix_inversion_test();
void copy_vertex_to_gsl_matrix(std::vector<Vertex>& vertex, gsl_matrix *matrix);
std::vector<double> leastSquarePlane(std::vector<std::vector<double>>& splinePoints);