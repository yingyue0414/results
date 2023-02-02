#pragma once

#include <math.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>
#include <omp.h>
#include <iomanip>
#include <random>
#include <ctime>       /* clock_t, clock, CLOCKS_PER_SEC */
#include "mesh.hpp"
#include "vector_math.hpp" 
#include "parameters.hpp"
#include "gsl_matrix_methods.hpp"
#include <gsl/gsl_matrix_double.h>

//mesh2surface / surface2mesh
void assignMesh2Surface(gsl_matrix *mesh2surface, std::vector<Vertex>& vertex, std::vector<Face>& face);

//boundary condition
int get_relatve_pt_periodic(int i, int n, int m);
void postprocess_ghost_periodic(Param& param, gsl_matrix *verticesOnMesh, std::vector<Vertex>& vertex);
void postprocess_ghost_free (Param& param, gsl_matrix *verticesOnMesh, std::vector<Vertex>& vertex, std::vector<Face>& face);
void next_step(std::vector<Vertex>& vertex, gsl_matrix *verticesProjSurface, double forceScaleConst, double randScaleConst, std::mt19937& gen, std::normal_distribution<>& normal_dist, bool isBoundaryPeriodic);