#pragma once

#include <math.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <ctime>       /* clock_t, clock, CLOCKS_PER_SEC */
#include "single_particle_diffusion.hpp"
#include "vector_math.hpp"
#include "parameters.hpp"
#pragma omp declare reduction(+: Energy : omp_out += omp_in ) initializer( omp_priv = omp_orig )
#pragma omp declare reduction(+: std::vector<Force>  : omp_out += omp_in ) initializer( omp_priv = omp_orig )

///////////////////////////////////////////////////////////////////////////
// head lines;
std::vector<Vertex> setVertex_Loop_scheme(double sidex, double sidey, double l);
std::vector<Face> setFace_Loop_scheme(double sidex, double sidey, double l);
void determine_Boundary_vertex_face(double sidex, double sidey, double l, std::vector<Vertex>& vertex, std::vector<Face>& face);
void determine_Ghost_vertex_face(double sidex, double sidey, double l, bool isBoundaryFixed, bool isBoundaryPeriodic, bool isBoundaryFree, std::vector<Vertex>& vertex, std::vector<Face>& face);
void determine_AdjacentFace_for_vertex(std::vector<Vertex>& vertex, std::vector<Face>& face);
void determine_AdjacentVertex_for_vertex(std::vector<Vertex>& vertex, std::vector<Face>& face);
int Find_NodeIndex(int node1, int node2, int node3, std::vector<Vertex>& vertex);
void determine_OneRingVertex_for_face(std::vector<Vertex>& vertex, std::vector<Face>& face);
void read_struture_vertex(std::vector<Vertex>& vertex, char* filename);
std::vector<std::vector<double>> setVMU(int GaussQuadratureN);
std::vector<double> setVMUcoefficient(int GaussQuadratureN);
std::vector<std::vector<double>> determine_ShapeFunctions(std::vector<double>& vwu);
void determine_SubdivisionMatrix(std::vector<std::vector<double>>& M, std::vector<std::vector<double>>& SM1, std::vector<std::vector<double>>& SM2, std::vector<std::vector<double>>& SM3, std::vector<std::vector<double>>& SM4);

//sphere
void setsphere_Loop_scheme(std::vector<Vertex>&  vertex, std::vector<Face>&  face, double& meanl, double r, double l);
void seticosahedron(std::vector<Vertex>& vertex, std::vector<Face>& face, double r);
std::vector<int> finddot(int a, int b, std::vector<Vertex>& vertex, std::vector<Face>& face);
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

