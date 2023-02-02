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

////////////////////////
//energy force declaration


///////////////////////////////////////////////////////////////////////////
// append functions;
void determine_IsInsertionPatch_for_face(std::vector<Face>& face, std::vector<std::vector<int>>& InsertionPatch);
void determine_SpontaneousCurvature_for_face(double perturbed_c0, double unperturbed_c0, std::vector<Face>& face);
void calculate_element_area_volume(std::vector<Vertex>& vertex, std::vector<Face>& face, int subDivideTimes, std::vector<double>& GaussQuadratureCoeff, std::vector<Shapefunctions>& ShapeFunctions, SubMatrix& subMatrix);
double sum_Membrane_Area(std::vector<Face>& face);
double sum_Membrane_Volume(std::vector<Face>& face);
void clear_forceONvertex_and_energyONface(std::vector<Vertex>& vertex, std::vector<Face>& face);
void update_PreviousCoord_for_vertex(std::vector<Vertex>& vertex);
void update_reference_from_CoordPrevious(std::vector<Vertex>& vertex);
void update_PreviousForce_for_vertex(std::vector<Vertex>& vertex);
void update_PreviousEnergy_for_face(std::vector<Face>& face);
std::vector<double> ForceScale(std::vector<std::vector<double>>& Force);
double calculate_mean_force(std::vector<Vertex>& vertex);
void manage_force_for_boundary_ghost_vertex(std::vector<Vertex>& vertex, std::vector<Face>& face, Param& param);
void update_vertex_from_NonlinearConjugateGradient_s0(std::vector<Vertex>& vertex, std::vector<Face>& face, double a, std::vector<std::vector<double>>& Force, Param& param);


/////////////////////////////////////////////////////////////////////////////
// calculations
void element_energy_force_regular(std::vector<std::vector<double>>& dots, Param& param, double c0, double& Hmean, std::vector<double>& normVector, double& Ebe, std::vector<std::vector<double>>& F_be, std::vector<std::vector<double>>& F_s, std::vector<std::vector<double>>& F_v, std::vector<double>& GaussQuadratureCoeff, std::vector<Shapefunctions>& ShapeFunctions);
void element_energy_force_irregular(std::vector<std::vector<double>>& Dots, Param& param, double c0, double& Hmean, std::vector<double>& normVector, double& E_bending, std::vector<std::vector<double>>& F_be, std::vector<std::vector<double>>& F_s, std::vector<std::vector<double>>& F_v, std::vector<double>& GaussQuadratureCoeff, std::vector<Shapefunctions>& ShapeFunctions, SubMatrix subMatrix);
void energy_force_regularization(std::vector<Vertex>& vertex, std::vector<Face>& face, Param& param);
void Energy_and_Force(std::vector<Vertex>& vertex, std::vector<Face>& face, Param& param, std::vector<double>& GaussQuadratureCoeff, std::vector<Shapefunctions>& ShapeFunctions, SubMatrix& subMatrix);
double lineSearch_for_StepSize_to_minimize_energy(double a0, std::vector<std::vector<double>>& s0, std::vector<Vertex>& vertex, std::vector<Face>& face, Param& param, std::vector<double>& GaussQuadratureCoeff, std::vector<Shapefunctions>& ShapeFunctions, SubMatrix& subMatrix);
void check_nodal_force(std::vector<Vertex>& vertex, std::vector<Face>& face, Param& param, std::vector<double>& GaussQuadratureCoeff, std::vector<Shapefunctions>& ShapeFunctions, SubMatrix& subMatrix);

/////////////////////////////////////////////////////////////////////////////
// spline points
bool moveVerticesBasedOnSpline(std::vector<std::vector<double>>& splinePoints,
                               std::vector<Vertex>& vertices,
                               double lbond);
bool getAverageVector(std::vector<std::vector<double>>& inVec2D, std::vector<double>& avgVec);
std::vector<int> getClosestVertexIndex(std::vector<std::vector<double>>& splinePoints,
                                  	   std::vector<Vertex>& vertices);
double getSquaredDistance(std::vector<double>& splinePoint, Vertex& vertex);
double calculateSplineEnergyForce(std::vector<std::vector<double>>& splinePoints,
                                std::vector<Vertex>& vertices,
                                std::vector<int> splinePoints_correspondingVertexIndex,
                                double lbond, double springConst, bool doLocalSearch, double splinePointsZcoordScaling);