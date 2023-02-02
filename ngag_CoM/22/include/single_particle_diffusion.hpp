#pragma once

#include <math.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <omp.h>
#include "vector_math.hpp"
#include "parameters.hpp"

// single particles diffusing on curved surface
#pragma omp declare reduction(+: std::vector<double> : omp_out += omp_in ) initializer( omp_priv = omp_orig )
#pragma omp declare reduction(+: std::vector<std::vector<double>> : omp_out += omp_in ) initializer( omp_priv = omp_orig )

///////////////////////////////////////////////////////////////////////////
// head lines
std::vector<double> TwelveShapeFunctions(std::vector<double>& vwu);
void determine_CrossFace_for_particle(const std::vector<double>& Direction, double dL, const std::vector<double>& Coord, const std::vector<double>& PlaneNorm, std::vector<Vertex>& vertex, std::vector<Face>& face, std::vector<CrossFace>& CrossFaces, bool& hitBoundary);
void single_particles_diffusing_on_surface(std::vector<Vertex>& vertex, std::vector<Face>& face, Param& param);

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////