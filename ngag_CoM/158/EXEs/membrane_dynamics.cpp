#include <math.h>
#include <stdlib.h>
#include <ctime>
#include <fstream>
#include <limits>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <omp.h>
#include <numeric> // for accumulate vector
#include <algorithm>    // std::min_element, std::max_element
#include <random>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "io.hpp"
#include "energy_force.hpp"
#include "vector_math.hpp"
#include "parameters.hpp"
#include "dynamics.hpp"
#include "gsl_matrix_methods.hpp"

using namespace std;

// basic parameters
double Radius   = 20.0;                     // sphere radius, nm
double l   = 1.0;                 // triangular side length, nm
double C0  = 1.0;                 // spontaneous curvature of insertion. Towards up is defined as positive
double c0  = 0.0;                      // spontaneous curvature of membrane
double ds  = 0.0;                      // insertion area
double ci  = 1.6;                  // area constraint coefficient
double miu = 1.0;                      // volume constraint target
double kc  = 1.0*4.17;                  // pN.nm  
double us  = 250.0;                  // pN/nm, area stretching modulus; 0.5*us*(ds)^2/s0;
double uv  = 1.0*kc;                  // coefficeint of the volume constraint, 0.5*uv*(dv)^2/v0; 0 for flat
double k   = (1.0e1)*kc;             // coefficient of the regulerization constraint, 
double K   = 1.0*k;                  // spring constant for insertion zones
double gama_shape = 0.2;
double gama_area = 0.2;
bool   isInsertionAreaConstraint = true;
double sigma = 0.0;              // 2*sigma is the lengthscale of decaying spontaneous curvature
bool   isAdditiveScheme = false; // additve scheme for the expansion of spontaneous curvature
int    GaussQuadratureN = 2; 
int    N   = 1E2;                      // total step of iteration
int    subDivideTimes = 5;       // subdivision times for the irregular patches
double CriterionForForce = 1.0e-2;
double CriterionForEnergy = 1.0e-3; 
double CriterionForArea = 1e-5;
double CriterionForRegularization = 1e-5;
bool   isGlobalConstraint = true;        // whether to use Global constraints at the beginning of the simulation
bool   isBoundaryFixed = false;
bool   isBoundaryPeriodic  = true;
bool   isBoundaryFree = false;
int    forceAreaConst = 0; // set to 0 if excluding area constraint term from force
int    forceVolumeConst = 0; // set to 0 if excluding volume constraint term from force

//maincode

void run_dynamics_flat() {

    //IO - input from parameter file
    //parameter assignment @deprecated; most parameters read from input.params
    Param param;
    param.kc = kc; param.us = us; param.k = k; param.K = K; 
    param.C0 = C0; param.c0 = c0; param.gama_shape = gama_shape; param.gama_area = gama_area; param.sigma = sigma; 
    param.GaussQuadratureN = GaussQuadratureN; param.subDivideTimes = subDivideTimes;
    param.isInsertionAreaConstraint =isInsertionAreaConstraint;
    param.isAdditiveScheme = isAdditiveScheme; param.isGlobalConstraint = isGlobalConstraint;
    param.l = l;
    param.isBoundaryFixed = isBoundaryFixed; 
    param.isBoundaryPeriodic = isBoundaryPeriodic; 
    param.isBoundaryFree = isBoundaryFree;
    param.usingNCG = true;
    param.meshpointOutput = true;
    import_file(param, "./input.params");//read param from input
    param.s0 = sqrt(3.0)/4.0*param.l*param.l;//2.0/4.0;

    vector<vector<double>> splinePoints;
    if (param.isEnergySplineIncluded) { //read spline point file if provided
        splinePoints = import_model_mesh(param.splinePointFileName);
    } 

    srand((unsigned)time(NULL)); 

    // build the triangular mesh plane. NOTE: ghost vertices and faces are included. 
    // All the boundary faces are ghost faces, which should be eliminated when output faces
    vector<Vertex> vertex = setVertex_Loop_scheme(param.sideX, param.sideY, param.l); // vertex position
    vector<Face> face = setFace_Loop_scheme(param.sideX, param.sideY, param.l); // face and its surrounding vertex
    determine_Boundary_vertex_face(param.sideX, param.sideY, param.l, vertex, face); // flag whether this vertex or face is on boundary
    determine_Ghost_vertex_face(param.sideX, param.sideY, param.l, isBoundaryFixed, isBoundaryPeriodic, isBoundaryFree, vertex, face); // flag whether this vertex or face is ghost
    determine_AdjacentFace_for_vertex(vertex, face); // find out what faces that have vertex_i, 6 or 3 or 2 or 1.
    determine_AdjacentVertex_for_vertex(vertex, face); // find out the nearby vertices for each vertex, 6 or less.
    determine_OneRingVertex_for_face(vertex, face); // find out the one-ring-vertices, 12 for flat surface with only regular patch.

    ///////////////////////////////////////
    int numvertex = vertex.size();
    int numface = face.size();
    // cout<<numvertex<<", "<<numface<<endl;
    // output the vertex and face matrix
    // add z-axis fluctutation to disgress from equilibrium state
    for (int i = 0; i < vertex.size(); i++) {
        double xi = vertex[i].Coord[0];
        double yi = vertex[i].Coord[1];
        if (xi * xi + yi * yi <= 100){
            vertex[i].Coord[2] = 10.0 - (xi * xi + yi * yi) / 10.0;
        }
        vertex[i].Coord[2] += 10.0;
    }
    

    ofstream outfile("face.csv");
    for (int i = 0; i < numface; i++) {
        outfile << face[i].AdjacentVertex[0] << ',' << face[i].AdjacentVertex[1] << ',' << face[i].AdjacentVertex[2] << '\n';
    }
    outfile.close();
    ofstream outfile1("vertex_begin.csv");
    for (int i = 0; i < numvertex; i++) {
        outfile1 << setprecision(16) << vertex[i].Coord[0] + 0.0 << ',' << vertex[i].Coord[1] + 0.0 << ',' << vertex[i].Coord[2] + 0.0 << '\n';
    }
    outfile1.close();

    ///////////////////////////////////////////////////////
    vector<vector<int> > InsertionPatch;// { {150, 151, 152, 153}};            
    determine_IsInsertionPatch_for_face(face, InsertionPatch);
    determine_SpontaneousCurvature_for_face(param.C0, param.c0, face);
    param.InsertionPatch = InsertionPatch;

    //////////////////////////////////////////////////////////
    // gauss_quadrature and shape functions
    vector<vector<double>> VWU = setVMU(param.GaussQuadratureN);
    vector<double> GaussQuadratureCoeff = setVMUcoefficient(param.GaussQuadratureN);
    vector<Shapefunctions> ShapeFunctions(VWU.size());
    for (int i = 0; i < VWU.size(); i++) {
        vector<double> vwu = VWU[i];
        ShapeFunctions[i].sf = determine_ShapeFunctions(vwu);
    }

    //@TODO: ~Dealing with irregular patch (unused in current flat membrane model)
    //irregular patch
    SubMatrix subMatrix;
    determine_SubdivisionMatrix(subMatrix.irregM, subMatrix.irregM1, subMatrix.irregM2, subMatrix.irregM3, subMatrix.irregM4);


    calculate_element_area_volume(vertex, face, subDivideTimes, GaussQuadratureCoeff, ShapeFunctions, subMatrix);
    double S0 = sum_Membrane_Area(face); // total area
    param.S0 = S0*ci;
    double V0 = sum_Membrane_Volume(face); // total volume

    ///////////////////////////////////////
    update_PreviousCoord_for_vertex(vertex);
    update_reference_from_CoordPrevious(vertex);

    //check whether the code is correct, especially whether the Force is correct!
    //check_nodal_force(vertex, face, param, GaussQuadratureCoeff, ShapeFunctions, subMatrix);
    Energy_and_Force(vertex, face, param, GaussQuadratureCoeff, ShapeFunctions, subMatrix);

    //1st energy force determines curvature and test to set spont curv = current curv
    for (int i = 0; i < face.size(); i++) {
        face[i].SpontCurvature = face[i].MeanCurvature * 20.0;
        //std::cout << "SptCurv@i=" << i << " :" << face[i].SpontCurvature<<endl;
    }
    std::cout << "======================================================" << endl;
    Energy_and_Force(vertex, face, param, GaussQuadratureCoeff, ShapeFunctions, subMatrix);
    update_PreviousCoord_for_vertex(vertex);
    update_PreviousForce_for_vertex(vertex);

    vector<double> MeanForce(N,0.0); 
    vector<double> AreaTotal(N,0.0); AreaTotal[0] = param.S;
    vector<Energy> energy_vec(N);

    vector<vector<double>> s0(vertex.size(), vector<double>(3,0.0));
    #pragma omp parallel for
    for (int i = 0; i < vertex.size(); i++) s0[i] = vertex[i].force.ForceTotal;
    vector<vector<double>> s1(vertex.size(), vector<double>(3,0.0)); 

    //generate a matrix the map between vertices
    //and their corresponding projection on limit surface
    gsl_matrix *mesh2surface = gsl_matrix_calloc(vertex.size(), vertex.size());
    //iterate vertices and calcualte mesh2surface matrix
    assignMesh2Surface(mesh2surface, vertex, face);
    //invert matrix to produce surface2mesh matrix
    gsl_matrix * surface2mesh = invert_a_matrix(mesh2surface);
    //print for TEST only
    //print_mat_contents(mesh2surface);
    //print_mat_contents(surface2mesh);

    //allocate matrix for vertices
    gsl_matrix * verticesOnMesh = gsl_matrix_alloc(vertex.size(), 3);
    copy_vertex_to_gsl_matrix(vertex, verticesOnMesh);
    //print_mat_contents(verticesOnMesh);
    gsl_matrix * verticesProjSurface = gsl_matrix_calloc(vertex.size(), 3);

    //dynamics parameter and random number generator
    double forceScaleConst = param.diffConst * param.timeStep / param.KbT;
    double randScaleConst = pow(2 * param.diffConst * param.timeStep, 0.5) * 0.0;//test override
    std::random_device rd{};//rd seed
    std::mt19937 gen{rd()};//rd generator
    std::normal_distribution<> normal_dist{0.0,1.0};

    // storing data in outfile
    if (param.meshpointOutput) {
        ofstream surfacepoint_csv("surfacepoint.csv");
        ofstream meshpoint_csv("meshpoint.csv");
        for (int i = 0; i < vertex.size(); i++) {
            surfacepoint_csv << setprecision(8) << vertex[i].Coord[0] + 0.0 << ',' << vertex[i].Coord[1] + 0.0 << ',' << vertex[i].Coord[2] + 0.0 << ',';
            meshpoint_csv << setprecision(8) << vertex[i].Coord[0] + 0.0 << ',' << vertex[i].Coord[1] + 0.0 << ',' << vertex[i].Coord[2] + 0.0 << ',';
        }
        surfacepoint_csv << '\n';
        meshpoint_csv << '\n';
        surfacepoint_csv.close();
        meshpoint_csv.close();
    }
    

    for (int iteration = 0; iteration < param.numIterations; iteration ++) {
        //1.control mesh to limit surface
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mesh2surface, verticesOnMesh,
                        0.0, verticesProjSurface);
        //verticesProjSurface = mesh2surface * verticesOnMesh;
        //print_mat_contents(verticesProjSurface);

        //2.next time step - calculate displacement on limit surface
        next_step(vertex, verticesProjSurface, forceScaleConst, randScaleConst,
                    gen, normal_dist, isBoundaryPeriodic);      

        //3.limit surface to control mesh: verticesOnMesh = surface2mesh * verticesProjSurface
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mesh2surface, verticesProjSurface,
                        0.0, verticesOnMesh);

        //4. postprocessing based on boundary condition
        //4.1. update ghost point in case of periodic boundary condition
        if (isBoundaryPeriodic) {
            postprocess_ghost_periodic(param, verticesOnMesh, vertex);
        }

        //4.2 in case of free boundary condition

        if (isBoundaryFree) {
            postprocess_ghost_free(param, verticesOnMesh, vertex, face);
        }

        //4.update values of vertex (vector of double) with verticesOnMesh (gsl matrix)
        if (param.meshpointOutput) {
            ofstream surfacepoint_csv_temp;
            ofstream meshpoint_csv_temp;
            surfacepoint_csv_temp.open("surfacepoint.csv", std::ios::app);
            meshpoint_csv_temp.open("meshpoint.csv", std::ios::app);
            for (int i = 0; i < vertex.size(); i++) {
                vector<double> realCoord(3);
                for (int j = 0; j < 3; j++) {
                    realCoord[j] = gsl_matrix_get(verticesProjSurface, i, j);
                }
                meshpoint_csv_temp << setprecision(8) << vertex[i].Coord[0] + 0.0 << ',' << vertex[i].Coord[1] + 0.0 << ',' << vertex[i].Coord[2] + 0.0 << ',';
                surfacepoint_csv_temp << setprecision(8) << realCoord[0] + 0.0 << ',' << realCoord[1] + 0.0 << ',' << realCoord[2] + 0.0 << ',';
            }
            meshpoint_csv_temp << '\n';
            surfacepoint_csv_temp << '\n';
            meshpoint_csv_temp.close();
            surfacepoint_csv_temp.close();
        }
        
        Energy_and_Force(vertex, face, param, GaussQuadratureCoeff, ShapeFunctions, subMatrix);
        cout<<"===================ITERATION:" << iteration << "===========================" << endl;

    }

    //currently running python3 code to convert to xyz file
    if (param.xyzOutput) {
        //@TODO add functionality
    }
}

int main() {
    //ProfilerStart("output.prof");
    run_dynamics_flat();
    //ProfilerStop();
}
