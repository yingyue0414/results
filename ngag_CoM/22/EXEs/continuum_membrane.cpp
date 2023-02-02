#include <math.h>
#include <stdlib.h>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <omp.h>
#include <numeric> // for accumulate vector
#include <algorithm>    // std::min_element, std::max_element

#include "io.hpp"
#include "energy_force.hpp"
#include "vector_math.hpp"
#include "parameters.hpp"
#include "gsl_matrix_methods.hpp"

using namespace std;
// basic parameters
double sideX  = 180.0;            // rectangle sidelength x, nm
double sideY  = sideX;            // rectangle sidelength y, nm
double Radius   = 20.0;                     // sphere radius, nm
double l   = 2.0;                 // triangular side length, nm
double C0  = 0.8;                 // spontaneous curvature of insertion. Towards up is defined as positive
double c0  = 0.0;                      // spontaneous curvature of membrane
double ds  = 0.0;                      // insertion area
double ci  = 1.0;                  // area constraint coefficient
double miu = 0.65;                      // volume constraint target
double kc  = 20*4.17;                  // pN.nm  
double us  = 250.0;                  // pN/nm, area stretching modulus; 0.5*us*(ds)^2/s0;
double uv  = 0.0;                  // coefficeint of the volume constraint, 0.5*uv*(dv)^2/v0; 0 for flat
double k   = (1.0e1)*kc;             // coefficient of the regulerization constraint, 
double K   = 1.0*k;                  // spring constant for insertion zones
double gama_shape = 0.2;
double gama_area = 0.2;
bool   isInsertionAreaConstraint = true;
double sigma = 0.0;              // 2*sigma is the lengthscale of decaying spontaneous curvature
bool   isAdditiveScheme = false; // additve scheme for the expansion of spontaneous curvature
int    GaussQuadratureN = 2; 
int    N   = 1E5;                      // total step of iteration
int    subDivideTimes = 5;       // subdivision times for the irregular patches
double CriterionForForce = 1.0e-2;
double CriterionForEnergy = 1.0e-3; 
double CriterionForArea = 1e-5;
double CriterionForRegularization = 1e-5;
bool   isGlobalConstraint = true;        // whether to use Global constraints at the beginning of the simulation
bool   isBoundaryFixed = false;
bool   isBoundaryPeriodic  = true;
bool   isBoundaryFree = false;




/////////////////////////////////////////////////////////////////////////////////
// main code



void run_flat() {

    srand((unsigned)time(NULL)); 

    // build the triangular mesh plane. NOTE: ghost vertices and faces are included. 
    // All the boundary faces are ghost faces, which should be eliminated when output faces
    vector<Vertex> vertex = setVertex_Loop_scheme(sideX, sideY, l); // vertex position
    vector<Face> face = setFace_Loop_scheme(sideX, sideY, l); // face and its surrounding vertex
    determine_Boundary_vertex_face(sideX, sideY, l, vertex, face); // flag whether this vertex or face is on boundary
    determine_Ghost_vertex_face(sideX, sideY, l, isBoundaryFixed, isBoundaryPeriodic, isBoundaryFree, vertex, face); // flag whether this vertex or face is ghost
    determine_AdjacentFace_for_vertex(vertex, face); // find out what faces that have vertex_i, 6 or 3 or 2 or 1.
    determine_AdjacentVertex_for_vertex(vertex, face); // find out the nearby vertices for each vertex, 6 or less.
    determine_OneRingVertex_for_face(vertex, face); // find out the one-ring-vertices, 12 for flat surface with only regular patch.

    //////////////////////////////////////////////////////
    ///////////////////////////////////////
    // read a structure file
    //char name[32] = "vertex_read.csv";
    //read_struture_vertex(vertex, name);
    ///////////////////////////////////////
    int numvertex = vertex.size();
    int numface = face.size();
    // cout<<numvertex<<", "<<numface<<endl;
    // output the vertex and face matrix
    ofstream outfile("face.csv");
    for (int i = 0; i < numface; i++) {
        outfile << face[i].AdjacentVertex[0] + 1 << ',' << face[i].AdjacentVertex[1] + 1 << ',' << face[i].AdjacentVertex[2] + 1 << '\n';
    }
    outfile.close();
    ofstream outfile1("vertex_begin.csv");
    for (int i = 0; i < numvertex; i++) {
        outfile1 << setprecision(16) << vertex[i].Coord[0] + 0.0 << ',' << vertex[i].Coord[1] + 0.0 << ',' << vertex[i].Coord[2] + 0.0 << '\n';
    }
    outfile1.close();
    ///////////////////////////////////////////////////////
    vector<vector<int> > InsertionPatch; //{{1418, 1419, 1420, 1421, 1422, 1423, 1424, 1425, 1426,
                                        //1484, 1485, 1486, 1487, 1488, 1489, 1490, 1491, 1492,
                                        //  1550, 1551, 1552, 1553, 1554, 1555, 1556, 1557, 1558}};        
    determine_IsInsertionPatch_for_face(face, InsertionPatch);
    determine_SpontaneousCurvature_for_face(C0, c0, face);
    //////////////////////////////////////////////////////////
    // gauss_quadrature and shape functions
    vector<vector<double>> VWU = setVMU(GaussQuadratureN);
    vector<double> GaussQuadratureCoeff = setVMUcoefficient(GaussQuadratureN);
    vector<Shapefunctions> ShapeFunctions(VWU.size());
    for (int i = 0; i < VWU.size(); i++) {
        vector<double> vwu = VWU[i];
        ShapeFunctions[i].sf = determine_ShapeFunctions(vwu);
    }
    SubMatrix subMatrix;
    determine_SubdivisionMatrix(subMatrix.irregM, subMatrix.irregM1, subMatrix.irregM2, subMatrix.irregM3, subMatrix.irregM4);
    //////////////////////////////////////////////////////////
    //calculate_element_area_volume(vertex, face, GaussQuadratureCoeff, ShapeFunctions); 
    // calculate the elemental area 
    calculate_element_area_volume(vertex, face, subDivideTimes, GaussQuadratureCoeff, ShapeFunctions, subMatrix);
    double S0 = sum_Membrane_Area(face); // total area
    double V0 = sum_Membrane_Volume(face); // total volume

    Param param;
    param.kc = kc; param.us = us; param.k = k; param.K = K; param.S0 = S0; param.V0 = V0;
    param.C0 = C0; param.c0 = c0; param.gama_shape = gama_shape; param.gama_area = gama_area; param.sigma = sigma; 
    param.GaussQuadratureN = GaussQuadratureN; param.subDivideTimes = subDivideTimes;
    param.isInsertionAreaConstraint =isInsertionAreaConstraint;
    param.isAdditiveScheme = isAdditiveScheme; param.isGlobalConstraint = isGlobalConstraint;
    param.l = l;
    param.isBoundaryFixed = isBoundaryFixed; 
    param.isBoundaryPeriodic = isBoundaryPeriodic; 
    param.isBoundaryFree = isBoundaryFree;
    param.usingNCG = true;
    param.InsertionPatch = InsertionPatch;
    
    
    //for testing only; add spline points to energy minimization
    //update the positions by moving the membrane upward (vertices)
    //override vertex_begin.csv afterwards
    double lbond = 12.0;
    double curr_zdisplacement = 0.0;
    param.lbond = lbond;

    import_file(param, "./input.params");//override param from input
    
    vector<double> planeParamsVec; //ax + by + c = z
    if (param.isEnergySplineIncluded) { //read spline point file if provided
        param.splinePoints = import_model_mesh(param.splinePointFileName);
        //move spline points upwards until the lower boundary is approximately z=0
        for (int i = 0; i < param.splinePoints.size(); i++){
            for (int j = 0; j < param.splinePoints[i].size(); j++){
                param.splinePoints[i][j] += 10.0;//move spline by 18.0
            }
        }
        moveVerticesBasedOnSpline(param.splinePoints, vertex, lbond); //move vertices upwards
        ofstream outfile2("vertex_begin.csv");
        for (int i = 0; i < numvertex; i++) {
            outfile2 << setprecision(16) << vertex[i].Coord[0] + 0.0 << ',' << vertex[i].Coord[1] + 0.0 << ',' << vertex[i].Coord[2] + 0.0 << '\n';
        }
        outfile2.close();
        //find closest vertex point and save in vector
        param.splinePoints_correspondingVertexIndex = getClosestVertexIndex(param.splinePoints, vertex);
    }
    
    param.s0 = sqrt(3.0)/4.0*param.l*param.l;//2.0/4.0;

    ///////////////////////////////////////
    update_PreviousCoord_for_vertex(vertex);
    update_reference_from_CoordPrevious(vertex);

    // check whether the code is correct, especially whether the Force is correct!
    //check_nodal_force(vertex, face, param, GaussQuadratureCoeff, ShapeFunctions, subMatrix);

    Energy_and_Force(vertex, face, param, GaussQuadratureCoeff, ShapeFunctions, subMatrix);
    /*if (param.isEnergySplineIncluded){
        param.energy.energySpline = calculateSplineEnergyForce(param.splinePoints, vertex, param.splinePoints_correspondingVertexIndex,
                            param.lbond, param.springConst, false, param.splinePointsZcoordScaling);
        param.energy.energyTotal += param.energy.energySpline;
    }*/
    update_PreviousCoord_for_vertex(vertex);
    update_PreviousForce_for_vertex(vertex);

    //single_particles_diffusing_on_surface(vertex, face, param); exit(0);

    vector<double> MeanForce(N,0.0); 
    vector<double> AreaTotal(N,0.0); AreaTotal[0] = param.S;
    vector<Energy> energy_vec(N);

    vector<vector<double>> s0(vertex.size(), vector<double>(3,0.0));
    #pragma omp parallel for
    for (int i = 0; i < vertex.size(); i++) s0[i] = vertex[i].force.ForceTotal;
    vector<vector<double>> s1(vertex.size(), vector<double>(3,0.0)); 
    
    vector<int> timesOffNCG(N,0);
    bool IsCriteriaSatisfied = false; 
    bool updateReference = true;
    int iteration = 0;
    double TrialStepSize = 0.0;
    while ( IsCriteriaSatisfied == false && iteration < N-1){
        // The step size a0 needs a trial value, which is determined by rule-of-thumb. 
        if ( iteration == 0 || iteration % 50 == 0 ){
            vector<double> Scale_s0 = ForceScale(s0);
            double MaxForceScale = *max_element(Scale_s0.begin(), Scale_s0.end());
//test
std::cout<<"max_force_scale="<<MaxForceScale<<", energy="<<param.energy.energyTotal<<endl; 
//exit(0);
            TrialStepSize = l / MaxForceScale; // a0 = 1;
        }else{
            TrialStepSize = TrialStepSize * 2e1;
        } 
        {
            double StepSize = lineSearch_for_StepSize_to_minimize_energy(TrialStepSize, s0, vertex, face, param, GaussQuadratureCoeff, ShapeFunctions, subMatrix);
            cout<<"step: "<<iteration<<", trial StepSize = "<<TrialStepSize<<", StepSize = "<<StepSize<<endl;
            if ( StepSize == -1 ){
               cout<<"step: "<<iteration<<". Note: no efficent step size a is found. Stop now! "<<endl;
               //printoutREF(vertexref);
               //printoutstuck(vertex0);
               break;
            }
            update_vertex_from_NonlinearConjugateGradient_s0(vertex, face, StepSize, s0, param); 
            if (iteration % 200 == 0 && param.springConst <= 1000.0){
                param.springConst *= 3;
            }
            // calculate the new Force and Energy
            Energy_and_Force(vertex, face, param, GaussQuadratureCoeff, ShapeFunctions, subMatrix);
            /*if (param.isEnergySplineIncluded){
                param.energy.energySpline = calculateSplineEnergyForce(param.splinePoints, vertex, param.splinePoints_correspondingVertexIndex,
                                 param.lbond, param.springConst, false, param.splinePointsZcoordScaling);
                param.energy.energyTotal += param.energy.energySpline;
            }*/
            update_PreviousCoord_for_vertex(vertex); /// update the previous position
            update_PreviousForce_for_vertex(vertex);
            update_PreviousEnergy_for_face(face);
            param.energyPrevious = param.energy;
            //////////////////////////////////////////////////////////////////////////////
            // calculate the direction s
            double NCGfactor0 = 0.0;
            #pragma omp parallel for reduction(+:NCGfactor0) 
            for (int i = 0; i < vertex.size(); i++){
                vector<double> Forcetmp = - vertex[i].forcePrevious.ForceTotal; 
                NCGfactor0 = NCGfactor0 + dot(Forcetmp, Forcetmp);
            }
            double NCGfactor = 0.0;
            #pragma omp parallel for reduction(+:NCGfactor) 
            for (int i = 0; i < vertex.size(); i++){
                vector<double> Forcetmp = - vertex[i].force.ForceTotal; 
                NCGfactor = NCGfactor + dot(Forcetmp, Forcetmp);
            }
            double peta1 = NCGfactor / NCGfactor0; if (param.usingNCG == false) peta1 = 0.0;
            #pragma omp parallel for 
            for (int i = 0; i < vertex.size(); i++){
                s1[i] = vertex[i].force.ForceTotal + peta1 * s0[i];  
            }
            //////////////////////////////////////////////////////////////////////////
            // update the direction for Nonlinear Conjugate gradient method
            s0 = s1;
            TrialStepSize = StepSize;
            
            // store the Energy and nodal Force
            energy_vec[iteration] = param.energy;
            MeanForce[iteration] = calculate_mean_force(vertex); 
            AreaTotal[iteration] = param.S;
            // update the reference 
            {
                if ( abs(energy_vec[iteration].energyTotal - energy_vec[iteration-1].energyTotal) < 1e-3 || abs(MeanForce[iteration]-MeanForce[iteration-1]) < 1e-3 ){
                   cout<<"update the reference structure! "<<endl;
                   update_reference_from_CoordPrevious(vertex);
                }
                //increase spring const
                /*if (iteration % 200 == 0 && param.springConst <= 1000.0){
                    param.springConst *= 3;
                }*/
            }

            ////////////////////////////////////////////////////////
            // to check whether NCG has been stucked for too many times consistently. if so, change the flag to not-using-NCG 
            if ( param.isNCGstucked == true ){
                timesOffNCG[iteration] = 1;
            }
            if ( iteration > 10 ){
                int NoffNCG = accumulate(timesOffNCG.begin()+iteration-9, timesOffNCG.begin()+iteration, 0);
                if ( NoffNCG > 8 ){
                    param.usingNCG = false;
                }else if ( NoffNCG < 3 ){
                    param.usingNCG = true;
                    param.isNCGstucked = false;
                }
            }
        }
        // output parameters     
        cout<<"step: "<< iteration <<". Energy= "<<energy_vec[iteration].energyTotal<<". meanF= "<<MeanForce[iteration]<<endl;
        //////////////////////////////////////////////////////////////////////////////////
        // check whether to stop. if the total Energy is flat for 100 simulation steps, then stop
        if ( iteration > 500 && abs((energy_vec[iteration].energyTotal-energy_vec[iteration-500].energyTotal)/500) < CriterionForEnergy && MeanForce[iteration] < CriterionForForce ){
            cout<<"The Energy is minimized. Stop now!"<<endl;
            IsCriteriaSatisfied = true;
            break;
        }
        //double Ere_vs_Etot = Energy[iteration].EnergyRegularization / Energy[iteration].EnergyTolt;

        //////////////////////////////////////////////////////////////////////////////////
        // output the vertex
        if ( iteration % 100 == 0 ){    
            int kk = iteration/100;
            char filename[20] = "vertex%d.csv";
            sprintf(filename,"vertex%d.csv",kk);
            ofstream outfile(filename);
            for (int j = 0; j < vertex.size(); j++) {
                outfile << vertex[j].Coord[0] << ',' << vertex[j].Coord[1] << ',' << vertex[j].Coord[2] << '\n';
            }
            outfile.close();
        }
        // output Energy and meanforce
        ofstream outfile2("EnergyForce.csv"); 
        for (int j = 0; j <= iteration; j++){
            if ( j == 0 ){
                outfile2 <<"Energy-Curvature, -Area, -Regularization, -Total ((pN.nm)); Mean Force (pN)"<< '\n';
            }
            outfile2 << energy_vec[j].energyCurvature << ", " << energy_vec[j].energyArea << ", " << energy_vec[j].energyRegularization << ", " << energy_vec[j].energyTotal<< ", " << MeanForce[j] << '\n';
        }
        outfile2.close();

        iteration ++;
    }
    // output the final structure
    ofstream outfile33("vertexfinal.csv");
    for (int j = 0; j < vertex.size(); j++) {
        outfile33 << setprecision(16) << vertex[j].CoordPrevious[0] << ',' << vertex[j].CoordPrevious[1] << ',' << vertex[j].CoordPrevious[2] << '\n';
    }
    outfile33.close();
}

void run_sphere() {

	srand((unsigned)time(NULL)); 

    // build the triangular mesh plane. NOTE: ghost vertices and faces are included. 
    // All the boundary faces are ghost faces, which should be eliminated when output faces
    double meanL = 0.0;     // the mean value of the triangular side length.
    vector<Vertex> vertex(12); // vertex 
    vector<Face> face(20);    // face and its surrounding vertex
    setsphere_Loop_scheme(vertex, face, meanL, Radius, l);

    determine_AdjacentFace_for_vertex(vertex, face); // find out what faces that have vertex_i, 6 or 5.
    determine_AdjacentVertex_for_vertex(vertex, face); // find out the nearby vertices for each vertex, 6 or 5.
    determine_OneRingVertex_for_face(vertex, face); // find out the one-ring-vertices, 12 or 11.

    //////////////////////////////////////////////////////
    ///////////////////////////////////////
    // read a structure file
    //char name[32] = "vertex_readR10prolate.csv";
    //char name[32] = "vertex_readR20.csv";
    //read_struture_vertex(vertex, name);
    ///////////////////////////////////////
    int numvertex = vertex.size();
    int numface = face.size();
    // cout<<numvertex<<", "<<numface<<endl;
    // output the vertex and face matrix
    
    ofstream outfile("face.csv");
    for (int i = 0; i < numface; i++) {
        outfile << face[i].AdjacentVertex[0] + 1 << ',' << face[i].AdjacentVertex[1] + 1 << ',' << face[i].AdjacentVertex[2] + 1 << '\n';
    }
    outfile.close();
    ofstream outfile1("vertex_begin.csv");
    for (int i = 0; i < numvertex; i++) {
        outfile1 << setprecision(16) << vertex[i].Coord[0] + 0.0 << ',' << vertex[i].Coord[1] + 0.0 << ',' << vertex[i].Coord[2] + 0.0 << '\n';
    }
    outfile1.close();
    
    ///////////////////////////////////////////////////////
    vector<vector<int> > InsertionPatch; //{ {6938, 912, 6941, 422},
                                        //  {6421, 652, 6424, 173} };            
    determine_IsInsertionPatch_for_face(face, InsertionPatch);
    determine_SpontaneousCurvature_for_face(C0, c0, face);
    //////////////////////////////////////////////////////////
    // gauss_quadrature and shape functions
    vector<vector<double>> VWU = setVMU(GaussQuadratureN);
    vector<double> GaussQuadratureCoeff = setVMUcoefficient(GaussQuadratureN);
    vector<Shapefunctions> ShapeFunctions(VWU.size());
    for (int i = 0; i < VWU.size(); i++) {
        vector<double> vwu = VWU[i];
        ShapeFunctions[i].sf = determine_ShapeFunctions(vwu);
    }
    SubMatrix subMatrix;
    determine_SubdivisionMatrix(subMatrix.irregM, subMatrix.irregM1, subMatrix.irregM2, subMatrix.irregM3, subMatrix.irregM4);
    //////////////////////////////////////////////////////////
    calculate_element_area_volume(vertex, face, subDivideTimes, GaussQuadratureCoeff, ShapeFunctions, subMatrix); // calculate the elemental area 
    double S0 = sum_Membrane_Area(face); // total area
    //double V0 = sum_Membrane_Volume(face); // total volume
    double RadiusReal = sqrt(S0/4.0/M_PI);
    double V0 = 4.0/3.0 * M_PI * pow(RadiusReal,3.0);
    cout<<"R = "<<RadiusReal<<". S0 = "<<S0*ci<<", V0 = "<<V0*miu<<endl;

    Param param;
    param.kc = kc; param.us = us; param.uv = uv; param.k = k; param.K = K; 
    param.C0 = C0; param.c0 = c0; param.gama_shape = gama_shape; param.gama_area = gama_area; param.sigma = sigma; 
    param.GaussQuadratureN = GaussQuadratureN; param.subDivideTimes = subDivideTimes;
    param.isInsertionAreaConstraint =isInsertionAreaConstraint;
    param.isAdditiveScheme = isAdditiveScheme; param.isGlobalConstraint = isGlobalConstraint;
    param.s0 = sqrt(3.0)/4.0*l*l;//2.0/4.0; // /insertionpatch.n_cols; 
    param.S0 = S0 * ci; 
    param.V0 = V0 * miu;
    param.Radius = RadiusReal; param.l = l;
    param.isBoundaryFixed = isBoundaryFixed; 
    param.isBoundaryPeriodic = isBoundaryPeriodic; 
    param.isBoundaryFree = isBoundaryFree;
    param.usingNCG = true;
    param.InsertionPatch = InsertionPatch;
    ///////////////////////////////////////
    update_PreviousCoord_for_vertex(vertex);
    update_reference_from_CoordPrevious(vertex);

    // check whether the code is correct, especially whether the force is correct!
    //check_nodal_force(vertex, face, param, GaussQuadratureCoeff, ShapeFunctions, subMatrix);
    Energy_and_Force(vertex, face, param, GaussQuadratureCoeff, ShapeFunctions, subMatrix);
    update_PreviousCoord_for_vertex(vertex);
    update_PreviousForce_for_vertex(vertex);
    //single_particles_diffusing_on_surface(vertex, face, param); exit(0);

    vector<double> MeanForce(N,0.0); 
    vector<double> AreaTotal(N,0.0); AreaTotal[0] = param.S;
    vector<Energy> Energy(N);
 
    vector<vector<double>> s0(vertex.size(), vector<double>(3,0.0));
    #pragma omp parallel for
    for (int i = 0; i < vertex.size(); i++) s0[i] = vertex[i].force.ForceTotal;
    vector<vector<double>> s1(vertex.size(), vector<double>(3,0.0)); 
    
    vector<int> timesOffNCG(N,0);
    bool IsCriteriaSatisfied = false; 
    bool updateReference = true;
    int iteration = 0;
    double TrialStepSize = 0.0;
    while ( IsCriteriaSatisfied == false && iteration < N-1){
        // The step size a0 needs a trial value, which is determined by rule-of-thumb. 
        if ( iteration == 0 || iteration % 50 == 0 ){
            vector<double> Scale_s0 = ForceScale(s0);
            double MaxForceScale = *max_element(Scale_s0.begin(), Scale_s0.end());
            TrialStepSize = l / MaxForceScale; // a0 = 1;
        }else{
            TrialStepSize = TrialStepSize * 2e1;
        }
        cout<<"beginning of step: "<< iteration <<endl;
        {
            double StepSize = lineSearch_for_StepSize_to_minimize_energy(TrialStepSize, s0, vertex, face, param, GaussQuadratureCoeff, ShapeFunctions, subMatrix);
            cout<<"step: "<<iteration<<", trial StepSize = "<<TrialStepSize<<", StepSize = "<<StepSize<<endl;
            if ( StepSize == -1 ){
               cout<<"step: "<<iteration<<". Note: no efficent step size a is found. Stop now! "<<endl;
               //printoutREF(vertexref);
               //printoutstuck(vertex0);
               break;
            }
            //update_vertex_from_NonlinearConjugateGradient_s0(vertex, face, StepSize, s0, param);
            #pragma omp parallel for 
            for (int i = 0; i < vertex.size(); i++) { vertex[i].Coord = vertex[i].CoordPrevious + StepSize * s0[i]; } 
            // calculate the new force and energy
            Energy_and_Force(vertex, face, param, GaussQuadratureCoeff, ShapeFunctions, subMatrix);
            update_PreviousCoord_for_vertex(vertex); /// update the previous position
            update_PreviousForce_for_vertex(vertex);
            update_PreviousEnergy_for_face(face);
            param.energyPrevious = param.energy;
            //////////////////////////////////////////////////////////////////////////////
            // calculate the direction s
            double NCGfactor0 = 0.0;
            vector<double> Forcetmp(3);
            #pragma omp parallel for reduction(+:NCGfactor0)
            for (int i = 0; i < vertex.size(); i++){
                negative(vertex[i].forcePrevious.ForceTotal, Forcetmp); 
                NCGfactor0 += dot(Forcetmp, Forcetmp);
            }
            double NCGfactor = 0.0;
            #pragma omp parallel for reduction(+:NCGfactor) 
            for (int i = 0; i < vertex.size(); i++){
                negative(vertex[i].force.ForceTotal, Forcetmp); 
                NCGfactor += dot(Forcetmp, Forcetmp);
            }
            double peta1 = NCGfactor / NCGfactor0; if (param.usingNCG == false) peta1 = 0.0;
            #pragma omp parallel for 
            for (int i = 0; i < vertex.size(); i++){
                s1[i] = vertex[i].force.ForceTotal + peta1 * s0[i];  
            }
            //////////////////////////////////////////////////////////////////////////
            // update the direction for Nonlinear Conjugate gradient method
            s0 = s1;
            TrialStepSize = StepSize;
            // store the energy and nodal force
            Energy[iteration] = param.energy;
            MeanForce[iteration] = calculate_mean_force(vertex); 
            AreaTotal[iteration] = param.S;
            // update the reference 
            if (iteration > 1){
               if ( abs(Energy[iteration].energyTotal - Energy[iteration-1].energyTotal) < 1e-3 || abs(MeanForce[iteration]-MeanForce[iteration-1]) < 1e-3 ){
                   cout<<"update the reference structure! "<<endl;
                   update_reference_from_CoordPrevious(vertex);
               }
            }
            ////////////////////////////////////////////////////////
            // to check whether NCG has been stucked for too many times consistently. if so, change the flag to not-using-NCG 
            if ( param.isNCGstucked == true ){
                timesOffNCG[iteration] = 1;
            }
            if ( iteration > 10 ){
                int NoffNCG = accumulate(timesOffNCG.begin()+iteration-9, timesOffNCG.begin()+iteration, 0);
                if ( NoffNCG > 8 ){
                    param.usingNCG = false;
                }else if ( NoffNCG < 3 ){
                    param.usingNCG = true;
                    param.isNCGstucked = false;
                }
            }
        }
        // output parameters     
        cout<<"step: "<< iteration <<". Energy= "<<Energy[iteration].energyTotal<<". meanF= "<<MeanForce[iteration]<<endl;
        //////////////////////////////////////////////////////////////////////////////////
        // check whether to stop. if the total energy is flat for 100 simulation steps, then stop
        if ( iteration > 500 && abs((Energy[iteration].energyTotal-Energy[iteration-500].energyTotal)/500) < CriterionForEnergy && MeanForce[iteration] < CriterionForForce ){
            cout<<"The energy is minimized. Stop now!"<<endl;
            IsCriteriaSatisfied = true;
            break;
        }
        //double Ere_vs_Etot = Energy[iteration].EnergyRegularization / Energy[iteration].EnergyTolt;

        //////////////////////////////////////////////////////////////////////////////////
        // output the vertex
        if ( iteration % 100 == 0 ){    
            int kk = iteration/100;
            char filename[20] = "vertex%d.csv";
            sprintf(filename,"vertex%d.csv",kk);
            ofstream outfile(filename);
            for (int j = 0; j < vertex.size(); j++) {
                outfile << vertex[j].Coord[0] << ',' << vertex[j].Coord[1] << ',' << vertex[j].Coord[2] << '\n';
            }
            outfile.close();
        }
        // output energy and meanforce
        ofstream outfile2("EnergyForce.csv"); 
        for (int j = 0; j <= iteration; j++){
            if ( j == 0 ){
                outfile2 <<"Energy-Curvature, -Area, -Regularization, -Total ((pN.nm)); Mean Force (pN)"<< '\n';
            }
            outfile2 << Energy[j].energyCurvature << ", " << Energy[j].energyArea << ", " << Energy[j].energyRegularization << ", " << Energy[j].energyTotal<< ", " << MeanForce[j] << '\n';
        }
        outfile2.close();

        iteration ++;
    }
    // output the final structure
    ofstream outfile33("vertexfinal.csv");
    for (int j = 0; j < vertex.size(); j++) {
        outfile33 << setprecision(16) << vertex[j].CoordPrevious[0] << ',' << vertex[j].CoordPrevious[1] << ',' << vertex[j].CoordPrevious[2] << '\n';
    }
    outfile33.close();
}

int main() {
    //ProfilerStart("output.prof");
    run_flat();
    //ProfilerStop();
}
