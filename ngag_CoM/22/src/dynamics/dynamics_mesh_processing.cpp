#include "dynamics.hpp"

using namespace std;

//calculate mesh2surface matrix
//surface = mesh2surface * mesh
void assignMesh2Surface(gsl_matrix *mesh2surface, vector<Vertex>& vertex, vector<Face>& face) {
    for (int i = 0; i < vertex.size(); i++) {
        //find adjacent face that is not ghost
        int indexAdjFace = -1; //means only ghost
        for (int jAdjFace = 0; jAdjFace < vertex[i].AdjacentFace.size(); jAdjFace++) {
            //avoid duplicate random access
            Face * adjF = &face[vertex[i].AdjacentFace[jAdjFace]];
            if ((!(adjF->IsGhost)) && (!(adjF->IsBoundary))) { //
                indexAdjFace = vertex[i].AdjacentFace[jAdjFace];
            }
        }
    if (indexAdjFace < 0) {

            gsl_matrix_set(mesh2surface, i, i, 1.0);

        } else {
            //there exists non-ghost / boundary face 
            //find vwu of vertex on face
            //!vwu assumed to be same sequence as defined in adjacent face
            Face * faceAdj = &(face[indexAdjFace]);
            vector<int> * adjVertex = &(faceAdj->AdjacentVertex);
            int vwuInd = -1;
            for (int j = 0; j < adjVertex->size(); j++) {
                if ((*adjVertex)[j] == i) {
                    vwuInd = j;
                }
            }
            //std::cout<< i << " ," << vwuAdj << endl;

            //use vwu to get sf
            vector<double> vwuAdj(3, 0.0);
            vwuAdj[vwuInd] = 1.0;
            vector<double> sfAdj = determine_ShapeFunctions(vwuAdj)[0]; //transposed sf

            //set mat(AB) (mesh2surface) to sf
            //@TODO: currently due to vwu not in order in adjacent vertex of face
            //using 0.5 / 0.0833333 directly for regular patches
            vector<int> * faceAdjOneRingVertex = &(faceAdj->OneRingVertex);
            vector<int> * iAdjVertex = &(vertex[i].AdjacentVertex);
            //regular patch
            /*
            for (int j = 0; j < 12; j++) {
                gsl_matrix_set(mesh2surface, i, (*faceAdjOneRingVertex)[j], sfAdj[j]);
            }
            */
            gsl_matrix_set(mesh2surface, i, i, 0.5);
            for (int j = 0; j < iAdjVertex->size(); j++) {
                gsl_matrix_set(mesh2surface, i, (*iAdjVertex)[j], 0.5/6);
            }
    
            
        }
        //check if ghost / boundary faces are correctly recongnized
        //std::cout << "vertex: " << i << ", indexAdjFace: " << indexAdjFace << endl;  
    }
}

//update x,y,z based on force term and random term
void next_step(vector<Vertex>& vertex, gsl_matrix *verticesProjSurface, double forceScaleConst, double randScaleConst,
                std::mt19937& gen, std::normal_distribution<>& normal_dist, bool isBoundaryPeriodic) {
    double disp = 0.0;
    double original = 0.0;
    for (int i = 0; i < vertex.size(); i++) {
        for (int j = 0; j < 3; j++) {
            disp = 0.0;
            // boundary condition
            if ((vertex[i].IsGhost || vertex[i].IsBoundary) && (!(isBoundaryPeriodic))) {
                disp = 0.0;
            } else {
                double forceterm = vertex[i].force.ForceCurvature[j] + vertex[i].force.ForceArea[j]; //get force term
                
                if (std::isnan(forceterm)) {
                    forceterm = 0.0;
                }
                double randomterm = normal_dist(gen); //get random term
                disp = forceScaleConst * forceterm + randScaleConst * randomterm;
            }

            //std::cout<<"index,"<<i<<","<<j<<",displacement"<<disp<<",B,"<<vertex[i].IsBoundary<<",G,"<<vertex[i].IsGhost<<disp[j]<<", ft: "<<forceterm<<", rt: "<<randomterm<<endl;
            original = gsl_matrix_get(verticesProjSurface, i, j);
            gsl_matrix_set(verticesProjSurface, i, j, original + disp);
        }
        //std::cout << "Force@("<<i<<"): "<< vertex[i].force.ForceCurvature[2]<<","<<vertex[i].force.ForceArea[2]<<std::endl;
    }
}