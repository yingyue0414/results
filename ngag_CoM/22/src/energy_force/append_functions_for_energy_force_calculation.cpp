#include "energy_force.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
// detailed functions

void determine_IsInsertionPatch_for_face(vector<Face>& face, vector<vector<int>>& InsertionPatch){
    #pragma omp parallel for
    for (int i = 0; i < InsertionPatch.size(); i++){
        for (int j = 0; j < InsertionPatch[0].size(); j++){
            int FaceIndex = InsertionPatch[i][j];
            face[FaceIndex].IsInsertionPatch = true;
        }
    }
}

void determine_SpontaneousCurvature_for_face(double perturbed_c0, double unperturbed_c0, vector<Face>& face){
    #pragma omp parallel for
    for (int i = 0; i < face.size(); i++){
        if ( face[i].IsBoundary == true ) 
            continue;

        if ( face[i].IsInsertionPatch == true ){
            face[i].SpontCurvature = perturbed_c0;
        }else{
            face[i].SpontCurvature = unperturbed_c0;
        }
    }
}

void calculate_element_area_volume(vector<Vertex>& vertex, vector<Face>& face, int subDivideTimes, vector<double>& GaussQuadratureCoeff, vector<Shapefunctions>& ShapeFunctions, SubMatrix& subMatrix){ 
    // five matrix used for subdivision of the irregular patch
    // M(17,11), M1(12,17), M2(12,17), M3(12,17), M4(11,17); 
    vector<vector<double>> M = subMatrix.irregM; vector<vector<double>> M1 = subMatrix.irregM1; 
    vector<vector<double>> M2 = subMatrix.irregM2; vector<vector<double>> M3 = subMatrix.irregM3; vector<vector<double>> M4 = subMatrix.irregM4;
    #pragma omp parallel for
    for (int i = 0; i < face.size(); i++){
        if ( face[i].IsBoundary == true ) // boundary faces won't contribute to the membrane area
            continue;
        
        double area = 0.0;
        double volume = 0.0;
        int numberOneRingVertex = face[i].OneRingVertex.size(); 
        if ( numberOneRingVertex == 12 ){ // regular patch 
            vector<vector<double>> dots(numberOneRingVertex,vector<double>(3)); //dots(12,3); 12 nodes
            for (int j = 0; j < numberOneRingVertex; j++){
                int NodeIndex = face[i].OneRingVertex[j];
                dots[j] = vertex[NodeIndex].Coord;
            }
            // Gaussian quadrature, 3 points
            for (int j = 0; j < GaussQuadratureCoeff.size(); j++){
                vector<vector<double>> sf = ShapeFunctions[j].sf;
                vector<double> x = sf[0] * dots;
                vector<double> a_1 = sf[1] * dots;
                vector<double> a_2 = sf[2] * dots;
                vector<double> a_3 = cross(a_1,a_2);
                double  sqa = norm(a_3); 
                double s = sqa; 
                vector<double> d = a_3 / sqa;
                area = area + 1.0/2.0 * GaussQuadratureCoeff[j] * s; 

                double v = 1.0/3.0 * s * dot(x,d);
                volume = volume + 1.0/2.0 * GaussQuadratureCoeff[j] * v; 
            }
        }else if ( numberOneRingVertex == 11 ) { // irregular patch
            vector<vector<double>> ori_dots(11,vector<double>(3,0.0));
            for ( int j = 0; j < 11; j++ ){
                int nodenum = face[i].OneRingVertex[j];
                ori_dots[j] = vertex[nodenum].Coord;
            }
            vector<vector<double>> temp (11,vector<double>(11,0.0)); for(int i=0; i<11; i++) temp[i][i] = 1.0; // identity matrix
            for ( int j = 0; j < subDivideTimes; j++){
                vector<vector<double>> newnodes17 = M * ori_dots; // 17 new nodes
                if (j != 0) {
                    temp = (M4*M) * temp;
                }
                vector<vector<double>> matrix = M*temp;

                vector<vector<double>> dots = M1*newnodes17;    // element 1
                for (int j = 0; j < GaussQuadratureCoeff.size(); j++){
                    vector<vector<double>> sf = ShapeFunctions[j].sf;
                    vector<double> x = sf[0] * dots;
                    vector<double> a_1 = sf[1] * dots;
                    vector<double> a_2 = sf[2] * dots;
                    vector<double> a_3 = cross(a_1,a_2);
                    double  sqa = norm(a_3); 
                    double s = sqa; 
                    vector<double> d = a_3 / sqa;
                    area = area + 1.0/2.0 * GaussQuadratureCoeff[j] * s; 
                    double v = 1.0/3.0 * s * dot(x,d);
                    volume = volume + 1.0/2.0 * GaussQuadratureCoeff[j] * v; 
                }
                dots = M2*newnodes17;    // element 2
                for (int j = 0; j < GaussQuadratureCoeff.size(); j++){
                    vector<vector<double>> sf = ShapeFunctions[j].sf;
                    vector<double> x = sf[0] * dots;
                    vector<double> a_1 = sf[1] * dots;
                    vector<double> a_2 = sf[2] * dots;
                    vector<double> a_3 = cross(a_1,a_2);
                    double  sqa = norm(a_3); 
                    double s = sqa; 
                    vector<double> d = a_3 / sqa;
                    area = area + 1.0/2.0 * GaussQuadratureCoeff[j] * s; 
                    double v = 1.0/3.0 * s * dot(x,d);
                    volume = volume + 1.0/2.0 * GaussQuadratureCoeff[j] * v; 
                }
                dots = M3*newnodes17;    // element 3
                for (int j = 0; j < GaussQuadratureCoeff.size(); j++){
                    vector<vector<double>> sf = ShapeFunctions[j].sf;
                    vector<double> x = sf[0] * dots;
                    vector<double> a_1 = sf[1] * dots;
                    vector<double> a_2 = sf[2] * dots;
                    vector<double> a_3 = cross(a_1,a_2);
                    double  sqa = norm(a_3); 
                    double s = sqa; 
                    vector<double> d = a_3 / sqa;
                    area = area + 1.0/2.0 * GaussQuadratureCoeff[j] * s; 
                    double v = 1.0/3.0 * s * dot(x,d);
                    volume = volume + 1.0/2.0 * GaussQuadratureCoeff[j] * v; 
                }

                vector<vector<double>> dots4 = M4 * newnodes17;   // element 4, still irregular patch
                ori_dots = dots4;
            }
        }

        face[i].ElementArea = area;
        face[i].ElementVolume = volume;
    }
}

double sum_Membrane_Area(vector<Face>& face){
    double sum = 0.0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < face.size(); i++){
        if ( face[i].IsBoundary == true ) // boundary faces won't contribute to the membrane area
            continue;
        sum += face[i].ElementArea;
    }
    return sum;
}
double sum_Membrane_Volume(vector<Face>& face){
    double sum = 0.0;
    #pragma omp parallel for reduction(+:sum)
    for (int i = 0; i < face.size(); i++){
        if ( face[i].IsBoundary == true ) // boundary faces won't contribute to the membrane area
            continue;
        sum += face[i].ElementVolume;
    }
    return sum;
}

void clear_forceONvertex_and_energyONface(vector<Vertex>& vertex, vector<Face>& face){
    #pragma omp parallel for
    for (int i = 0; i < vertex.size(); i++){
        Force Forcetmp;
        vertex[i].force = Forcetmp;
    }
    #pragma omp parallel for
    for (int i = 0; i < face.size(); i++){
        Energy Energytmp;
        face[i].energy = Energytmp;
    }
}

void update_PreviousCoord_for_vertex(vector<Vertex>& vertex){
    #pragma omp parallel for
    for (int i = 0; i < vertex.size(); i++){
        vertex[i].CoordPrevious = vertex[i].Coord;
    }
}

void update_reference_from_CoordPrevious(vector<Vertex>& vertex){
    #pragma omp parallel for 
    for (int i = 0; i < vertex.size(); i++){
        vertex[i].ReferenceCoord = vertex[i].CoordPrevious;
    }
}

void update_PreviousForce_for_vertex(vector<Vertex>& vertex){
    #pragma omp parallel for
    for (int i = 0; i < vertex.size(); i++){
        vertex[i].forcePrevious = vertex[i].force;
    }
}

void update_PreviousEnergy_for_face(vector<Face>& face){
    #pragma omp parallel for
    for (int i = 0; i < face.size(); i++){
        face[i].energyPrevious = face[i].energy;
    }
}

vector<double> ForceScale(vector<vector<double>>& Force){
    int num = Force.size();
    vector<double> out(num);
    #pragma omp parallel for
    for (int i = 0; i < num; i++){
        out[i] = norm( Force[i]);
    }
    return out;
}

double calculate_mean_force(vector<Vertex>& vertex){
    vector<double> forcescale(vertex.size(),0.0);
    #pragma omp parallel for
    for (int i = 0; i < vertex.size(); i++){
        forcescale[i] = norm( vertex[i].force.ForceTotal );
        //forcescale[i] = norm( vertex[i].Force.ForceCurvature );
        //forcescale[i] = norm( vertex[i].Force.ForceArea );
        //forcescale[i] = norm( vertex[i].Force.ForceVolume );
        //forcescale[i] = norm( vertex[i].Force.ForceRegularization );
    }
    double sum = 0.0;
    #pragma omp parallel for reduction (+:sum)
    for (int i = 0; i < vertex.size(); i++){
        sum += forcescale[i];
    }
    return sum / vertex.size();
}


void manage_force_for_boundary_ghost_vertex(vector<Vertex>& vertex, vector<Face>& face, Param& param){
    double sidex = param.sideX; double sidey = param.sideY; double l = param.l;
    int n = round(sidex/l); double dx = sidex/n; // x axis division
    double a = dx;
    double dy = sqrt(3.0)/2.0 * a; 
    int m = round(sidey/dy);      // y axis division 
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; }
    
    int vertexnum = (n+1)*(m+1);
    int facenum = m*n*2;
    if ( param.isBoundaryFixed == true ){
        Force zeroForce;
        #pragma omp parallel for 
        for (int i = 0; i < facenum; i++){
            if ( face[i].IsBoundary == false )
                continue;
            int node1 = face[i].AdjacentVertex[0];
            int node2 = face[i].AdjacentVertex[1];
            int node3 = face[i].AdjacentVertex[2];
            vertex[node1].force = zeroForce;
            vertex[node2].force = zeroForce;
            vertex[node3].force = zeroForce;
        }
    }else if ( param.isBoundaryPeriodic == true ){
        Force zeroForce;
        #pragma omp parallel for 
        for (int i = 0; i < vertexnum; i++){
            if ( vertex[i].IsGhost == true ){
                vertex[i].force = zeroForce;
            }
        }
    }else if ( param.isBoundaryFree == true ){   
        Force zeroForce;
        #pragma omp parallel for 
        for ( int i = 0; i < n+1; i++ ){
            int j = 0;
            int index = (n+1)*j + i;
            vertex[index].force = zeroForce;
            index = (n+1)*m + i;
            vertex[index].force = zeroForce;
        }
        #pragma omp parallel for
        for ( int j = 0; j < m+1; j++ ){
            int i = 0;
            int index = (n+1)*j + i;
            vertex[index].force = zeroForce;
            index = (n+1)*j + n;
            vertex[index].force = zeroForce;
        }
    }
}


void update_vertex_from_NonlinearConjugateGradient_s0(vector<Vertex>& vertex, vector<Face>& face, double a, vector<vector<double>>& Force, Param& param){
    double sidex = param.sideX; double sidey = param.sideY; double l = param.l;
    int n = round(sidex/l); double dx = sidex/n; // x axis division
    double aa = dx;
    double dy = sqrt(3.0)/2.0 * aa; 
    int m = round(sidey/dy);      // y axis division 
    if ( pow(-1.0,m) < 0.0 ) { m = m + 1; }
    
    int vertexnum = (n+1)*(m+1);
    int facenum = m*n*2;

    // update the vertex position Coord.
    #pragma omp parallel for 
    for (int i = 0; i < vertex.size(); i++){
        vertex[i].Coord = vertex[i].CoordPrevious + a * Force[i];
    }
    // deal with the boundar/ghost vertex
    if ( param.isBoundaryFixed == true ){
        #pragma omp parallel for 
        for (int i = 0; i < facenum; i++){
            if ( face[i].IsBoundary != true )
                continue;
            int node1 = face[i].AdjacentVertex[0];
            int node2 = face[i].AdjacentVertex[1];
            int node3 = face[i].AdjacentVertex[2];
            vertex[node1].Coord = vertex[node1].CoordPrevious;
            vertex[node2].Coord = vertex[node2].CoordPrevious;
            vertex[node3].Coord = vertex[node3].CoordPrevious;
        }
    }else if ( param.isBoundaryPeriodic == true ){
        #pragma omp parallel for 
        for (int i = 0; i < n; i++){
            int index0 = n*0 + i;
            int index1 = n*1 + i;
            int index2 = n*2 + i;
            int index00 = n*(m-6) + i;
            int index11 = n*(m-5) + i;
            int index22 = n*(m-4) + i;
            vertex[index0].Coord = vertex[index0].CoordPrevious + ( vertex[index00].Coord-vertex[index00].CoordPrevious );
            vertex[index1].Coord = vertex[index1].CoordPrevious + ( vertex[index11].Coord-vertex[index11].CoordPrevious );
            vertex[index2].Coord = vertex[index2].CoordPrevious + ( vertex[index22].Coord-vertex[index22].CoordPrevious );
        }
        #pragma omp parallel for 
        for (int i = 0; i < n; i++){
            int index0 = n*(m-2) + i;
            int index1 = n*(m-1) + i;
            int index2 = n*(m-0) + i;
            int index00 = n*4 + i;
            int index11 = n*5 + i;
            int index22 = n*6 + i;
            vertex[index0].Coord = vertex[index0].CoordPrevious + ( vertex[index00].Coord-vertex[index00].CoordPrevious );
            vertex[index1].Coord = vertex[index1].CoordPrevious + ( vertex[index11].Coord-vertex[index11].CoordPrevious );
            vertex[index2].Coord = vertex[index2].CoordPrevious + ( vertex[index22].Coord-vertex[index22].CoordPrevious );
        }
    }else if ( param.isBoundaryFree == true ){
        int index1, index2, index3, index4;
        // left side
        for (int j = 2; j < m ; j++ ){
            if ( pow(-1.0,j) > 0.0 ){
                index1 = (n+1)*j + 0;
                index2 = (n+1)*j + 1;
                index3 = (n+1)*(j-1) + 1;
                index4 = (n+1)*(j-1) + 2;
            }else{
                index1 = (n+1)*j + 0;
                index2 = (n+1)*j + 1;
                index3 = (n+1)*(j-1) + 0;
                index4 = (n+1)*(j-1) + 1;
            }
            vertex[index1].Coord = vertex[index2].Coord + ( vertex[index3].Coord - vertex[index4].Coord );
        }
        index1 = (n+1)*1 + 0; 
        index2 = (n+1)*2 + 0;
        index3 = (n+1)*1 + 1;
        index4 = (n+1)*2 + 1;
        vertex[index1].Coord = vertex[index2].Coord + ( vertex[index3].Coord - vertex[index4].Coord );
        // right side
        index1 = (n+1)*1 + n; 
        index2 = (n+1)*2 + n-1;
        index3 = (n+1)*1 + n-1;
        index4 = (n+1)*2 + n-2;
        vertex[index1].Coord = vertex[index2].Coord + ( vertex[index3].Coord - vertex[index4].Coord );
        for (int j = 2; j < m ; j++ ){
            if ( pow(-1.0,j) > 0.0 ){
                index1 = (n+1)*j + n;
                index2 = (n+1)*j + n-1;
                index3 = (n+1)*(j-1) + n;
                index4 = (n+1)*(j-1) + n-1;
            }else{
                index1 = (n+1)*j + n;
                index2 = (n+1)*j + n-1;
                index3 = (n+1)*(j-1) + n-1;
                index4 = (n+1)*(j-1) + n-2;
            }
            vertex[index1].Coord = vertex[index2].Coord + ( vertex[index3].Coord - vertex[index4].Coord );
        }
        // bottom 
        for (int i = 0; i < n; i++){
            index1 = (n+1)*0 + i;
            index2 = (n+1)*1 + i;
            index3 = (n+1)*1 + i+1;
            index4 = (n+1)*2 + i;
            vertex[index1].Coord = vertex[index2].Coord + ( vertex[index3].Coord - vertex[index4].Coord );
        }
        index1 = (n+1)*0 + n; 
        index2 = (n+1)*0 + n-1;
        index3 = (n+1)*1 + n;
        index4 = (n+1)*1 + n-1;
        vertex[index1].Coord = vertex[index2].Coord + ( vertex[index3].Coord - vertex[index4].Coord );
        // top 
        for (int i = 1; i < n+1; i++){
            index1 = (n+1)*m + i;
            index2 = (n+1)*(m-1) + i;
            index3 = (n+1)*(m-1) + i+1;
            index4 = (n+1)*(m-2) + i;
            vertex[index1].Coord = vertex[index2].Coord + ( vertex[index3].Coord - vertex[index4].Coord );
        }
        index1 = (n+1)*m + n; 
        index2 = (n+1)*m + n-1;
        index3 = (n+1)*(m-1) + n;
        index4 = (n+1)*(m-1) + n-1;
        vertex[index1].Coord = vertex[index2].Coord + ( vertex[index3].Coord - vertex[index4].Coord );
    } 
}


