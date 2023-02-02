#include "energy_force.hpp"

using namespace std;





///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void element_energy_force_regular(vector<vector<double>>& dots, Param& param, double c0, double& Hmean, vector<double>& normVector, double& Ebe, vector<vector<double>>& F_be, vector<vector<double>>& F_s, vector<vector<double>>& F_v, vector<double>& GaussQuadratureCoeff, vector<Shapefunctions>& ShapeFunctions) {
    // F_be is the Force related to the curvature; F_s is the Force related to the area-constraint; F_v is the Force related to the volume-constraint.
    // All these forces should be 0 as input.
    // initialize output parameters
    Ebe = 0.0;
    Hmean = 0.0;
    //////////////////////////////////////////////////////////////
    double kc = param.kc;
    double us = param.us/param.S0;
    double S0 = param.S0;
    double S  = param.S;
    double uv = param.uv/param.V0;
    double V0 = param.V0;
    double V  = param.V;
    /////////////////////////////////////////////////////////////
    //double definition
    double sqa = 0.0;
    double sqa_1 = 0.0;
    double sqa_2 = 0.0;

    double H_curv = 0.0;
    double ebe = 0.0;
    //vector definition
    
   
    vector<double> x(3);
    vector<double> a_1(3);
    vector<double> a_2(3);
    vector<double> a_11(3);
    vector<double> a_22(3);
    vector<double> a_12(3);
    vector<double> a_21(3);
    vector<double> xa(3);
   
    vector<double> xa_1(3);
    vector<double> xa_2(3);

    vector<double> a_3(3);
    vector<double> a_31(3);
    vector<double> a_32(3);
    vector<double> d(3);
    vector<double> d_1(3);
    vector<double> d_2(3);
    vector<double> a1(3);
    vector<double> a2(3);
    vector<double> a11(3);
    vector<double> a12(3);
    vector<double> a21(3);
    vector<double> a22(3);

    vector<double> n1_be(3);
    vector<double> n2_be(3);
    vector<double> m1_be(3);
    vector<double> m2_be(3);
    vector<double> n1_cons(3);
    vector<double> n2_cons(3);
    vector<double> n1_conv(3);
    vector<double> n2_conv(3);

    //2d vectors
    vector<vector<double>> f_be(12, vector<double>(3,0.0)); 
    vector<vector<double>> f_cons(12, vector<double>(3,0.0)); 
    vector<vector<double>> f_conv(12, vector<double>(3,0.0)); 
    vector<vector<double>> da1(3, vector<double>(3,0.0));
    vector<vector<double>> da2(3, vector<double>(3,0.0));

    //TEMP vectors / double
    vector<double> tmp_f(3);
    vector<double> tmp_l(3);
    vector<double> tmp_sum(3);
    double tmp_const_f = 0.0;
    double tmp_const_l = 0.0;
    
    // Gaussian quadrature, second-order or 3 points.
    for (int i = 0; i < GaussQuadratureCoeff.size(); i++) {
        vector<vector<double>> sf = ShapeFunctions[i].sf; //12 shape functions
        // a_1,2,3 covariant vectors; a1,2 contravariant vectors;
        // a_11: a_1 differential to v; a_12: a_1 differential to w;
        rowvec_matrix_multiplication(sf[0], dots, x);
        rowvec_matrix_multiplication(sf[1], dots, a_1);
        rowvec_matrix_multiplication(sf[2], dots, a_2);
        rowvec_matrix_multiplication(sf[3], dots, a_11);
        rowvec_matrix_multiplication(sf[4], dots, a_22);
        rowvec_matrix_multiplication(sf[5], dots, a_12);
        rowvec_matrix_multiplication(sf[6], dots, a_21);
        cross(a_1, a_2, xa);
        sqa = norm(xa);

        //vector<double> xa_1 = cross(a_11,a_2) + cross(a_1,a_21);
        a_cross_b_plus_c_cross_d(a_11, a_2, a_1, a_21, tmp_f, tmp_l, xa_1);
        //vector<double> xa_2 = cross(a_12,a_2) + cross(a_1,a_22);
        a_cross_b_plus_c_cross_d(a_12, a_2, a_1, a_22, tmp_f, tmp_l, xa_2);

        sqa_1 = 1.0/sqa * dot(xa, xa_1);
        sqa_2 = 1.0/sqa * dot(xa, xa_2);
        //vector<double> a_3 = xa/sqa;
        const_division(xa, sqa, a_3);
        //vector<double> a_31 = 1.0/sqa/sqa * (xa_1*sqa - xa*sqa_1);
        const_multiplication(xa_1, sqa, tmp_f);
        const_multiplication(xa, sqa_1, tmp_l);
        subtraction(tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_sum, 1.0/sqa/sqa, a_31);
        //vector<double> a_32 = 1.0/sqa/sqa * (xa_2*sqa - xa*sqa_2);
        const_multiplication(xa_2, sqa, tmp_f);
        const_multiplication(xa, sqa_2, tmp_l);
        subtraction(tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_sum, 1.0/sqa/sqa, a_32);
        
        d = a_3;
        d_1 = a_31;
        d_2 = a_32;

        //vector<double> a1 = cross(a_2,a_3)/sqa;
        cross(a_2, a_3, tmp_f);
        const_division(tmp_f, sqa, a1);
        //vector<double> a2 = cross(a_3,a_1)/sqa;
        cross(a_3, a_1, tmp_f);
        const_division(tmp_f, sqa, a2);
        //vector<double> a11 = 1.0/sqa/sqa *( (cross(a_21,a_3)+cross(a_2,a_31))*sqa - cross(a_2,a_3)*sqa_1 );
        a_cross_b_plus_c_cross_d(a_21, a_3, a_2, a_31, tmp_f, tmp_l, tmp_sum);
        cross(a_2, a_3, tmp_f);
        const_multiplication(tmp_f, sqa_1, tmp_l);
        const_multiplication(tmp_sum, sqa, tmp_f);
        subtraction(tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_sum, 1.0/sqa/sqa, a11);
        //vector<double> a12 = 1.0/sqa/sqa *( (cross(a_22,a_3)+cross(a_2,a_32))*sqa - cross(a_2,a_3)*sqa_2 );
        a_cross_b_plus_c_cross_d(a_22, a_3, a_2, a_32, tmp_f, tmp_l, tmp_sum);
        cross(a_2, a_3, tmp_f);
        const_multiplication(tmp_f, sqa_2, tmp_l);
        const_multiplication(tmp_sum, sqa, tmp_f);
        subtraction(tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_sum, 1.0/sqa/sqa, a12);
        //vector<double> a21 = 1.0/sqa/sqa *( (cross(a_31,a_1)+cross(a_3,a_11))*sqa - cross(a_3,a_1)*sqa_1 );
        a_cross_b_plus_c_cross_d(a_31, a_1, a_3, a_11, tmp_f, tmp_l, tmp_sum);
        cross(a_3, a_1, tmp_f);
        const_multiplication(tmp_f, sqa_1, tmp_l);
        const_multiplication(tmp_sum, sqa, tmp_f);
        subtraction(tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_sum, 1.0/sqa/sqa, a21);
        //vector<double> a22 = 1.0/sqa/sqa *( (cross(a_32,a_1)+cross(a_3,a_12))*sqa - cross(a_3,a_1)*sqa_2 );
        a_cross_b_plus_c_cross_d(a_32, a_1, a_3, a_12, tmp_f, tmp_l, tmp_sum);
        cross(a_3, a_1, tmp_f);
        const_multiplication(tmp_f, sqa_2, tmp_l);
        const_multiplication(tmp_sum, sqa, tmp_f);
        subtraction(tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_sum, 1.0/sqa/sqa, a22);

        H_curv = 0.5 * ( dot(a1,d_1) + dot(a2,d_2) );
        //vector<double> n1_be = -kc*(2.0*H_curv-c0)*(dot(a1,a1)*d_1 + dot(a1,a2)*d_2) + kc*0.5*pow(2.0*H_curv-c0,2.0)*a1;
        tmp_const_f = -kc*(2.0*H_curv-c0);
        tmp_const_l = kc*0.5*pow(2.0*H_curv-c0,2);
        const_multiplication(d_1, dot(a1,a1), tmp_f);
        const_multiplication(d_2, dot(a1,a2), tmp_l);
        addition(tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_sum, tmp_const_f, tmp_f);
        const_multiplication(a1, tmp_const_l, tmp_l);
        addition(tmp_f, tmp_l, n1_be);

        //vector<double> n2_be = -kc*(2.0*H_curv-c0)*(dot(a2,a1)*d_1 + dot(a2,a2)*d_2) + kc*0.5*pow(2.0*H_curv-c0,2.0)*a2;
        const_multiplication(d_1, dot(a2,a1), tmp_f);
        const_multiplication(d_2, dot(a2,a2), tmp_l);
        addition(tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_sum, tmp_const_f, tmp_f);
        const_multiplication(a2, tmp_const_l, tmp_l);
        addition(tmp_f, tmp_l, n2_be);

        //vector<double> m1_be = kc*(2.0*H_curv-c0)*a1;
        const_multiplication(a1, -tmp_const_f, m1_be);

        //vector<double> m2_be = kc*(2.0*H_curv-c0)*a2;
        const_multiplication(a2, -tmp_const_f, m2_be);

        //vector<double> n1_cons = us * (S-S0)*a1;
        //vector<double> n2_cons = us * (S-S0)*a2;
        tmp_const_f = us * (S-S0);
        const_multiplication(a1, tmp_const_f, n1_cons);
        const_multiplication(a2, tmp_const_f, n2_cons);

        //vector<double> n1_conv = 1.0/3.0*uv*(V-V0)*(dot(x,d)*a1 - dot(x,a1)*d);
        //vector<double> n2_conv = 1.0/3.0*uv*(V-V0)*(dot(x,d)*a2 - dot(x,a2)*d);
        tmp_const_f = 1.0/3.0*uv*(V-V0);
        const_multiplication(a1, dot(x,d), tmp_f);
        const_multiplication(d, dot(x,a1), tmp_l);
        subtraction(tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_sum, tmp_const_f, n1_conv);

        const_multiplication(a2, dot(x,d), tmp_f);
        const_multiplication(d, dot(x,a2), tmp_l);
        subtraction(tmp_f, tmp_l, tmp_sum);
        const_multiplication(tmp_sum, tmp_const_f, n2_conv);

        ebe = 0.5 * kc * sqa * pow(2.0*H_curv-c0, 2);    // bending Energy

        //std::cout << "(2H, c0, 2h-c0, ebe)" << 2.0*H_curv << ", " << c0 << ", " << (2.0*H_curv-c0) <<", " << ebe << endl;

        for (int j = 0; j < 12; j++) {
            vector<vector<double>> da1 = -sf[3][j]*kron(a1,d) - sf[1][j]*kron(a11,d) - sf[1][j]*kron(a1,d_1) - sf[6][j]*kron(a2,d) - sf[2][j]*kron(a21,d) - sf[2][j]*kron(a2,d_1);
            vector<vector<double>> da2 = -sf[5][j]*kron(a1,d) - sf[1][j]*kron(a12,d) - sf[1][j]*kron(a1,d_2) - sf[4][j]*kron(a2,d) - sf[2][j]*kron(a22,d) - sf[2][j]*kron(a2,d_2);
            vector<double> tempf_be = n1_be*sf[1][j] + m1_be*da1 + n2_be*sf[2][j] + m2_be*da2;
            f_be[j] = - tempf_be * sqa; // the Force is the negative derivative 
            vector<double> tempf_cons = n1_cons*sf[1][j] + n2_cons*sf[2][j];
            f_cons[j] = - tempf_cons*sqa; // the Force is the negative derivative 
            vector<double> tempf_conv = n1_conv*sf[1][j] + n2_conv*sf[2][j] + 1.0/3.0*uv*(V-V0)*d*sf[0][j];
            f_conv[j] = - tempf_conv*sqa;
        }
        Hmean += 1.0/2.0*GaussQuadratureCoeff[i]*H_curv;
        Ebe   += 1.0/2.0*GaussQuadratureCoeff[i]*ebe;
        F_be  += 1.0/2.0*GaussQuadratureCoeff[i]*f_be;
        F_s   += 1.0/2.0*GaussQuadratureCoeff[i]*f_cons;
        F_v   += 1.0/2.0*GaussQuadratureCoeff[i]*f_conv;
        normVector += 1.0/2.0*GaussQuadratureCoeff[i]*d;
    }
}

void element_energy_force_irregular(vector<vector<double>>& Dots, Param& param, double c0, double& Hmean, vector<double>& normVector, 
double& E_bending, vector<vector<double>>& F_be, vector<vector<double>>& F_s, vector<vector<double>>& F_v, vector<double>& GaussQuadratureCoeff, vector<Shapefunctions>& ShapeFunctions, SubMatrix subMatrix) {
    // F_be is the curvature Force; F_s is the area-constraint Force; F_v is the volume-constraint Force.
    // initialize the output parameters
    E_bending = 0.0;
    // five matrix used for subdivision of the irregular patch
    // M(17,11), M1(12,17), M2(12,17), M3(12,17), M4(11,17); 
    vector<vector<double>> M = subMatrix.irregM; vector<vector<double>> M1 = subMatrix.irregM1; 
    vector<vector<double>> M2 = subMatrix.irregM2; vector<vector<double>> M3 = subMatrix.irregM3; vector<vector<double>> M4 = subMatrix.irregM4;
    /////////////////////////////////////////////////////////////////
    // bending Energy E_bending and cunchu2_constraint Force
    vector<vector<double>> ori_dots = Dots; // (11,3)
    vector<vector<double>> temp (11,vector<double>(11,0.0)); for(int i=0; i < 11; i++) temp[i][i] = 1.0; // identity matrix
    for (int j = 0; j < param.subDivideTimes; j++) {
        vector<vector<double>> newnodes17 = M * ori_dots; // 17 new nodes

        if (j != 0) {
            temp = (M4*M) * temp;
        }
        vector<vector<double>> matrix = M*temp;

        vector<vector<double>> dots = M1*newnodes17;    // element 1
        vector<vector<double>> f_be(12,vector<double>(3,0.0)); // f1(12,3)
        vector<vector<double>> f_s = f_be; 
        vector<vector<double>> f_v = f_be;
        vector<double> normvector(3,0.0); double hmean = 0.0;
        double ebe = 0.0; 
        element_energy_force_regular(dots, param, c0, hmean, normvector, ebe, f_be, f_s, f_v, GaussQuadratureCoeff, ShapeFunctions);
        vector<vector<double>> m1m = transpose(M1*matrix);
        F_be = F_be + m1m*f_be;
        F_s  = F_s  + m1m*f_s;
        F_v  = F_v  + m1m*f_v;
        E_bending = E_bending + ebe;
        Hmean = Hmean + 1.0/pow(4.0,j+1.0) * hmean; // 1/4^(j+1) is the weight-coefficient to calculate the mean value
        normVector = normVector + 1.0/pow(4.0,j+1.0) * normvector;

        dots = M2*newnodes17;    // element 2
        f_be = f_be * 0.0;
        f_s = f_s * 0.0;
        f_v = f_v * 0.0;
        hmean = 0.0; normvector = normvector * 0.0;
        ebe = 0.0;
        element_energy_force_regular(dots, param, c0, hmean, normvector, ebe, f_be, f_s, f_v, GaussQuadratureCoeff, ShapeFunctions);
        vector<vector<double>> m2m = transpose(M2*matrix);
        F_be = F_be + m2m*f_be;
        F_s  = F_s  + m2m*f_s;
        F_v  = F_v  + m2m*f_v;
        E_bending = E_bending + ebe;
        Hmean = Hmean + 1.0/pow(4.0,j+1.0) * hmean; 
        normVector = normVector + 1.0/pow(4.0,j+1.0) * normvector;

        dots = M3*newnodes17;    // element 3
        f_be = f_be * 0.0;
        f_s = f_s * 0.0;
        f_v = f_v * 0.0;
        hmean = 0.0; normvector = normvector * 0.0;
        ebe = 0.0;
        element_energy_force_regular(dots, param, c0, hmean, normvector, ebe, f_be, f_s, f_v, GaussQuadratureCoeff, ShapeFunctions);
        vector<vector<double>> m3m = transpose(M3*matrix);
        F_be = F_be + m3m*f_be;
        F_s  = F_s  + m3m*f_s;
        F_v  = F_v  + m3m*f_v;
        E_bending = E_bending + ebe;
        Hmean = Hmean + 1.0/pow(4.0,j+1.0) * hmean; 
        normVector = normVector + 1.0/pow(4.0,j+1.0) * normvector;

        vector<vector<double>> dots4 = M4 * newnodes17;   // element 4, still irregular patch
        ori_dots = dots4;
    }
}

// rgularization Energy and forces
void energy_force_regularization(vector<Vertex>& vertex, vector<Face>& face, Param& param){
    double k  = param.k;
    double K  = param.K;
    int deformnumber_shape = 0;
    int deformnumber_area = 0;
    vector<vector<double>> fre(vertex.size(), vector<double>(3,0.0));

    #pragma omp parallel for reduction(+:fre)
    for (int i = 0; i < face.size(); i++){    
        bool IsInsertionPatch = face[i].IsInsertionPatch;
        double Ere = 0.0;

        int node0 = face[i].AdjacentVertex[0]; // three nodes of this face element
        int node1 = face[i].AdjacentVertex[1];
        int node2 = face[i].AdjacentVertex[2];
        vector<double> vector0 = vertex[node0].Coord - vertex[node1].Coord;  double side0 = norm(vector0);
        vector<double> vector1 = vertex[node1].Coord - vertex[node2].Coord;  double side1 = norm(vector1);
        vector<double> vector2 = vertex[node2].Coord - vertex[node0].Coord;  double side2 = norm(vector2);
        double s = (side0 + side1 + side2)/2.0;
        double area = sqrt( s*(s-side0)*(s-side1)*(s-side2) );
        double meanside = 1.0/3.0*(side0 + side1+ side2);
        double gama = 1.0/pow(meanside,2.0) * ( pow(side0-meanside,2.0) + pow(side1-meanside,2.0) + pow(side2-meanside,2.0) );
        vector<double> vectorold0 = vertex[node0].ReferenceCoord - vertex[node1].ReferenceCoord;  double sideold0 = norm(vectorold0);
        vector<double> vectorold1 = vertex[node1].ReferenceCoord - vertex[node2].ReferenceCoord;  double sideold1 = norm(vectorold1);
        vector<double> vectorold2 = vertex[node2].ReferenceCoord - vertex[node0].ReferenceCoord;  double sideold2 = norm(vectorold2);
        //std::cout<<"v0"<<vector0.size()<<"::"<<vectorold0.size()<<endl;
        s = (sideold0 + sideold1 + sideold2)/2.0;
        double areaold = sqrt( s*(s-sideold0)*(s-sideold1)*(s-sideold2) ); //double areaold = S0/face.n_rows;
            
        bool isDeformShape = false;
        if ( gama > param.gama_shape && param.usingRpi == true ){
            isDeformShape = true;
        }
        bool isDeformArea = false;
        double a0 = areaold; 
        if ( abs(area-a0)/a0 >= param.gama_area && param.usingRpi == true ){
            isDeformArea = true;
        }
            
        if ( isDeformShape == false && isDeformArea == false ){
            Ere = k/2.0*(pow(side0-sideold0,2.0) + pow(side1-sideold1,2.0) + pow(side2-sideold2,2.0));
            fre[node0] += k*( (side0-sideold0)*(-vector0/side0) + (side2-sideold2)*(vector2/side2) );
            fre[node1] += k*( (side1-sideold1)*(-vector1/side1) + (side0-sideold0)*(vector0/side0) );
            fre[node2] += k*( (side2-sideold2)*(-vector2/side2) + (side1-sideold1)*(vector1/side1) );
        }else if ( isDeformArea == true ){
            deformnumber_area ++;
            double meanside = sqrt( 4.0*areaold/sqrt(3.0) );
            Ere = k/2.0*(pow(side0-meanside,2.0) + pow(side1-meanside,2.0) + pow(side2-meanside,2.0));
            fre[node0] += k*( (side0-meanside)*(-vector0/side0) + (side2-meanside)*(vector2/side2) );
            fre[node1] += k*( (side1-meanside)*(-vector1/side1) + (side0-meanside)*(vector0/side0) );
            fre[node2] += k*( (side2-meanside)*(-vector2/side2) + (side1-meanside)*(vector1/side1) );
        }else if ( isDeformShape == true && isDeformArea == false ){
            deformnumber_shape ++;
            double meansideold = sqrt( 4.0*area/sqrt(3.0) );
            Ere = k/2.0*(pow(side0-meansideold,2.0) + pow(side1-meansideold,2.0) + pow(side2-meansideold,2.0));
            fre[node0] += k*( (side0-meansideold)*(-vector0/side0) + (side2-meansideold)*(vector2/side2) );
            fre[node1] += k*( (side1-meansideold)*(-vector1/side1) + (side0-meansideold)*(vector0/side0) );
            fre[node2] += k*( (side2-meansideold)*(-vector2/side2) + (side1-meansideold)*(vector1/side1) );
        }
        // store Energy
        face[i].energy.energyRegularization = Ere;
    }
    ///////////////////////////////////////////////////
    // store Force
    #pragma omp parallel for
    for (int i = 0; i < vertex.size(); i++){
        vertex[i].force.ForceRegularization = fre[i];
    }
    param.DeformationList.DeformShape     = deformnumber_shape;
    param.DeformationList.DeformArea      = deformnumber_area;
    param.DeformationList.Undeform        = face.size() - deformnumber_shape - deformnumber_area;
}

void Energy_and_Force(vector<Vertex>& vertex, vector<Face>& face, Param& param, vector<double>& GaussQuadratureCoeff, vector<Shapefunctions>& ShapeFunctions, SubMatrix& subMatrix) {
    // update the total area and volume;
    calculate_element_area_volume(vertex, face, param.subDivideTimes, GaussQuadratureCoeff, ShapeFunctions, subMatrix);
    param.S = sum_Membrane_Area(face);    
    param.V = sum_Membrane_Volume(face); // total volume   
    /////////////////////////////////////////////
    // bending Force and Energy, constraint Force and Energy
    clear_forceONvertex_and_energyONface(vertex, face);
    
    vector<Force> force_at_node(vertex.size());
    #pragma omp parallel for reduction(+:force_at_node) 
    for ( int i = 0; i < face.size(); i++) {
        //std::cout << "Face index: " << i << endl;
        if ( face[i].IsBoundary == true )
            continue;
        
        int n = face[i].OneRingVertex.size();
        // one ring vertices
        vector<vector<double>> dots(n,vector<double>(3,0.0)); 
        for (int j = 0; j < n; j++) {
            int nodenum = face[i].OneRingVertex[j];
            dots[j] = vertex[nodenum].Coord;
        }         
        double c0 = face[i].SpontCurvature; // spontaneous curvature of each patch,
        double ebe = 0.0;    // curvature Energy of this element;
        double h_mean = 0.0; // mean curvature of this element;
        vector<double> normVector(3,0.0);
        vector<vector<double>> F_be(n,vector<double>(3,0.0)); // bending or curvature term
        vector<vector<double>> F_consA(n,vector<double>(3,0.0)); // area term
        vector<vector<double>> F_consV(n,vector<double>(3,0.0)); // volume term
        // regular patch
        if ( face[i].OneRingVertex.size() == 12 ){ 
            element_energy_force_regular(dots, param, c0, h_mean, normVector, ebe, F_be, F_consA, F_consV, GaussQuadratureCoeff, ShapeFunctions);
        // irregular patch
        }else if ( face[i].OneRingVertex.size() == 11 ){ 
            element_energy_force_irregular(dots, param, c0, h_mean, normVector, ebe, F_be, F_consA, F_consV, GaussQuadratureCoeff, ShapeFunctions, subMatrix);
        }
        face[i].energy.energyCurvature = ebe;
        face[i].MeanCurvature = h_mean;
        face[i].normVector = normVector;  
        for (int j = 0; j < n; j++) {
            int nodenum = face[i].OneRingVertex[j];
            Force forceTmp;
            forceTmp.ForceCurvature = F_be[j];
            forceTmp.ForceArea      = F_consA[j];
            forceTmp.ForceVolume    = F_consV[j];
            force_at_node[nodenum] += forceTmp;
            //std::cout << "Force at node: " << F_be[j][0] << F_be[j][1]<< F_be[j][2] << F_consA[j][0] << F_consV[j][0] << endl;
        }
    }
    //#pragma omp parallel for 
    for (int i = 0; i < vertex.size(); i++){
        vertex[i].force = force_at_node[i];
    }
    ///////////////////////////////////////////////////////////////////
    // regularization Force and Energy
    energy_force_regularization(vertex, face, param);
    //////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////
    
    // On each vertex, the total Force is the sum of internal, constraint and regularization Force.
    #pragma omp parallel for
    for (int i = 0; i < vertex.size(); i++){
        for (int j = 0; j < 3; j++) {
            if (std::isnan(vertex[i].force.ForceVolume[j])) {
                vertex[i].force.ForceVolume[j] = 0.0;
            }
        }
        vertex[i].force.ForceTotal = vertex[i].force.ForceCurvature + vertex[i].force.ForceArea + vertex[i].force.ForceVolume + vertex[i].force.ForceThickness + vertex[i].force.ForceTilt + vertex[i].force.ForceRegularization;    
    }
    
    // On each face, the total Energy is the sum of internal, constraint and regularization Energy.
    #pragma omp parallel for
    for (int i = 0; i < face.size(); i++){
        face[i].energy.energyTotal = face[i].energy.energyCurvature + face[i].energy.energyArea + face[i].energy.energyVolume + face[i].energy.energyThickness + face[i].energy.energyTilt + face[i].energy.energyRegularization;
        //std::cout<<"face="<<i<<", energy_curv/area/vol/thic/tilt/reg="<<face[i].energy.energyCurvature << face[i].energy.energyArea<< face[i].energy.energyVolume <<face[i].energy.energyThickness <<face[i].energy.energyTilt <<face[i].energy.energyRegularization<<endl;
    }
    //////////////////////////////////////////////////////////////////
    // total Energy to store in param.energy
    //calculate_total_energy_save_in_paramEnergy(face, param);
    Energy tmpEnergy;
    #pragma omp parallel for reduction(+:tmpEnergy)
    for (int i = 0; i < face.size(); i++){
        tmpEnergy += face[i].energy;
    }
    double us = param.us; double S0 = param.S0; double S  = param.S;
    if (param.isGlobalConstraint == true) {
        if (S0 == 0.0) {
            tmpEnergy.energyArea = 0.0;
        } else {
            tmpEnergy.energyArea = 0.5*us/S0*pow(S-S0,2.0);
        }
    }
    double uv = param.uv; double V0 = param.V0; double V  = param.V;
    if (V0 == 0.0) {
        tmpEnergy.energyVolume = 0.0;
    } else {
        tmpEnergy.energyVolume = 0.5*uv/V0*pow(V-V0,2.0);
    }
    tmpEnergy.energyTotal = tmpEnergy.energyCurvature + tmpEnergy.energyArea + tmpEnergy.energyVolume + tmpEnergy.energyThickness + tmpEnergy.energyTilt + tmpEnergy.energyRegularization;
    param.energy = tmpEnergy;

    /////////////////////////////////////////////////////////////////////////
    if (param.isEnergySplineIncluded){
        param.energy.energySpline = calculateSplineEnergyForce(param.splinePoints, vertex, param.splinePoints_correspondingVertexIndex,
                                 param.lbond, param.springConst, false, param.splinePointsZcoordScaling);
                param.energy.energyTotal += param.energy.energySpline;
    }
    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    manage_force_for_boundary_ghost_vertex(vertex, face, param);
    /*
    for (int i = 0; i < vertex.size(); i++) {
        std::cout << "force_curve @ " << i << " = "
                                      << vertex[i].force.ForceCurvature[0] << " , "
                                      << vertex[i].force.ForceCurvature[1] << " , "
                                      << vertex[i].force.ForceCurvature[2] << " , " << endl;
        std::cout << "force_area @ " << i << " = "
                                      << vertex[i].force.ForceArea[0] << " , "
                                      << vertex[i].force.ForceArea[1] << " , "
                                      << vertex[i].force.ForceArea[2] << " , " << endl;
        std::cout << "force_vol @ " << i << " = "
                                      << vertex[i].force.ForceVolume[0] << " , "
                                      << vertex[i].force.ForceVolume[1] << " , "
                                      << vertex[i].force.ForceVolume[2] << " , " << endl;
        std::cout << "force_thic @ " << i << " = "
                                      << vertex[i].force.ForceThickness[0] << " , "
                                      << vertex[i].force.ForceThickness[1] << " , "
                                      << vertex[i].force.ForceThickness[2] << " , " << endl;                              
        std::cout << "force_tilt @ " << i << " = "
                                      << vertex[i].force.ForceTilt[0] << " , "
                                      << vertex[i].force.ForceTilt[1] << " , "
                                      << vertex[i].force.ForceTilt[2] << " , " << endl;
        std::cout << "force_reg @ " << i << " = "
                                      << vertex[i].force.ForceRegularization[0] << " , "
                                      << vertex[i].force.ForceRegularization[1] << " , "
                                      << vertex[i].force.ForceRegularization[2] << " , " << endl;                                                            
        std::cout << "force_total @ " << i << " = "
                                      << vertex[i].force.ForceTotal[0] << " , "
                                      << vertex[i].force.ForceTotal[1] << " , "
                                      << vertex[i].force.ForceTotal[2] << " , " << endl;
        
    }*/
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
double lineSearch_for_StepSize_to_minimize_energy(double a0, vector<vector<double>>& s0, vector<Vertex>& vertex, vector<Face>& face, Param& param, vector<double>& GaussQuadratureCoeff, vector<Shapefunctions>& ShapeFunctions, SubMatrix& subMatrix){
    double a = a0;
    int    numvertex = vertex.size();
    double E0 = param.energy.energyTotal;
    double E1;
    //param.isNCGstucked = false;
    param.usingRpi = true;
    
    double NCGfactor0 = 0.0;
    #pragma omp parallel for reduction(+:NCGfactor0) 
    for (int i = 0; i < vertex.size(); i++){
        NCGfactor0 += dot(-vertex[i].forcePrevious.ForceTotal, s0[i]);
    }
    double c1 = 1e-4; // two factors for the nonlinear conjugate gradient method
    double c2 = 0.1;
    bool   isCriterionSatisfied = false;
    while ( isCriterionSatisfied == false ) {
        a = a * 0.8;        
        update_vertex_from_NonlinearConjugateGradient_s0(vertex, face, a, s0, param); // use CoordPrevious to calculate the new coord
        Energy_and_Force(vertex, face, param, GaussQuadratureCoeff, ShapeFunctions, subMatrix);
   /**     if (param.isEnergySplineIncluded){
            param.energy.energySpline = calculateSplineEnergyForce(param.splinePoints, vertex, param.splinePoints_correspondingVertexIndex,
                                param.lbond, param.springConst, false, param.splinePointsZcoordScaling);
            param.energy.energyTotal += param.energy.energySpline;
        }*/
        E1 = param.energy.energyTotal;
        
        if ( param.usingNCG == true && param.isNCGstucked == false ){
            double NCGfactor = 0.0;
            #pragma omp parallel for reduction(+:NCGfactor) 
            for (int i = 0; i < vertex.size(); i++){
                NCGfactor += dot(-vertex[i].force.ForceTotal, s0[i]);
            }
            //if ( E1 <= E0 + c1*a*NCGfactor0 && NCGfactor >= c2*NCGfactor0 ){ // Wolfe conditions
            if ( E1 <= E0 + c1*a*NCGfactor0 && abs(NCGfactor) <= c2*abs(NCGfactor0) ){ // strong Wolfe conditions
                break;
            }
            if ( a < 1.0e-15 ) {
                cout<<"Now change the NCG WolfeConditions to simple line search method!"<<endl;
                //WolfeConditions = false;
                a = a0;
                update_reference_from_CoordPrevious(vertex);
                param.isNCGstucked = true;
                param.usingRpi = false;
                // change s0 to nodalForce
                #pragma omp parallel for 
                for (int i = 0; i < vertex.size(); i++){
                    s0[i] = vertex[i].forcePrevious.ForceTotal;
                }
            }
        }else{
            param.isNCGstucked = true;
            if ( E1 < E0 ) {
                //cout<<"Simple method, a = "<<a<<", Ebe = "<<energy1(0)<<", Econs = "<<energy1(1)<<", Ebar = "<<energy1(2)<<", Ere = "<<energy1(3)<<", Etot = "<<energy1(5)<<endl;
                break;
            }
            if ( a < 1.0e-15 ) {
                 cout<<"Note: cannot find an efficient samll a to minimize the Energy even with the simple method!"<<endl;
                 return -1;
            }
        }
        
    }
    return a;
}

/*
void printoutGlobalNCG(mat vertex1) {
    int numvertex = vertex1.n_rows;
    ofstream outfile("vertex_Global_NCG.csv");
    for (int j = 0; j < numvertex; j++) {
        outfile << vertex1(j,0)+0.0 << ',' << vertex1(j,1)+0.0 << ',' << vertex1(j,2)+0.0 << '\n';
    }
    outfile.close();
}
void printoutGlobalGD(mat vertex1) {
    int numvertex = vertex1.n_rows;
    ofstream outfile("vertex_Global_GD.csv");
    for (int j = 0; j < numvertex; j++) {
        outfile << vertex1(j,0)+0.0 << ',' << vertex1(j,1)+0.0 << ',' << vertex1(j,2)+0.0 << '\n';
    }
    outfile.close();
}
void printoutstuck(mat vertex1) {
    int numvertex = vertex1.n_rows;
    ofstream outfile("vertex_stuck.csv");
    for (int j = 0; j < numvertex; j++) {
        outfile << vertex1(j,0)+0.0 << ',' << vertex1(j,1)+0.0 << ',' << vertex1(j,2)+0.0 << '\n';
    }
    outfile.close();
}

void printoutREF(mat vertex1){
    int numvertex = vertex1.n_rows;
    ofstream outfile("vertex_reference.csv");
    for (int j = 0; j < numvertex; j++) {
        outfile << vertex1(j,0)+0.0 << ',' << vertex1(j,1)+0.0 << ',' << vertex1(j,2)+0.0 << '\n';
    }
    outfile.close();
}

void printoutforce(mat vertex1){
    int numvertex = vertex1.n_rows;
    ofstream outfile("Force.csv");
    for (int j = 0; j < numvertex; j++) {
        outfile << vertex1(j,0)+0.0 << ',' << vertex1(j,1)+0.0 << ',' << vertex1(j,2)+0.0 << '\n';
    }
    outfile.close();
}

void printout_spontaneouscurvature(rowvec spontcurv){
    ofstream outfile11("spont.csv");
    for (int i = 0; i < spontcurv.n_cols; i++) {
        outfile11 << setprecision(16) << spontcurv(i)+0.0 << '\n';
    }
    outfile11.close(); 
}
*/

void check_nodal_force(vector<Vertex>& vertex, vector<Face>& face, Param& param, vector<double>& GaussQuadratureCoeff, vector<Shapefunctions>& ShapeFunctions, SubMatrix& subMatrix){
    cout<<"check if the nodal Force is correct: "<<endl;
    Energy_and_Force(vertex, face, param, GaussQuadratureCoeff, ShapeFunctions, subMatrix);
    double E = param.energy.energyTotal;
    int numvertex = vertex.size();
    vector<vector<double>> ForceReal(numvertex,vector<double>(3,0.0)); 
    double dx = 1.0e-8;
    for (int i = 0; i < numvertex; i++) {
        if (vertex[i].IsGhost == true) 
            continue;

        vector<Vertex> vertextry = vertex;
        vertextry[i].Coord[0] = vertex[i].Coord[0] + dx;
        Energy_and_Force(vertextry, face, param, GaussQuadratureCoeff, ShapeFunctions, subMatrix);
        double Etmpx = param.energy.energyTotal;
        double fx = - (Etmpx - E)/dx;

        vertextry = vertex;
        vertextry[i].Coord[1] = vertex[i].Coord[1] + dx;
        Energy_and_Force(vertextry, face, param, GaussQuadratureCoeff, ShapeFunctions, subMatrix);
        double Etmpy = param.energy.energyTotal;
        double fy = - (Etmpy - E)/dx;

        vertextry = vertex;
        vertextry[i].Coord[2] = vertex[i].Coord[2] + dx;
        Energy_and_Force(vertextry, face, param, GaussQuadratureCoeff, ShapeFunctions, subMatrix);
        double Etmpz = param.energy.energyTotal;
        double fz = - (Etmpz - E)/dx;

        ForceReal[i][0] = fx;
        ForceReal[i][1] = fy;
        ForceReal[i][2] = fz;
        double Freal = norm(ForceReal[i]);
        double F = norm(vertex[i].force.ForceTotal);
        double dF = norm(ForceReal[i] - vertex[i].force.ForceTotal);
        cout<<i<<std::setprecision(15)<<". Freal = "<<Freal<<", F = "<<F<<", dF = "<<dF<<endl;
    }
}
