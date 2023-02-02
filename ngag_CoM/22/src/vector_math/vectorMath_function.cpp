#include "vector_math.hpp"

using namespace std;

/*
 * This cpp defines the void output version of vectorMath operators which saves the 
 * output to a result pointer input instead of allocating new memory space for the
 * result. Using these methods would be recommended for max efficiency, but do note
 * that you need to define the CORRECT dimensions for you result (in comparison you
 * do not have to when using operators defined in vectorMath_operator.cpp).
 * 
 * - Y Ying (Dec 6, 2022)
 */

void addition(const vector<double>& v1, const vector<double>& v2, vector<double>& v_sum) {
    if ( v1.size() != v2.size() ){
        cout<<"Wrong! Can't add two vectors that are not of the same dimention!"<<endl;
        exit(0);
    }
    for (int i = 0; i < v1.size(); i++){
        v_sum[i] = v1[i] + v2[i];
    }
}

void subtraction(const vector<double>& v1, const vector<double>& v2, vector<double>& v_diff){
    if ( v1.size() != v2.size() ){
        cout<<"Wrong! Can't subtract two vectors that are not of the same dimention!"<<endl;
        exit(0);
    }
    for (int i = 0; i < v1.size(); i++){
        v_diff[i] = v1[i] - v2[i];
    }
}
void negative(const vector<double>& v1, vector<double>& v_neg){
    for (int i = 0; i < v1.size(); i++){
        v_neg[i] = - v1[i];
    }
}
void const_multiplication(const vector<double>& v1, const double num, vector<double>& v_result){
    for (int i = 0; i < v1.size(); i++){
        v_result[i] = v1[i] * num;
    }
}
void const_division(const vector<double>& v1, const double num, vector<double>& tmp){
    for (int i = 0; i < v1.size(); i++){
        tmp[i] = v1[i] / num;
    }
}
void cross(const vector<double>& m1, const vector<double>& m2, vector<double>& tmp){ // only works for 3D
    tmp[0] = m1[1] * m2[2] - m1[2] * m2[1];
    tmp[1] = m1[2] * m2[0] - m1[0] * m2[2];
    tmp[2] = m1[0] * m2[1] - m1[1] * m2[0];
}

void addition(const vector<vector<double>>& m1, const vector<vector<double>>& m2, vector<vector<double>>& tmp){
    if ( m1.size() != m2.size() || m1[0].size() != m2[0].size() ){
        cout<<"Wrong! Can't add two matrix that are not the same dimention!"<<endl;
        exit(0);
    }
    for (int i = 0; i < m1.size(); i++){
        for (int j = 0; j < m1[0].size(); j++){
            tmp[i][j] = m1[i][j] + m2[i][j];
        }
    }
}
void subtraction(const vector<vector<double>>& m1, const vector<vector<double>>& m2, vector<vector<double>>& tmp){
    if ( m1.size() != m2.size() || m1[0].size() != m2[0].size() ){
        cout<<"Wrong! Can't subtract two matrix that are not the same dimention!"<<endl;
        exit(0);
    }
    for (int i = 0; i < m1.size(); i++){
        for (int j = 0; j < m1[0].size(); j++){
            tmp[i][j] = m1[i][j] - m2[i][j];
        }
    }
}
void multiplication(const vector<vector<double>>& m1, const vector<vector<double>>& m2, vector<vector<double>>& tmp){
    if ( m1[0].size() != m2.size() ){
        cout<<"Wrong! Wrong dimentions when multiplying two matrix."<<endl;
        exit(0);
    }
    double sum = 0.0;
    for (int i = 0; i < m1.size(); i++){
        for (int j = 0; j < m2[0].size(); j++){
            sum = 0.0;
            for (int k = 0; k < m2.size(); k++){
                sum += m1[i][k] * m2[k][j];
            }
            tmp[i][j] = sum;
        }
    }
}
void getUnitVector(const vector<double>& v1, vector<double>& v_unit){
    double distance2 = 0.0;
    for (int i = 0; i < v1.size(); i++){
        distance2 += v1[i] * v1[i];
    }
    const_division(v1, distance2, v_unit);
}
void const_multiplication(const vector<vector<double>>& m1, const double num, vector<vector<double>>& tmp){
    for (int i = 0; i < m1.size(); i++){
        tmp[i] = m1[i] * num;
    }
}
void const_division(const vector<vector<double>>& m1, const double num, vector<vector<double>>& tmp){
    for (int i = 0; i < m1.size(); i++){
        tmp[i] = m1[i] / num;
    }
}
void rowvec_matrix_multiplication(const vector<double>& v1, const vector<vector<double>>& m1, vector<double>& tmp){
    if ( v1.size() != m1.size() ){
        cout<<"Wrong! Wrong dimentions when multiplying a vector by a matrix."<<endl;
        exit(0);
    }
    double sum = 0.0;
    for (int i = 0; i < m1[0].size(); i++){
        sum = 0.0;
        for (int j = 0; j < v1.size(); j++){
            sum += v1[j] * m1[j][i];
        }
        tmp[i] = sum;
    }
}

void kron (const vector<double>& v1, const vector<double>& v2, vector<vector<double>>& m_result) {
    for (int i = 0; i < v1.size(); i++){
        for (int j = 0; j < v2.size(); j++) {
            m_result[i][j] = v1[i] * v2[j];
        }
    }
}

//Complex method

//calculate a x b + c x d assuming a,b,c,d are all vector(3)
//using 3 tmp vectors (vector(3))
void a_cross_b_plus_c_cross_d(  const vector<double>& a,
                                const vector<double>& b,
                                const vector<double>& c,
                                const vector<double>& d,
                                vector<double>& tmp_f,
                                vector<double>& tmp_l,
                                vector<double>& v_result) {
    
    cross(a, b, tmp_f);
    cross(c, d, tmp_l);
    addition(tmp_f, tmp_l, v_result);
}

//bring down the second dimension of shapefunction matrix
vector<vector<double>> down_2nd_dimension_sf(const vector<Shapefunctions>& sf, int j) {
    const int k_max = sf[0].sf[0].size();
    vector<vector<double>> sf_mat (sf.size(), vector<double>(k_max));
    for (int i = 0; i < sf.size(); i++) {
        for (int k = 0; k < k_max; k++) {
            sf_mat[i][k] = sf[i].sf[j][k];
        }
    }
    return sf_mat;
}