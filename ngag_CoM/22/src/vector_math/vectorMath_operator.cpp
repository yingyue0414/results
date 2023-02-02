
// vector math
#include "vector_math.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////
vector<double> operator+(const vector<double>& v1, const vector<double>& v2){
    if ( v1.size() != v2.size() ){
        cout<<"Wrong! Can't add two vectors that are not the same dimention!"<<endl;
        exit(0);
    }
    vector<double> tmp(v1.size());
    for (int i = 0; i < v1.size(); i++){
        tmp[i] = v1[i] + v2[i];
    }
    return tmp;
}
void operator+=(vector<double>& v1, const vector<double>& v2){
    if ( v1.size() != v2.size() ){
        cout<<"Wrong! Can't add two vectors that are not the same dimention!"<<endl;
        exit(0);
    }
    for (int i = 0; i < v1.size(); i++){
        v1[i] += v2[i];
    }
}
vector<double> operator-(const vector<double>& v1, const vector<double>& v2){
    if ( v1.size() != v2.size() ){
        cout<<"Wrong! Can't minus two vectors that are not the same dimention!"<<endl;
        exit(0);
    }
    vector<double> tmp(v1.size());
    for (int i = 0; i < v1.size(); i++){
        tmp[i] = v1[i] - v2[i];
    }
    return tmp;
}
void operator-=(vector<double>& v1, const vector<double>& v2){
    if ( v1.size() != v2.size() ){
        cout<<"Wrong! Can't plus two vectors that are not the same dimention!"<<endl;
        exit(0);
    }
    for (int i = 0; i < v1.size(); i++){
        v1[i] -= v2[i];
    }
}
vector<double> operator-(const vector<double>& v1){
    vector<double> tmp(v1.size());
    for (int i = 0; i < v1.size(); i++){
        tmp[i] = - v1[i];
    }
    return tmp;
}
vector<double> operator*(const vector<double>& v1, const double num){
    vector<double> tmp(v1.size());
    for (int i = 0; i < v1.size(); i++){
        tmp[i] = v1[i] * num;
    }
    return tmp;
}
vector<double> operator*(const double num, const vector<double>& v1){
    vector<double> tmp(v1.size());
    for (int i = 0; i < v1.size(); i++){
        tmp[i] = v1[i] * num;
    }
    return tmp;
}
vector<double> operator/(const vector<double>& v1, const double num){
    vector<double> tmp(v1.size());
    for (int i = 0; i < v1.size(); i++){
        tmp[i] = v1[i] / num;
    }
    return tmp;
}
double dot(const vector<double>& v1, const vector<double>& v2){
    if ( v1.size() != v2.size() ){
        cout<<"Wrong! Can't dot_product two vectors that are not the same dimention!"<<endl;
        exit(0);
    }
    double tmp = 0.0;;
    for (int i = 0; i < v1.size(); i++){
        tmp += v1[i] * v2[i];
    }
    return tmp;
}
vector<double> cross(const vector<double>& m1, const vector<double>& m2){ // only works for 3D
    vector<double> temp(3);
    temp[0] = m1[1] * m2[2] - m1[2] * m2[1];
    temp[1] = m1[2] * m2[0] - m1[0] * m2[2];
    temp[2] = m1[0] * m2[1] - m1[1] * m2[0];
    return temp;
}
double norm(const vector<double>& v1){
    double sum = 0.0;
    for (int i = 0; i < v1.size(); i++){
        sum += v1[i] * v1[i];
    }
    return sqrt(sum);
}

vector<vector<double>> operator+(const vector<vector<double>>& m1, const vector<vector<double>>& m2){
    if ( m1.size() != m2.size() || m1[0].size() != m2[0].size() ){
        cout<<"Wrong! Can't plus two matrix that are not the same dimention!"<<endl;
        exit(0);
    }
    vector<vector<double>> tmp (m1.size(), vector<double>(m1[0].size()));
    for (int i = 0; i < m1.size(); i++){
        for (int j = 0; j < m1[0].size(); j++){
            tmp[i][j] = m1[i][j] + m2[i][j];
        }
    }
    return tmp;
}

void operator+=(vector<vector<double>>& m1, const vector<vector<double>>& m2){
    if ( m1.size() != m2.size() || m1[0].size() != m2[0].size() ){
        cout<<"Wrong! Can't plus two matrix that are not the same dimention!"<<endl;
        exit(0);
    }
    for (int i = 0; i < m1.size(); i++){
        m1[i] += m2[i];
    }
}

vector<vector<double>> operator-(const vector<vector<double>>& m1, const vector<vector<double>>& m2){
    if ( m1.size() != m2.size() || m1[0].size() != m2[0].size() ){
        cout<<"Wrong! Can't minus two matrix that are not the same dimention!"<<endl;
        exit(0);
    }
    vector<vector<double>> tmp (m1.size(), vector<double>(m1[0].size()));
    for (int i = 0; i < m1.size(); i++){
        for (int j = 0; j < m1[0].size(); j++){
            tmp[i][j] = m1[i][j] - m2[i][j];
        }
    }
    return tmp;
}
vector<vector<double>> operator*(const vector<vector<double>>& m1, const vector<vector<double>>& m2){
    if ( m1[0].size() != m2.size() ){
        cout<<"Wrong! Wrong dimentions when to time two matrix."<<endl;
        exit(0);
    }
    vector<vector<double>> tmp (m1.size(), vector<double>(m2[0].size()));
    for (int i = 0; i < m1.size(); i++){
        for (int j = 0; j < m2[0].size(); j++){
            double sum = 0.0;
            for (int k = 0; k < m2.size(); k++){
                sum += m1[i][k] * m2[k][j];
            }
            tmp[i][j] = sum;
        }
    }
    return tmp;
}
vector<vector<double>> operator*(const vector<vector<double>>& m1, const double num){
    vector<vector<double>> tmp (m1.size(), vector<double>(m1[0].size()));
    for (int i = 0; i < m1.size(); i++){
        tmp[i] = m1[i] * num;
    }
    return tmp;
}
vector<vector<double>> operator*(const double num, const vector<vector<double>>& m1){
    vector<vector<double>> tmp (m1.size(), vector<double>(m1[0].size()));
    for (int i = 0; i < m1.size(); i++){
        tmp[i] = m1[i] * num;
    }
    return tmp;
}
vector<vector<double>> operator/(const vector<vector<double>>& m1, const double num){
    vector<vector<double>> tmp (m1.size(), vector<double>(m1[0].size()));
    for (int i = 0; i < m1.size(); i++){
        tmp[i] = m1[i] / num;
    }
    return tmp;
}
vector<double> operator*(const vector<double>& v1, const vector<vector<double>>& m1){
    if ( v1.size() != m1.size() ){
        cout<<"Wrong! Wrong dimentions when to time a vector and matrix."<<endl;
        exit(0);
    }
    vector<double> tmp (m1[0].size());
    for (int i = 0; i < m1[0].size(); i++){
        double sum = 0.0;
        for (int j = 0; j < v1.size(); j++){
            sum += v1[j] * m1[j][i];
        }
        tmp[i] = sum;
    }
    return tmp;
}
double norm(const vector<vector<double>>& m1){
    if ( m1[0].size() > 1 ){
        cout<<"Can't calculate the norm for a matrix!"<<endl;
        exit(0);
    }
    double sum = 0.0;
    for (int i = 0; i < m1.size(); i++){
        sum += m1[i][0] * m1[i][0];
    }
    return sqrt(sum);
}
vector<vector<double>> transpose(const vector<double>& v1){
    vector<vector<double>> tmp(v1.size(),vector<double>(1));
    for (int i = 0; i < v1.size(); i++){
        tmp[i][0] = v1[i];
    }
    return tmp;
}
vector<vector<double>> transpose(const vector<vector<double>>& m1){
    vector<vector<double>> tmp(m1[0].size(),vector<double>(m1.size()));
    for (int i = 0; i < m1[0].size(); i++){
        for (int j = 0; j < m1.size(); j++){
            tmp[i][j]= m1[j][i];
        }
    }
    return tmp;
}
vector<vector<double>> kron(const vector<double>& v1, const vector<double>& v2){ // kron( v1', v2 ); note v1 needs to be transposed first.
    vector<vector<double>> tmp(v1.size(),vector<double>(v2.size()));
    for (int i = 0; i < v1.size(); i++){
        tmp[i] = v1[i] * v2;
    }
    return tmp;
}

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
// struct force
Force operator+(const Force& F1, const Force& F2){
    Force Ftmp;
    Ftmp.ForceCurvature = F1.ForceCurvature + F2.ForceCurvature; 
    Ftmp.ForceArea      = F1.ForceArea      + F2.ForceArea; 
    Ftmp.ForceVolume    = F1.ForceVolume    + F2.ForceVolume; 
    Ftmp.ForceThickness = F1.ForceThickness + F2.ForceThickness; 
    Ftmp.ForceTilt   = F1.ForceTilt      + F2.ForceTilt; 
    Ftmp.ForceRegularization = F1.ForceRegularization + F2.ForceRegularization; 
    Ftmp.ForceTotal     = F1.ForceTotal     + F2.ForceTotal; 
    return Ftmp;
}
void operator+=(Force& F1, const Force& F2){
    F1.ForceCurvature += F2.ForceCurvature; 
    F1.ForceArea      += F2.ForceArea; 
    F1.ForceVolume    += F2.ForceVolume; 
    F1.ForceThickness += F2.ForceThickness; 
    F1.ForceTilt      += F2.ForceTilt; 
    F1.ForceRegularization += F2.ForceRegularization; 
    F1.ForceTotal     += F2.ForceTotal; 
}
vector<Force> operator+(const vector<Force>& F1, const vector<Force>& F2){
    if ( F1.size() != F2.size() ){
        cout<<"Wrong! Can't plus two vectors that are not the same dimention!"<<endl;
        exit(0);
    }
    vector<Force> Ftmp(F1.size());
    for (int i = 0; i < F1.size(); i++){
        Ftmp[i] = F1[i] + F2[i];
    }
    return Ftmp;
}
void operator+=(vector<Force>& F1, const vector<Force>& F2){
    if ( F1.size() != F2.size() ){
        cout<<"Wrong! Can't plus two vectors that are not the same dimention!"<<endl;
        exit(0);
    }
    for (int i = 0; i < F1.size(); i++){
        F1[i] += F2[i];
    }
}
////////////////////////////////////////////
// struct energy
Energy operator+(const Energy& E1, const Energy& E2){
    Energy Etmp;
    Etmp.energyCurvature = E1.energyCurvature + E2.energyCurvature;
    Etmp.energyArea      = E1.energyArea      + E2.energyArea;
    Etmp.energyVolume    = E1.energyVolume    + E2.energyVolume;
    Etmp.energyThickness = E1.energyThickness + E2.energyThickness;
    Etmp.energyTilt      = E1.energyTilt      + E2.energyTilt;
    Etmp.energyRegularization = E1.energyRegularization + E2.energyRegularization;
    Etmp.energyTotal     = E1.energyTotal     + E2.energyTotal;
    return Etmp;
}
void operator+=(Energy& E1, const Energy& E2){
    E1.energyCurvature += E2.energyCurvature;
    E1.energyArea      += E2.energyArea;
    E1.energyVolume    += E2.energyVolume;
    E1.energyThickness += E2.energyThickness;
    E1.energyTilt      += E2.energyTilt;
    E1.energyRegularization += E2.energyRegularization;
    E1.energyTotal     += E2.energyTotal;
}
vector<Energy> operator+(const vector<Energy>& E1, const vector<Energy>& E2){
    if ( E1.size() != E2.size() ){
        cout<<"Wrong! Can't plus two vectors that are not the same dimention!"<<endl;
        exit(0);
    }
    vector<Energy> Etmp(E1.size());
    for (int i = 0; i < E1.size(); i++){
        Etmp[i] = E1[i] + E2[i];
    }
    return Etmp;
}
void operator+=(vector<Energy>& E1, const vector<Energy>& E2){
    if ( E1.size() != E2.size() ){
        cout<<"Wrong! Can't plus two vectors that are not the same dimention!"<<endl;
        exit(0);
    }
    for (int i = 0; i < E1.size(); i++){
        E1[i] += E2[i];
    }
}

//////////////////////////////////////////
// solve the matrix 
/*
vector<double> solve_AtimeXeqb(vector<vector<double>>& A, vector<vector<double>>& b){
    if (A.size() != A[0].size() || A.size() != b.size() || b[0].size() != 1){ // A is not a square matrix, or A,b are not consistant in dimension
        cout<<"Wrong! Can't solve this equations, because of wrong dimentions!"<<endl;
        exit(0);
    }
    // augment matrix.
    vector<vector<double>> Ab(A.size(), A[0].size() + 1);
    for (int i = 0; i < Ab.size(); i++){
        for (int j = 0; j < Ab[0].size(); j++){
            if (j < A[0].size()){
                Ab[i][j] = A[i][j];
            }else{
                Ab[i][j] = b[i][0];
            }
        }
    }
    // identity matrix on the left side
    for (int i = 0; i < Ab.size(); i++){
        Ab[i] = Ab[i] / Ab[i][i];
        for (int j = i+1; j < Ab.size(); j++){
            Ab[j] = Ab[j] - 
        }
    }
}
*/



/////////////////////////////////////////////////////////////////////////////////
// define gauss randoms
double gr() {
    static double v1, v2, s;
    double x;
    if (phase == 0) {
        do {
            double u1 = (double)rand() / RAND_MAX;
            double u2 = (double)rand() / RAND_MAX;
            v1 = 2 * u1 - 1;
            v2 = 2 * u2 - 1;
            s = v1 * v1 + v2 * v2;
        } while (s >= 1 || s == 0);
        x = v1 * sqrt(-2 * log(s) / s);
    } else {
        x = v2 * sqrt(-2 * log(s) / s);
    }
    phase = 1 - phase;
    return x;
}
