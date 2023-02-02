#pragma once

#include <iostream>
#include <vector>
#include <math.h>
#include <parameters.hpp>


//operator
////////////////////////////////////////////////////
// vector calculation
std::vector<double> operator+(const std::vector<double>& v1, const std::vector<double>& v2);
void operator+=(std::vector<double>& v1, const std::vector<double>& v2);
std::vector<double> operator-(const std::vector<double>& v1, const std::vector<double>& v2);
void operator-=(std::vector<double>& v1, const std::vector<double>& v2);
std::vector<double> operator-(const std::vector<double>& v1);
std::vector<double> operator*(const std::vector<double>& v1, const double num);
std::vector<double> operator*(const double num, const std::vector<double>& v1);
std::vector<double> operator/(const std::vector<double>& v1, const double num);
double dot(const std::vector<double>& v1, const std::vector<double>& v2);
std::vector<double> cross(const std::vector<double>& m1, const std::vector<double>& m2);
double norm(const std::vector<double>& v1);

std::vector<std::vector<double>> operator+(const std::vector<std::vector<double>>& m1, const std::vector<std::vector<double>>& m2);

void operator+=(std::vector<std::vector<double>>& m1, const std::vector<std::vector<double>>& m2);

std::vector<std::vector<double>> operator-(const std::vector<std::vector<double>>& m1, const std::vector<std::vector<double>>& m2);
std::vector<std::vector<double>> operator*(const std::vector<std::vector<double>>& m1, const std::vector<std::vector<double>>& m2);
std::vector<std::vector<double>> operator*(const std::vector<std::vector<double>>& m1, const double num);
std::vector<std::vector<double>> operator*(const double num, const std::vector<std::vector<double>>& m1);
std::vector<std::vector<double>> operator/(const std::vector<std::vector<double>>& m1, const double num);
std::vector<double> operator*(const std::vector<double>& v1, const std::vector<std::vector<double>>& m1);
double norm(const std::vector<std::vector<double>>& m1);
std::vector<std::vector<double>> transpose(const std::vector<double>& v1);
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& m1);
std::vector<std::vector<double>> kron(const std::vector<double>& v1, const std::vector<double>& v2);

/////////////////////////////////////////////////////
/////////////////////////////////////////////////////
// struct force
Force operator+(const Force& F1, const Force& F2);
void operator+=(Force& F1, const Force& F2);
std::vector<Force> operator+(const std::vector<Force>& F1, const std::vector<Force>& F2);
void operator+=(std::vector<Force>& F1, const std::vector<Force>& F2);
////////////////////////////////////////////
// struct energy
Energy operator+(const Energy& E1, const Energy& E2);
void operator+=(Energy& E1, const Energy& E2);
std::vector<Energy> operator+(const std::vector<Energy>& E1, const std::vector<Energy>& E2);
void operator+=(std::vector<Energy>& E1, const std::vector<Energy>& E2);
//////////////////////////////////////////
// solve the matrix 



/////////////////////////////////////////////////////////////////////////////////
// define gauss randoms
static int phase = 0;
double gr();


//function
void addition(const std::vector<double>& v1, const std::vector<double>& v2, std::vector<double>& v_sum);
void subtraction(const std::vector<double>& v1, const std::vector<double>& v2, std::vector<double>& v_diff);
void negative(const std::vector<double>& v1, std::vector<double>& v_neg);
void getUnitVector(const std::vector<double>& v1, std::vector<double>& v_unit);
void const_multiplication(const std::vector<double>& v1, const double num, std::vector<double>& v_result);
void const_division(const std::vector<double>& v1, const double num, std::vector<double>& tmp);
void cross(const std::vector<double>& m1, const std::vector<double>& m2, std::vector<double>& temp);
void addition(const std::vector<std::vector<double>>& m1, const std::vector<std::vector<double>>& m2, std::vector<std::vector<double>>& tmp);
void subtraction(const std::vector<std::vector<double>>& m1, const std::vector<std::vector<double>>& m2, std::vector<std::vector<double>>& tmp);
void multiplication(const std::vector<std::vector<double>>& m1, const std::vector<std::vector<double>>& m2, std::vector<std::vector<double>>& tmp);
void const_multiplication(const std::vector<std::vector<double>>& m1, const double num, std::vector<std::vector<double>>& tmp);
void const_division(const std::vector<std::vector<double>>& m1, const double num, std::vector<std::vector<double>>& tmp);
void rowvec_matrix_multiplication(const std::vector<double>& v1, const std::vector<std::vector<double>>& m1, std::vector<double>& tmp);
void kron (const std::vector<double>& v1, const std::vector<double>& v2, std::vector<std::vector<double>>& m_result);

//complex functions
void a_cross_b_plus_c_cross_d(  const std::vector<double>& a,
                                const std::vector<double>& b,
                                const std::vector<double>& c,
                                const std::vector<double>& d,
                                std::vector<double>& tmp_f,
                                std::vector<double>& tmp_l,
                                std::vector<double>& v_result);
std::vector<std::vector<double>> down_2nd_dimension_sf(const std::vector<Shapefunctions>& sf, int j);