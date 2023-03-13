//
// Created by hughm on 13/03/2023.
//

#ifndef ASSIGNMENT1_VECTOR_UTILS_H
#define ASSIGNMENT1_VECTOR_UTILS_H

//
// Created by hughm on 13/03/2023.
//

#include <iostream>
#include <vector>
#include <cassert>
#include <functional>

#include "bspline.hpp"

using function_1D  = std::function<double(double)>;
using list_of_vecs = std::vector<std::vector<double>>;

//=======================================================
//Utility Functions
void printv(const std::vector<double>& a);

//Overload vector operations to make direct products easier
//V-V Multiplication
std::vector<double> operator*=(std::vector<double> & a, const std::vector<double> & b);
std::vector<double> operator*(std::vector<double> a, const std::vector<double> & b);

//V-V Addition
std::vector<double> operator+=(std::vector<double> & a, const std::vector<double> & b);
std::vector<double> operator+(std::vector<double> a, const std::vector<double> & b);

//V-D Multiplication
std::vector<double> operator*(double & a, std::vector<double> v);
std::vector<double> operator*(std::vector<double> v, double & a){return(a*v);};

//Conversions & Integrations
double vint(const std::vector<double>& a);
std::vector<double> make_rgrid(double rmin = 0.001, double rmax = 100, int n_grid = 101);
std::vector<double> potvec_from_func(const function_1D& V, const std::vector<double>& rgrid);

//=======================================================
//B-Spline Evaluation
list_of_vecs generate_splines(int nsplines, const std::vector<double> & rgrid, int k_spline=7);
list_of_vecs generate_spline_diffs(int nsplines, const std::vector<double> & rgrid, int k_spline=7);

#endif //ASSIGNMENT1_VECTOR_UTILS_H
