//
// Created by hughm on 18/03/2023.
//

#ifndef ASSIGNMENT1_SPLINE_EVAL_H
#define ASSIGNMENT1_SPLINE_EVAL_H

#include <vector>

using list_of_vecs = std::vector<std::vector<double>>;
list_of_vecs generate_splines(int nsplines, const std::vector<double> & rgrid, int k_spline = 7);
list_of_vecs generate_spline_diffs(int nsplines, const std::vector<double> & rgrid, int k_spline = 7, double dr=-1);

#endif //ASSIGNMENT1_SPLINE_EVAL_H
