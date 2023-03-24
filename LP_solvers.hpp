//
// Created by hughm on 12/03/2023.
//

#include "matrix.hpp"

#ifndef ASSIGNMENT1_LP_SOLVERS_H
#define ASSIGNMENT1_LP_SOLVERS_H

namespace lapack{
    using namespace matrix;
    struct MatrixAndVector {
        //Init Function
        MatrixAndVector(unt dimension): matrix(dimension), vector(dimension) {}

        sqmatrix matrix;
        std::vector<double> vector;
    };

    MatrixAndVector LP_Eig_A(sqmatrix A);

    MatrixAndVector LP_Eig_AB(sqmatrix A, sqmatrix B);

}

#endif //ASSIGNMENT1_LP_SOLVERS_H
