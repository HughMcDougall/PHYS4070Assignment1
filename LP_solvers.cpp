//
// Wrapper functions for lapack routines
// HM Mar 23
//

#include <cassert>
#include <iostream>
#include "LP_solvers.hpp"
#include "matrix.hpp"

namespace lapack{
    using namespace matrix;

    //---------------------------------
    //Pull in LAPACK functions
    //---------------------------------
    extern "C"
    int dsyev_(
            char *jobz, char *uplo, int *dimension, double *A,
           int *dimension2, double *out_eigen_values, double *workspace_array,
           int *workspace_size, int *info
           );

    extern "C"
    int dsygv_(
            int *itype,
            char *jobz,
            char *uplo,
            int *dimension,
            double *A,
            int *dimensionA,
            double *B,
            int *dimensionB,
            double *out_eigen_values,
            double *workspace_array,
            int *workspace_size,
            int *info
        );

    MatrixAndVector LP_Eig_A(sqmatrix A){
        /// Solves the eigenvector system [Ax=kx] for x and k
        /// Returns as struct "MatrixAndVector out" with solutions
        /// Stored in out.vector and out.matrix

        //Creat dummy storage for outputs
        MatrixAndVector result(A.N());

        //Set flags
        char jobz{'V'};
        char uplo{'U'};
        int dimension = static_cast<int>(A.N());
        int lwork = 6 * dimension;
        std::vector<double> work(static_cast<std::size_t>(lwork));
        int info;

        //Overwrites A with results
        dsyev_(&jobz,
               &uplo,
               &dimension,
               A.data(),
               &dimension,
               result.vector.data(),
               work.data(),
               &lwork, &info
               );


        //Move result to output
        result.matrix = std::move(A);

        return result;

    }

    MatrixAndVector LP_Eig_AB(sqmatrix A, sqmatrix B){
        /// Solves the eigenvector system [Ax=kBx] for x and k
        /// Returns as struct "MatrixAndVector out" with solutions
        /// Stored in out.vector and out.matrix

        assert(A.size()==B.size() && "A and B must be same size in LP_Eig_AB()");

        //Creat dummy storage for outputs
        MatrixAndVector result(A.N());


        //Set flags
        int itype = 1;
        char jobz{'V'};
        char uplo{'U'};
        int dimension = static_cast<int>(A.N());
        int lwork = 6 * dimension;
        std::vector<double> work(static_cast<std::size_t>(lwork));
        int info;

        //Peform lapack routine:
        dsygv_(&itype,
                &jobz,
                &uplo,
                &dimension,
                A.data(),
                &dimension,
                B.data(),
               &dimension,
                result.vector.data(),
                work.data(),
                &lwork,
                &info
                );

        //Move result to output
        result.matrix = std::move(A);

        return result;

    }
}