//
// Created by hughm on 18/03/2023.
//

//Imports
//  System
#include <vector>
#include <functional>

//  Project Local
#include "LP_solvers.hpp"
#include "bspline.hpp"
#include "matrix.hpp"
#include "vector_utils.hpp"
#include "Potentials_And_Solvers.hpp"
//=======================================================
//Aliases & Holder Structs
using function_1D  = std::function<double(double)>;
using list_of_vecs = std::vector<std::vector<double>>;
//=======================================================
//Physically meaningful functions
double V_hydrogen(double r, int Z, int l){
    /// A simple hydrogen-like potential
    return -Z/r + l*(l+1)/r/r/2;
}

double V_green(double r, int Z,  double d, double h ){
    /// Green's function, a first estimate of the electron-electron potential term from mean field theory
    return (Z-1) / r * h * (exp(r/d)-1 ) / (1+h*(exp(r/d)-1));
}
energy_and_waves solve_energies(std::vector<double> & V, const list_of_vecs & splines, const list_of_vecs & spline_diff){
    /// Function to get energies and energy eigenstates of a potential
    /// Takes voltage vector 'V' and bspline + bspline_diff evaluations

    // Sanity checks for inputs using assert
    assert(splines.size()==spline_diff.size() && "Mismatch in spline and spline_diff numbers");
    int nsplines = splines.size();
    for (int i=0; i<nsplines; i++){
        assert(V.size()==splines[i].size() && V.size() == spline_diff[i].size() && "Length of all splines and voltage must match");
    }
    int ngrid = V.size();

    //===========
    // MAKE MATRICES
    matrix::sqmatrix H(nsplines); //hamiltonian matrix
    matrix::sqmatrix B(nsplines); //bspline diagonalization matrix

    //Loop over upper triangular terms
    for (int i=0; i<nsplines; i++){
        for (int j=i; j<nsplines; j++){

            //Do integrals. No 'dr' integral element as it cancels out on either side of the eigenstate eqn
            H(i,j) = 0.5 * vint(spline_diff[i] * spline_diff[j]) + vint(splines[i] * V * splines[j]);
            B(i,j) = vint(splines[i] * splines[j]);

            //Fill lower triangular terms
            if (i!=j){
                H(j,i) += H(i,j);
                B(j,i) += B(i,j);
            }
        }
    }

    //===========
    // USE LAPACK TO GET EIGENSTATES
    lapack::MatrixAndVector solutions = lapack::LP_Eig_AB(H,B);

    //===========
    // ASSEMBLE INTO WAVEFUNCTIONS AND OUTPUT
    energy_and_waves out(nsplines); // empty object

    out.energies = std::move(solutions.vector); //Save eigenstates as energies
    for (int i=0; i<nsplines; i++){ //Recombine bsplines using coefficients to get wavefunctions
        out.waves[i].resize(ngrid);
        for (int j=0; j<nsplines; j++){
            out.waves[i]    +=  solutions.matrix(i,j) * splines[j];
        }
    }
    return out;
}