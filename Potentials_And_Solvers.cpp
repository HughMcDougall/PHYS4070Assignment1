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
#include "calculateYK.hpp"

//=======================================================
//Aliases & Holder Structs
using function_1D  = std::function<double(double)>;
using list_of_vecs = std::vector<std::vector<double>>;
//=======================================================
//Potentials
double V_hydrogen(double r, int Z, int l){
    /// A simple hydrogen-like potential
    return -Z/r + l*(l+1)/r/r/2;
}

double V_green(double r, int Z,  double d, double h ){
    /// Green's function, a first estimate of the electron-electron potential term from mean field theory
    return (Z-1) / r * h * (exp(r/d)-1 ) / (1+h*(exp(r/d)-1));
}

//=======================================================
//Solvers

energy_and_waves solve_energies(std::vector<double> & V, const list_of_vecs & splines, const list_of_vecs & spline_diff, double dr){
    /// Function to get energies and energy eigenstates of a potential
    /// Takes voltage vector 'V' and bspline + bspline_diff evaluations

    // Sanity checks for inputs using assert
    assert(splines.size()==spline_diff.size() && "Mismatch in spline and spline_diff numbers");
    int nsplines = splines.size();
    for (int i=0; i<nsplines; i++){
        assert(V.size()==splines[i].size() && V.size() == spline_diff[i].size() && "Length of all splines and voltage must match");
    }

    //===========
    // MAKE MATRICES
    matrix::sqmatrix H = make_H(V,splines,spline_diff); //hamiltonian matrix
    matrix::sqmatrix B = make_B(splines); //bspline diagonalization matrix

    //===========
    // USE LAPACK TO GET EIGENSTATES
    lapack::MatrixAndVector solutions = lapack::LP_Eig_AB(H,B);

    //===========
    // ASSEMBLE INTO WAVEFUNCTIONS AND OUTPUT
    energy_and_waves out(nsplines); // empty object

    out.energies = std::move(solutions.vector); //Save eigenstates as energies
    out.waves    = mix_modes(solutions.matrix, splines, dr);

    //===========
    return out;
}

matrix::sqmatrix make_H(const std::vector<double> & V, const list_of_vecs & waves, const list_of_vecs & waves_diffs, double dr) {
    /// Makes the non-exchange hamiltonian matrix for a set of waves and their derivatives
    /// Assumes dr
    int nwaves = waves.size();
    matrix::sqmatrix H(nwaves);

    for (int i=0; i<nwaves; i++){
        for (int j=i; j<nwaves; j++){

            //Do integrals to get H elements
            H(i,j) = 0.5 * vint(waves_diffs[i] * waves_diffs[j], dr) + vint(waves[i] * V * waves[j], dr);

            //Fill lower triangular terms
            if (i!=j){
                H(j,i) = H(i,j);
            }
        }
    }
    return(H);
}

matrix::sqmatrix make_B(const list_of_vecs & waves, double dr) {
    /// Makes the mode mixing matrix for a set of basis waves
    int nwaves = waves.size();
    matrix::sqmatrix B(nwaves);

    for (int i=0; i<nwaves; i++){
        for (int j=i; j<nwaves; j++){

            //Do integrals to get H elements
            B(i,j) = vint(waves[i] * waves[j], dr);

            //Fill lower triangular terms
            if (i!=j){
                B(j,i) = B(i,j);
            }
        }
    }
    return(B);
}

list_of_vecs mix_modes(const matrix::sqmatrix & mat, const list_of_vecs & old_modes, double dr){
    /// A function to quickly combine a coefficient matrix and a set of basis functions into a new set of solutions

    int nmodes= old_modes.size();
    int ngrid = old_modes[0].size();
    list_of_vecs out(nmodes);

    for (int i=0; i<nmodes; i++){

        // Set correct size
        out[i].resize(ngrid);

        // Do dot product between basis vectors and coefficient matrix
        for (int j=0; j<nmodes; j++){
            out[i]    +=  mat(i,j) * old_modes[j];
        }

        // If spatial element provided, normalize output
        if(dr!=0) {
            out[i] /= pow(vint(out[i]*out[i], dr),0.5);
        }

    }

    return out;
}

std::vector<energy_and_waves> hartree(const std::vector<double> & Vnuc_s,const std::vector<double> & Vnuc_l, std::vector<double> & rgrid,
                         const energy_and_waves & initial_solution_sorbital, const energy_and_waves & initial_solution_lorbital,
                         int maxits, double tol){
    /// Performs the Hartree method itteratively to get a self consistent solution for the s orbitals and applies to the l orbitals
    /// 'previous_solutions' for first itteration should already be an orthogonal basis from solve_energies()

    // Asserts and number of evals
    assert (initial_solution_sorbital.waves.size() == initial_solution_lorbital.waves.size() && "Initial solutions in hartree() must have same wave number");
    int nwaves = initial_solution_sorbital.waves.size();

    int ngrid  = rgrid.size();
    assert (Vnuc_s.size()==Vnuc_l.size() && "Potential vectors in hartree() must be same length as rgrid");
    for (int i=0; i<nwaves; i++){
        assert(initial_solution_sorbital.waves[i].size() ==  ngrid && initial_solution_lorbital.waves[i].size() == ngrid && "All initial sol waves in hartree() must be save dimension as rgrid");
    }

    assert(maxits>0 && tol<1 && "Incorrect hartree convergence limits in hartree()");

    //------------------------------------------

    //storage for Output, two save spaces for solution storing struc
    std::vector<energy_and_waves> out(2, energy_and_waves(nwaves));
    out[0].energies=vcopy(initial_solution_sorbital.energies);

    //Establish number of evals
    double dr  = (rgrid.back()-rgrid.front()) / (ngrid-1);

    //Make storage spaces for solutions
    list_of_vecs waves(nwaves);
    list_of_vecs waves_new(nwaves);
    list_of_vecs waves_diffs(nwaves);

    //l_orbitals
    list_of_vecs wavel(nwaves);
    list_of_vecs wavel_diffs(nwaves);

    for (int i=0; i<nwaves; i++){
        waves[i] =vcopy(initial_solution_sorbital.waves[i]);
        waves_new[i].resize(ngrid);
        waves_diffs[i].resize(ngrid);

        //Don't itterate l orbitals in hartree so don't need wavel_new
        wavel[i].resize(ngrid);
        wavel_diffs[i] = vdiff(initial_solution_lorbital.waves[i],dr);
    }

    //Make vectors and matrices for use in loop
    std::vector<double> V(ngrid);
    std::vector<double> y0_1s1s(ngrid);
    matrix::sqmatrix H(nwaves);
    matrix::sqmatrix B(nwaves);

    //===========
    // S ORBITAL
    //--------------------------------
    //Do Hartree Itterations to get s orbitals
    int ittno = 0;
    bool converged = false;
    lapack::MatrixAndVector solutions(nwaves); // object for storing outputs of lapack functions

    //Main Hartree Loop for solving s orbital
    while (ittno<maxits && not converged){
        ittno+=1;

        //Get y^{0}_{1s1s} and wave derivatives for populating H matrix
        for (int i=0; i<nwaves; i++){ waves_diffs[i] = vdiff(waves[i], dr);} // Get derivatives
        y0_1s1s = YK::ykab(0, waves[0],waves[0],rgrid);
        V = Vnuc_s + 2 * y0_1s1s;

        //Loop over upper triangular terms to populate H matrix
        H = make_H(V, waves, waves_diffs, dr);
        B = make_B(waves, dr);

        //Solve using lapack dsyev (basis should already be orthogonal)
        solutions = lapack::LP_Eig_AB(H,B);

        //Recombine old waves into new solutions
        for (int i=0; i<nwaves; i++){
            waves_new[i]*=0; // Zero out old elements
            for (int j=0; j<nwaves; j++){
                waves_new[i]    +=  solutions.matrix(i,j) * waves[j];
            }
        }

        // Check convergence based on energies
        std::cout<<ittno<<"\t";
        converged=true;
        for( double el : solutions.vector/out[0].energies-1 ){
            if (abs(el)>tol){
                converged=false;
            }
            std::cout<< el<<", ";
        }
        std::cout<< "\n";



        for (int i=0; i<nwaves; i++){waves[i] = waves[i] *(vcopy(waves_new[i])+waves[i]);} // update 'old' wave to new wave. BROKEN. WAS USING COPY
        out[0].energies=vcopy(solutions.vector);
    } //Itterations

    //Save S-orbital results
    out[0].energies = std::move(solutions.vector); //Save eigenstates as energies
    out[0].waves    = waves;

    //===========
    // L ORBITAL
    // Use convergent 1s orbital to calculate l orbitals
    y0_1s1s = YK::ykab(0, waves[0],waves[0],rgrid);
    V = Vnuc_l + 2 * y0_1s1s;
    for (int i=0; i<nwaves; i++){
        for (int j=i; j<nwaves; j++){

            //Do integrals. No 'dr' integral element as it cancels out on either side of the eigenstate eqn
            H(i,j) = 0.5 * vint(wavel_diffs[i] * wavel_diffs[j],dr) + vint(wavel[i] * V * wavel[j], dr);

            //Fill lower triangular terms
            if (i!=j){
                H(j,i) = H(i,j);
            }
        }
    }
    solutions = lapack::LP_Eig_A(H);

    //Recombine old waves into new solutions
    for (int i=0; i<nwaves; i++){
        for (int j=0; j<nwaves; j++){
            wavel[i]    +=  solutions.matrix(i,j) * initial_solution_lorbital.waves[j];
        }
    }

    //Save l-orbital results
    out[1].energies = std::move(solutions.vector); //Save eigenstates as energies
    out[1].waves    = wavel;

    return out;




}