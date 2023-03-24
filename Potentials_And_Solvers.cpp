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
#include "spline_eval.hpp"

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
//Compact procedures / matrix construction

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


matrix::sqmatrix make_H_exchange(const list_of_vecs & waves, const std::vector<double> & P1s_wave, int l, const std::vector<double> & rgrid){
    /// Constructs the hamiltonian matrix term that encodes the exchange potential

    int nwaves = waves.size();
    int ngrid  = rgrid.size();
    double dr = (rgrid[ngrid] - rgrid[0]) / (ngrid-1);
    matrix::sqmatrix H(nwaves);

    //Get angular coefficient
    double lam;
    if(l==0){lam=1/2;} else { lam=1/6;}

    //Calculate RHS of exchange potential terms
    list_of_vecs Vex_bj(nwaves);
    for (int i=0; i<nwaves; i++){
        Vex_bj[i].resize(waves[0].size());
        Vex_bj[i] = -2*lam*P1s_wave * YK::ykab(l, P1s_wave, waves[i], rgrid);

    }

    for (int i=0; i<nwaves; i++){
        for (int j=i; j<nwaves; j++){

            //Do integrals to get H elements
            H(i,j) = vint(waves[i]*Vex_bj[j], dr);

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


energy_and_waves single_step(const matrix::sqmatrix & H, energy_and_waves system, double dr){
    /// Performs a single update on an ortho-normal set of waves
    // [MISSINGNO] Safety checks

    // Use linalg to get solutions
    lapack::MatrixAndVector solutions = lapack::LP_Eig_A(H);

    // mix_modes
    system.energies = std::move(solutions.vector);
    system.waves    = mix_modes(solutions.matrix, system.waves, dr);

    // return
    return(system);
}

std::vector<energy_and_waves> hartree(std::vector<double> Vnuc_s, std::vector<double> Vnuc_l, std::vector<double> rgrid,
                                      energy_and_waves s_sols, energy_and_waves l_sols,
                                      const  list_of_vecs & bsplines, const list_of_vecs & bsplines_diff,
                                      int maxits, double tol, int ens_to_check){


}


std::vector<energy_and_waves> hartree_fock(std::vector<double> Vnuc_s, std::vector<double> Vnuc_l, std::vector<double> rgrid,
                              energy_and_waves s_sols, energy_and_waves l_sols,
                              int maxits, double tol, int ens_to_check){

    // [MISSINGNO] Safety checks and add const refs

    //===============================

    //Setup
    int ngrid = rgrid.size();
    int nwaves = s_sols.energies.size();
    double dr = (rgrid[ngrid] - rgrid[0]) / (ngrid-1);

    list_of_vecs bsplines       =  generate_splines(nwaves, rgrid);
    list_of_vecs bsplines_diff  =  generate_spline_diffs(nwaves, rgrid);

    //Objects to itterate in loop
    //list_of_vecs & waves = s_sols.waves;
    list_of_vecs wave_diffs(nwaves);
    //std::vector<double> & P1s_wave = waves[0];

    matrix::sqmatrix H(nwaves);
    std::vector<double> V(ngrid);
    std::vector<double> Vdir(ngrid);

    //Convergence loop properties
    int ittno;
    double e_diff;
    std::vector<double> e_old(ens_to_check);

    //===============================
    //S orbital calculation
    std::cout << "Beginning S orbital itterations \n";
    ittno=0;
    e_diff=tol*2;
    for (int i=0; i<ens_to_check; i++){ e_old[i]=s_sols.energies[i];}

    while(ittno<maxits){
        e_diff=0;
        //Get derivatives of waves
        for(int i=0; i<nwaves; i++){
            //wave_diffs[i]= vdiff(s_sols.waves[i],dr);
        }

        //Calculate Direct potential
        Vdir = YK::ykab(0, s_sols.waves[0],s_sols.waves[0],rgrid) * 2;
        V = Vnuc_s + Vdir;

        //Make Hamiltonian
        //H = make_H(V, waves, wave_diffs, dr);// + make_H_exchange(s_sols.waves, s_sols.waves[0], 0, rgrid);

        //Do Update
        //s_sols = single_step(H, s_sols, dr);

        s_sols = solve_energies(V, bsplines, bsplines_diff, dr);

        //Check changes in energies
        e_diff=0;
        for (int i=0; i<ens_to_check; i++){
            e_diff=std::max(e_diff, fabs(e_old[i]/s_sols.energies[i]-1));
            e_old[i]=s_sols.energies[i];
        }

        //Advance and output
        ittno+=1;
        std::cout<<"ittno:\t"<<ittno<<"\t Fractional change:\t"<<e_diff<<"\n";
        printv(e_old);

        }

    //===============================
    //L orbital calculation
    std::cout << "Beginning L orbital itterations \n";
    //waves = l_sols.waves;
    ittno = 0;
    e_diff=tol*2;
    for (int i=0; i<ens_to_check; i++){ e_old[i]=l_sols.energies[i];}
    Vdir = YK::ykab(0, s_sols.waves[0],s_sols.waves[0],rgrid)  * 2;
    while(ittno<maxits){
        e_diff=0;
        //Get derivatives of waves
        std::cout << ittno << "\n";
        for(int i=0; i<nwaves; i++){wave_diffs[i]= vdiff(l_sols.waves[i],dr);}

        //Calculate potential
        V = Vnuc_l + Vdir;

        //Make Hamiltonian
        H = make_H(V, l_sols.waves, wave_diffs, dr);// + make_H_exchange(l_sols.waves, s_sols.waves[0], 1, rgrid);

        //Do Update
        l_sols = single_step(H, l_sols, dr);

        //[MISSINGNO] - energy checks

        ittno+=1;

        //Check changes in energies
        e_diff=0;
        for (int i=0; i<ens_to_check; i++){
            e_diff=std::max(e_diff, fabs(e_old[i]/l_sols.energies[i]-1));
            e_old[i]=l_sols.energies[i];
        }

        //Advance and output
        ittno+=1;
        std::cout<<"ittno:\t"<<ittno<<"\t Fractional change:\t"<<e_diff<<"\n";
    }

    //===========
    std::vector<energy_and_waves> out = {s_sols,l_sols};
    return out;
}

