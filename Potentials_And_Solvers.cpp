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
    double dr = (rgrid[ngrid-1] - rgrid[0]) / (ngrid-1);
    matrix::sqmatrix H(nwaves);

    //Get angular coefficient
    double lam;
    if(l==0){ lam=1.0/2.0; } else { lam=1.0/6.0;}

    //Calculate RHS of exchange potential terms
    list_of_vecs Vex_bj(nwaves);
    for (int j=0; j<nwaves; j++){
        Vex_bj[j].resize(ngrid);
        Vex_bj[j] = YK::ykab(l, P1s_wave, waves[j], rgrid) * P1s_wave * -2.0 * lam;
    }

    for (int i=0; i<nwaves; i++){
        for (int j=0; j<nwaves; j++){
            //Do integrals to get H elements
            H(i,j) = vint(waves[i]*Vex_bj[j], dr);
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

std::vector<energy_and_waves> hartree(const std::vector<double> & rgrid, const std::vector<double> & Vnuc_s, const std::vector<double> & Vnuc_l,
                                      energy_and_waves solutions_s, energy_and_waves solutions_l,
                                      const list_of_vecs & bsplines, const list_of_vecs & bsplines_diff,
                                      int maxits, double tol, int ens_to_check
                                      ){

    //MISSINGNO - Safety checks

    //Grid properties
    int ngrid = rgrid.size();
    double dr = (rgrid[ngrid-1] - rgrid[0]) / (ngrid-1);

    //Voltage vectors
    std::vector<double> V(ngrid);
    std::vector<double> Vdir(ngrid); //For storing e-e interaction

    //Energy convergence properties
    int ittno = 0;
    double echange = tol*10;
    std::vector<double> e_old(ens_to_check);
    for (int i=0; i<ens_to_check; i++){ e_old[i]=solutions_s.energies[i]; }

    std::cout<< "Doing S-Orbital iterations in hartree(): \n \n";
    std::cout<<"ittno"<<"\t"<<"e_1s"<<"\t"<<"e_2s"<<"\t"<<"echange"<<"\n";
    while (ittno < maxits && echange>tol){
        //Get Y^{0}_{1s}{1s}
        Vdir = YK::ykab(0, solutions_s.waves[0],solutions_s.waves[0],rgrid)*2.0;

        //re-calculate potential and use to get solution
        V      = Vnuc_s + Vdir;
        solutions_s = solve_energies(V, bsplines, bsplines_diff, dr);

        //Check changes in energies
        echange=0;
        for (int i=0; i<ens_to_check; i++){
            echange=std::max(echange, fabs(e_old[i]/solutions_s.energies[i]-1));
            e_old[i]=solutions_s.energies[i];
        }

        ittno+=1;
        std::cout<<ittno<<"\t"<<solutions_s.energies[0]<<"\t"<<solutions_s.energies[1]<<"\t"<<echange<<"\n";
    }
    //Use convergent 1s orbital to solve l orbital
    V      = Vnuc_l + Vdir;
    solutions_l = solve_energies(V, bsplines, bsplines_diff, dr);

    //========================================
    return std::vector<energy_and_waves> {solutions_s, solutions_l};
}

std::vector<energy_and_waves> hartree_fast(const std::vector<double> & rgrid, const std::vector<double> & Vnuc_s, const std::vector<double> & Vnuc_l,
                                      energy_and_waves solutions_s, energy_and_waves solutions_l,
                                      int maxits, double tol, int ens_to_check){

    //MISSINGNO - Safety checks

    //Grid properties
    int nwaves = solutions_s.waves.size();
    int ngrid = rgrid.size();
    double dr = (rgrid[ngrid-1] - rgrid[0]) / (ngrid-1);

    //For storing derivatives
    list_of_vecs wavediffs(nwaves);
    for (int i=0; i<nwaves; i++){wavediffs[i].resize(ngrid);}

    //Voltage vectors
    std::vector<double> V(ngrid);
    std::vector<double> Vdir(ngrid); //For storing e-e interaction
    matrix::sqmatrix H(nwaves);

    //Energy convergence properties
    int ittno = 0;
    double echange = tol*10;
    std::vector<double> e_old(ens_to_check);
    for (int i=0; i<ens_to_check; i++){ e_old[i]=solutions_s.energies[i]; }

    std::cout<< "Doing S-Orbital iterations in hartree_fast(): \n \n";
    std::cout<<"ittno"<<"\t"<<"e_1s"<<"\t"<<"e_2s"<<"\t"<<"echange"<<"\n";
    while (ittno < maxits && echange>tol){

        for (int i=0;i<nwaves;i++){ wavediffs[i] = vdiff(solutions_s.waves[i],dr);}

        //Get Y^{0}_{1s}{1s}
        Vdir = YK::ykab(0, solutions_s.waves[0],solutions_s.waves[0],rgrid)*2.0;
        V      = Vnuc_s + Vdir;


        H = make_H(V, solutions_s.waves,wavediffs, dr);
        solutions_s = single_step(H,solutions_s, dr);

        //Check changes in energies
        echange=0;
        for (int i=0; i<ens_to_check; i++){
            echange  = std::max(echange, fabs(e_old[i]/solutions_s.energies[i]-1) );
            e_old[i] = solutions_s.energies[i];
        }

        ittno+=1;
        std::cout<<ittno<<"\t"<<solutions_s.energies[0]<<"\t"<<solutions_s.energies[1]<<"\t"<<echange<<"\n";
    }
    //Use convergent 1s orbital to solve l orbital
    V      = Vnuc_l + Vdir;
    for (int i=0;i<nwaves;i++){ wavediffs[i] = vdiff(solutions_l.waves[i],dr);}
    H = make_H(V, solutions_l.waves,wavediffs, dr);
    solutions_l = single_step(H,solutions_l, dr);

    //========================================
    return std::vector<energy_and_waves> {solutions_s, solutions_l};
}

std::vector<energy_and_waves> hartree_fock(const std::vector<double> & rgrid, const std::vector<double> & Vnuc_s, const std::vector<double> & Vnuc_l,
                                                energy_and_waves solutions_s, energy_and_waves solutions_l,
                                                const list_of_vecs & bsplines, const list_of_vecs & bsplines_diff,
                                                int maxits, double tol, int ens_to_check){
    //MISSINGNO - Safety checks

    //==================================================================
    //Grid properties
    int nwaves = solutions_s.waves.size();
    int ngrid = rgrid.size();
    double dr = (rgrid[ngrid-1] - rgrid[0]) / (ngrid-1);

    //For storing derivatives
    lapack::MatrixAndVector lapack_sols(nwaves);

    //Voltage vectors
    std::vector<double> V(ngrid);
    std::vector<double> Vdir(ngrid); //For storing e-e interaction
    matrix::sqmatrix H(nwaves);
    matrix::sqmatrix B = make_B(bsplines,dr); //Remains constant over all itterations

    //==================================================================

    std::cout<< "Doing S-Orbital iterations in hartree_fock(): \n \n";
    std::cout<<"ittno"<<"\t"<<"e_1s"<<"\t"<<"e_2s"<<"\t"<<"echange"<<"\n";
    int ittno = 0;
    double echange = tol*10;
    std::vector<double> e_old(ens_to_check);
    for (int i=0; i<ens_to_check; i++){ e_old[i]=solutions_s.energies[i]; }

    while (ittno < maxits && echange>tol){

        //Get Y^{0}_{1s}{1s} and update potential:
        Vdir   = YK::ykab(0, solutions_s.waves[0],solutions_s.waves[0],rgrid) * 2.0;
        V      = Vnuc_s + Vdir;

        //Make Hamiltonian and update solutions
        H = make_H(V, bsplines,bsplines_diff, dr) + make_H_exchange(bsplines, solutions_s.waves[0], 0, rgrid);

        lapack_sols = lapack::LP_Eig_AB(H,B);
        solutions_s.energies = std::move(lapack_sols.vector);
        solutions_s.waves    = mix_modes(lapack_sols.matrix, bsplines, dr);

        //Check changes in energies
        echange=0;
        for (int i=0; i<ens_to_check; i++){
            echange=std::max(echange, fabs(e_old[i]/solutions_s.energies[i]-1));
            e_old[i]=solutions_s.energies[i];
        }

        ittno+=1;
        std::cout<<ittno<<"\t"<<solutions_s.energies[0]<<"\t"<<solutions_s.energies[1]<<"\t"<<echange<<"\n";
    }

    //==================================================================

    std::cout<< "Doing P-Orbital iterations in hartree_fock(): \n \n";
    std::cout<<"ittno"<<"\t"<<"e_2s"<<"\t"<<"echange"<<"\n";
    ittno=0;
    echange = tol*2;
    for (int i=0; i<ens_to_check; i++){ e_old[i]=solutions_l.energies[i]; }
    while (ittno < maxits && echange>tol){

        //Vdir is based on S orbital, and so doesn't change
        V      = Vnuc_l + Vdir;

        //Make Hamiltonian and update solutions
        H = make_H(V, bsplines,bsplines_diff, dr) + make_H_exchange(bsplines, solutions_s.waves[0], 1, rgrid);

        lapack_sols = lapack::LP_Eig_AB(H,B);
        solutions_l.energies = std::move(lapack_sols.vector);
        solutions_l.waves    = mix_modes(lapack_sols.matrix, bsplines, dr);

        //Check changes in energies
        echange=0;
        for (int i=0; i<ens_to_check; i++){
            echange = std::max(echange, fabs(e_old[i]/solutions_l.energies[i]-1));
            e_old[i]=solutions_l.energies[i];
        }

        ittno+=1;
        std::cout<<ittno<<"\t"<<solutions_l.energies[0]<<"\t"<<echange<<"\n";
    }

    //========================================
    return std::vector<energy_and_waves> {solutions_s, solutions_l};

}

std::vector<energy_and_waves> hartree_fock_fast(const std::vector<double> & rgrid, const std::vector<double> & Vnuc_s, const std::vector<double> & Vnuc_l,
                                           energy_and_waves solutions_s, energy_and_waves solutions_l,
                                           int maxits, double tol, int ens_to_check){

    //MISSINGNO - Safety checks

    //==================================================================
    //Grid properties
    int nwaves = solutions_s.waves.size();
    int ngrid = rgrid.size();
    double dr = (rgrid[ngrid-1] - rgrid[0]) / (ngrid-1);

    //For storing derivatives
    list_of_vecs wavediffs(nwaves);
    for (int i=0; i<nwaves; i++){wavediffs[i].resize(ngrid);}

    //Voltage vectors
    std::vector<double> V(ngrid);
    std::vector<double> Vdir(ngrid); //For storing e-e interaction
    matrix::sqmatrix H(nwaves);

    //==================================================================

    std::cout<< "Doing S-Orbital iterations in hartree_fock_fast(): \n \n";
    std::cout<<"ittno"<<"\t"<<"e_2s"<<"\t"<<"echange"<<"\n";
    int ittno = 0;
    double echange = tol*10;
    std::vector<double> e_old(ens_to_check);
    for (int i=0; i<ens_to_check; i++){ e_old[i]=solutions_s.energies[i]; }

    while (ittno < maxits && echange>tol){

        for (int i=0;i<nwaves;i++){ wavediffs[i] = vdiff(solutions_s.waves[i],dr);}

        //Get Y^{0}_{1s}{1s} and update potential:
        Vdir = YK::ykab(0, solutions_s.waves[0],solutions_s.waves[0],rgrid) * 2.0;
        V      = Vnuc_s + Vdir;

        //Make Hamiltonian and update solutions
        H = make_H(V, solutions_s.waves,wavediffs, dr) + make_H_exchange(solutions_s.waves, solutions_s.waves[0], 0, rgrid);
        solutions_s = single_step(H,solutions_s, dr);

        //Check changes in energies
        echange=0;
        for (int i=0; i<ens_to_check; i++){
            echange=std::max(echange, fabs(e_old[i]/solutions_s.energies[i]-1));
            e_old[i]=solutions_s.energies[i];
        }

        ittno+=1;
        std::cout<<ittno<<"\t"<<solutions_s.energies[0]<<"\t"<<solutions_s.energies[1]<<"\t"<<echange<<"\n";
    }

    //==================================================================

    std::cout<< "Doing P-Orbital iterations in hartree_fock_fast(): \n \n";
    std::cout<<"ittno"<<"\t"<<"e_2p"<<"\t"<<"echange"<<"\n";
    ittno=0;
    echange = tol*2;
    for (int i=0; i<ens_to_check; i++){ e_old[i]=solutions_l.energies[i]; }
    while (ittno < maxits && echange>tol){

        for (int i=0;i<nwaves;i++){ wavediffs[i] = vdiff(solutions_l.waves[i],dr);}

        //Vdir is based on S orbital, and so doesn't change
        V      = Vnuc_l + Vdir;

        //Make Hamiltonian and update solutions
        H = make_H(V, solutions_l.waves,wavediffs, dr) + make_H_exchange(solutions_l.waves, solutions_s.waves[0], 1, rgrid);
        solutions_l = single_step(H,solutions_l, dr);

        //Check changes in energies
        echange=0;
        for (int i=0; i<ens_to_check; i++){
            echange = std::max(echange, fabs(e_old[i]/solutions_l.energies[i]-1));
            e_old[i]=solutions_l.energies[i];
        }

        ittno+=1;
        std::cout<<ittno<<"\t"<<solutions_l.energies[0]<<"\t"<<echange<<"\n";
    }

    //========================================
    return std::vector<energy_and_waves> {solutions_s, solutions_l};
}
