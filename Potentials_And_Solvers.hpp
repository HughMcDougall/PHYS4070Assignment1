//
// Created by hughm on 18/03/2023.
//

#include <vector>
#ifndef ASSIGNMENT1_POTENTIALS_AND_SOLVERS_H
#define ASSIGNMENT1_POTENTIALS_AND_SOLVERS_H
//=======================================================
//Aliases & Holder Structs
using list_of_vecs = std::vector<std::vector<double>>;

double V_hydrogen(double r, int Z, int l);
double V_green(double r, int Z,  double d, double h );
matrix::sqmatrix make_H(const std::vector<double> & V, const list_of_vecs & waves, const list_of_vecs & waves_diffs, double dr = 0);
matrix::sqmatrix make_H_exchange(const list_of_vecs & waves, const std::vector<double> & P1s_wave, int l, const std::vector<double> & rgrid);
matrix::sqmatrix make_B(const list_of_vecs & waves, double dr=0);
list_of_vecs mix_modes(const matrix::sqmatrix & mat, const list_of_vecs & old_modes, double dr=0);

struct energy_and_waves{
    /// Holder class for our solutions
    energy_and_waves(int n): N(n), energies(n), waves(n){}
    int N;
    std::vector<double> energies;
    list_of_vecs waves;
};


energy_and_waves solve_energies(std::vector<double> & V, const list_of_vecs & splines, const list_of_vecs & spline_diff, double dr=0);
energy_and_waves single_step(const matrix::sqmatrix & H, energy_and_waves system, double dr=0);

std::vector<energy_and_waves> hartree(const std::vector<double> & rgrid,  const std::vector<double> & Vnuc_s, const std::vector<double> & Vnuc_l,
                                      energy_and_waves solutions_s, energy_and_waves solutions_l,
                                      const list_of_vecs & bsplines, const list_of_vecs & bsplines_diff,
                                      int maxits = 40, double tol = 0.01, int ens_to_check = 5
                                      );

std::vector<energy_and_waves> hartree_fast(const std::vector<double> & rgrid,  const std::vector<double> & Vnuc_s, const std::vector<double> & Vnuc_l,
                                      energy_and_waves solutions_s, energy_and_waves solutions_l,
                                      int maxits = 40, double tol = 0.01, int ens_to_check = 5
);

std::vector<energy_and_waves> hartree_fock(const std::vector<double> & rgrid,  const std::vector<double> & Vnuc_s, const std::vector<double> & Vnuc_l,
                                                energy_and_waves solutions_s, energy_and_waves solutions_l,
                                                const list_of_vecs & bsplines, const list_of_vecs & bsplines_diff,
                                                int maxits = 40, double tol = 0.01, int ens_to_check = 5
);


std::vector<energy_and_waves> hartree_fock_fast(const std::vector<double> & rgrid,  const std::vector<double> & Vnuc_s, const std::vector<double> & Vnuc_l,
                                           energy_and_waves solutions_s, energy_and_waves solutions_l,
                                           int maxits = 40, double tol = 0.01, int ens_to_check = 5
);
#endif //ASSIGNMENT1_POTENTIALS_AND_SOLVERS_H
