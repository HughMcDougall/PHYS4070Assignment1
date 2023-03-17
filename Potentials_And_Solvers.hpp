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
struct energy_and_waves{
    /// Holder class for our solutions
    energy_and_waves(int n): N(n), energies(n), waves(n){}
    int N;
    std::vector<double> energies;
    list_of_vecs waves;
};
energy_and_waves solve_energies(std::vector<double> & V, const list_of_vecs & splines, const list_of_vecs & spline_diff);

#endif //ASSIGNMENT1_POTENTIALS_AND_SOLVERS_H
