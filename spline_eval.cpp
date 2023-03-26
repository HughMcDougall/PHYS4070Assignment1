//
// Functions for easily generating splines. Contains wrappers for bspline.cpp
// HM Mar 23
//

#include "spline_eval.hpp"
//Imports
//  System
#include <vector>
#include <functional>

//  Project Local
#include "LP_solvers.hpp"
#include "bspline.hpp"
#include "matrix.hpp"
#include "vector_utils.hpp"

using list_of_vecs = std::vector<std::vector<double>>;
//=======================================================
//B-Spline Evaluation
list_of_vecs generate_splines(int nsplines, const std::vector<double> & rgrid, int k_spline){
    /// Wrapper function for generating splines quickly
    /// Takes nsplines, evenly spaced sorted grid of radii, and order of splines (defaults to 7)
    //Sanity Checks
    assert(nsplines && "Too few splines for spline order. Require nsplines + 3> 2 * k_spline");
    assert(rgrid.size()> 2 && "Invalid radial grid in generating bsplines");

    double rmin = rgrid.front();
    double rmax = rgrid.back();
    int n_grid = rgrid.size();

    //Output Storage (vec of vecs)
    list_of_vecs out(nsplines);
    //Init Bspline object. '+3' because we skip first 2 and last 1 bplsines
    BSpline bspl(k_spline, nsplines+3, rmin, rmax);

    //For each spline:
    double r;
    for (int n=0;   n<nsplines;  n++){
        out[n].resize(n_grid); //Set to correct size

        //For each gridpoint
        for (int k=0;   k<n_grid; k++){
            r         = rgrid[k];
            out[n][k] = bspl.b(n+2, r); // Eval and store spline. '+2' because we skip the first two splines
        }
    }
    return(out);
}

list_of_vecs generate_spline_diffs(int nsplines, const std::vector<double> & rgrid, int k_spline, double dr){
    /// Wrapper function for generating splines derivatives quickly
    /// Takes nsplines, evenly spaced sorted grid of radii, and order of splines (defaults to 7)

    //Sanity Checks
    assert(nsplines && "Too few splines for spline order. Require nsplines + 3> 2 * k_spline");
    assert(rgrid.size()> 2 && "Invalid radial grid in generating bsplines");
    double rmin = rgrid.front();
    double rmax = rgrid.back();
    int n_grid = rgrid.size();

    //If no dr provided, estimate from ngrid
    if (dr==-1){dr = (rmax-rmin) / (n_grid-1);}

    //Output Storage (vec of vecs)
    list_of_vecs out(nsplines);
    //Init Bspline object. '+3' because we skip first 2 and last 1 bplsines
    BSpline bspl(k_spline, nsplines+3, rmin, rmax);

    //For each spline:
    double r;
    for (int n=0;   n<nsplines;  n++){
        out[n].resize(n_grid); //Set to correct size

        //For each gridpoint
        for (int k=0;   k<n_grid; k++){
            r = rgrid[k];
            // Eval and store spline. '+2' because we skip the first two splines
            out[n][k] = (bspl.b(n+2, r+dr/2) - bspl.b(n+2, r-dr/2)) / dr;
        }
    }
    return(out);
}

