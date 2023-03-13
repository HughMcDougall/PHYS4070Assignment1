//
// Created by hughm on 13/03/2023.
//

#include <iostream>
#include <vector>
#include <cassert>
#include <functional>

#include "bspline.hpp"

using function_1D  = std::function<double(double)>;
using list_of_vecs = std::vector<std::vector<double>>;

//=======================================================
//Utility Functions
void printv(const std::vector<double>& a){
    for (int i=0; i<a.size(); i++){
        std::cout << a[i] << ", ";
    }
    std::cout<<"\n";
}

//Overload vector operations to make direct products easier
//V-V Multiplication
std::vector<double> operator*=(std::vector<double> & a, const std::vector<double> & b){
    assert(a.size()==b.size() && "Tried to multiply two vectors with different lengths");
    for (int i=0; i<a.size(); i++){
        a[i]*=b[i];
    }
    return(a);
}
std::vector<double> operator*(std::vector<double> a, const std::vector<double> & b){return a*=b;}

//V-V Addition
std::vector<double> operator+=(std::vector<double> & a, const std::vector<double> & b){
    assert(a.size()==b.size() && "Tried to add two vectors with different lengths");
    for (int i=0; i<a.size(); i++){
        a[i]+=b[i];
    }
    return(a);
}
std::vector<double> operator+(std::vector<double> a, const std::vector<double> & b){return a+=b;}

//V-D Multiplication
std::vector<double> operator*(double & a, std::vector<double> v){
    for (int i=0; i<v.size(); i++){
        v[i]*=a;
    }
    return(v);
}
std::vector<double> operator*(std::vector<double> v, double & a){return(a*v);}

//Conversions & Integrations
//Function to integrate over a vector
double vint(const std::vector<double>& a){
    //Integrates over a vector. For integration purposes. Uses simspsons rule;
    int n = a.size();
    if (n%2==0){std::cerr<<"Warning! Integrating vector with even number of points. Last point dropped. May be imprecise";}
    double out = 0;

    //Peform weighted summation in keeping /w simpsons rule
    out += a.front();
    for (int i=1; i<n; i+=2){ //Odd Indices
        out+=a[i]*4;
    }
    for (int i=2; i<n; i+=2){ //Even Indices
        out+=a[i]*2;
    }
    out += a.back();

    out /= 3;

    return out;
}

std::vector<double> make_rgrid(double rmin = 0.001, double rmax = 100, int n_grid = 101){
    /// Creates a linspace of radial grid points. Saves time on passing the same args to every function
    double dr = (rmax-rmin) / (n_grid-1);
    double r = rmin;
    std::vector<double> out(n_grid);

    for (int i=0; i<n_grid; i++){
        out[i]=r;
        r+=dr;
    }

    return out;
}

std::vector<double> potvec_from_func(const function_1D& V, const std::vector<double>& rgrid){
    ///Converts a 1D potential function into a std::vector for quick-swapping of potential functions
    int n_grid = rgrid.size();
    std::vector<double> out (n_grid);

    for (int i=0; i<n_grid; i++){
        out[i]=V(rgrid[i]);
    }

    return out;
}
//=======================================================
//B-Spline Evaluation
list_of_vecs generate_splines(int nsplines, const std::vector<double> & rgrid, int k_spline = 7){
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

list_of_vecs generate_spline_diffs(int nsplines, const std::vector<double> & rgrid, int k_spline = 7){
    //Sanity Checks
    assert(nsplines && "Too few splines for spline order. Require nsplines + 3> 2 * k_spline");
    assert(rgrid.size()> 2 && "Invalid radial grid in generating bsplines");
    double rmin = rgrid.front();
    double rmax = rgrid.back();
    int n_grid = rgrid.size();

    double dr = (rmax-rmin) / (n_grid-1);

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