/*
 * PHYS4070 Assignment 01
 * Bsc Phys Hon Semester 1 2023
 * Hugh McDougall- 43202007
*/

//Imports
//  System
#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <fstream>

//  Project Local
#include "LP_solvers.hpp"
#include "bspline.hpp"
#include "matrix.hpp"
// #include "vector_utils.hpp"

//=======================================================
//Aliases & Holder Structs
using function_1D  = std::function<double(double)>;
using list_of_vecs = std::vector<std::vector<double>>;

struct energy_and_waves{
    /// Holder class for our solutions
    energy_and_waves(int n): N(n), energies(n), waves(n){}
    int N;
    std::vector<double> energies;
    list_of_vecs waves;
};
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

//=======================================================
//Physically meaningful functions
double V_hydrogen(double r, int Z, int l){
    return -Z/r + l*(l+1)/r/r/2;
}

energy_and_waves solve_energies(std::vector<double> & V, const list_of_vecs & splines, const list_of_vecs & spline_diff){
    /// Function to get energies and energy eigenstates of a potential

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


//=======================================================
//Main runtime
int main(int argc, char *argv[]) {
    /// Read Inputs
    std::cout << "Program Starting\n";
    int nsplines, ngrid;
    double rmin, rmax;
    int z, l;

    // Read inputs. Scope to allow to collapse easily in editor
    {
        if (argc > 1){
            nsplines = std::stoi(argv[1]);
        } else{
            nsplines = 60;
        }

        if (argc > 2){
            ngrid = std::stoi(argv[2]);
        } else{
            ngrid = 1001;
        }

        if (argc > 3){
            rmin = std::stoi(argv[3]);
        } else{
            rmin = 0.001;
        }

        if (argc > 4){
            rmax = std::stoi(argv[4]);
        } else{
            rmax  = 100;
        }

    }


    if (true){
    //-----------------------------------------
    //Make grid, b-splines & Diffs
    std::cout << "Generating " << nsplines << " splines...\n";

    std::vector<double> rgrid = make_rgrid(rmin, rmax, ngrid);

    list_of_vecs bsplines       =  generate_splines(nsplines, rgrid);
    list_of_vecs bsplines_diff  =  generate_spline_diffs(nsplines, rgrid);

    //Save bsline graphs
    {
        std::ofstream spline_file("./outputs/splines.txt");
        std::ofstream spline_diffile("./outputs/spline_diffs.txt");
        for (int i = 0; i < ngrid; i++) {

            spline_file << rgrid[i] << "\t ";
            spline_diffile << rgrid[i] << "\t ";

            for (int j = 0; j < nsplines; j++) {
                spline_file << bsplines[j][i] << "\t ";
                spline_diffile << bsplines_diff[j][i] << "\t ";
            }
            spline_file << "\n";
            spline_diffile << "\n";
        }
        spline_file.close();
        spline_diffile.close();
        }

    std::cout << "Done.\n";

    //-----------------------------------------
    //Make voltage vector
    int Z = 1;
    int l = 0;

    std::cout << "Generating Voltage potential for Z= " << Z << " and l= " << l << "\n";
    auto Vfunc = [Z,l](double r){ return V_hydrogen(r, Z, l); }; // Make potential func
    std::vector<double> V = potvec_from_func(Vfunc, rgrid);

    int n=0;

    //Save voltage graph
    {
        std::ofstream potfile("./outputs/potential.txt");
        for (int i=0; i<ngrid; i++){
            potfile << rgrid[i] << "\t ";
            potfile << V[i];
            potfile << "\n";
        }

        potfile.close();
    }

    std::cout << "Done.\n";

    //-----------------------------------------
    //eval energies
    std::cout << "Calculating energy eigenstates...\n";

    energy_and_waves energy_sols = solve_energies(V, bsplines, bsplines_diff);

    std::cout << "Calculated energies:\t";
    printv(energy_sols.energies);

    //Save solutions
    {
        std::ofstream wavefile("./outputs/waves.txt");
        std::ofstream energyfile("./outputs/energies.txt");
        for (int i = 0; i < ngrid; i++) {

            wavefile << rgrid[i] << "\t ";

            for (int j = 0; j < nsplines; j++) {
                wavefile << energy_sols.waves[j][i] << "\t ";
            }
            wavefile << "\n";
        }

        for (int j = 0; j < nsplines; j++) {
            energyfile << j+1 << "\t";
            energyfile << energy_sols.energies[j] << "\n";
        }

        energyfile.close();
        wavefile.close();
    }

    std::cout << "Done.\n";
    //-----------------------------------------
    std::cout << "\t\t Program Done.\n";
    } else{
        /* This is a test section to make sure the matrix solving is behaving.
         * This shouldn't run at all. Delete before final ship
         * */

        std::vector<double> a{1,1,1,1,1,1,2};
        std::vector<double> b{2,2,4,5,6,7,3};
        std::vector<double> c{3,-1,-1,-1,-1,-1,1};

        printv(a);
        printv(a*b);
        printv(a+c);
        printv(a);

        matrix::sqmatrix A(2);
        A(0,0)=1;
        A(1,1)=2;
        A.print();
        A = 0.5*A;
        std::cout<<"0.5*A = \n";
        A.print();

    }

    return 0;
}
