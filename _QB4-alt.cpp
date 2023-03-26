/*
 * PHYS4070 Assignment 01
 * Bsc Phys Hon Semester 1 2023
 * Hugh McDougall - 43202007
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
#include "vector_utils.hpp"
#include "spline_eval.hpp"
#include "Potentials_And_Solvers.hpp"
#include "calculateYK.hpp"

//=======================================================
//Aliases & Holder Structs
using function_1D  = std::function<double(double)>;
using list_of_vecs = std::vector<std::vector<double>>;

//=======================================================
//Main runtime
int main(int argc, char *argv[]) {
    /// Read Inputs
    /// Order of inputs : nsplines,   ngrid,   rmin,   rmax,   save_outputs,   output_dir       maxits  tol     ens_to_check
    /// Default inputs  : 60,         5001,   0.001,  100,     False,           "./outputs_B4"  40      1E-6    5

    std::cout << "------------------\n";
    std::cout << "------------------\n";
    std::cout << "Starting Question B4\n";
    std::cout << "------------------\n";
    std::cout << "------------------\n";
    int nsplines, ngrid;
    double rmin, rmax;
    int maxits, ens_to_check;
    double tol;
    bool save_outputs;
    std::string output_dir;

    // Read inputs. Scoped to allow to collapse easily in editor
    {
        if (argc > 1){
            nsplines = std::stoi(argv[1]);
        } else{
            nsplines = 60;
        }

        if (argc > 2){
            ngrid = std::stoi(argv[2]);
        } else{
            ngrid = 10001;
        }

        if (argc > 3){
            rmin = std::stof(argv[3]);
        } else{
            rmin = 0.00001;
        }

        if (argc > 4){
            rmax = std::stof(argv[4]);
        } else{
            rmax  = 100;
        }

        if (argc > 5){
            maxits = std::stoi(argv[5]);
        } else{
            maxits  = 100;
        }

        if (argc > 6){
            tol = std::stof(argv[6]);
        } else{
            tol  = 0.000001;
        }

        if (argc > 7){
            ens_to_check = std::stoi(argv[7]);
        } else{
            ens_to_check  = 5;
        }

        if (argc > 8){
            save_outputs = std::stoi(argv[8]);
        } else{
            save_outputs  = false;
        }

        if (argc > 9){
            output_dir = argv[9];
        } else{
            output_dir  = "./outputs/B4/";
        }

        std::cout << "Parameters:\n";
        std::cout << "nsplines = \t"    << nsplines << "\n";
        std::cout << "ngrid = \t"       << ngrid << "\n";
        std::cout << "rmin = \t"    << rmin << "\n";
        std::cout << "rmax = \t"       << rmax << "\n";
        std::cout << "maxits = \t"    << maxits << "\n";
        std::cout << "tol = \t"       << tol << "\n";
        std::cout << "ens_to_check = \t"       << ens_to_check << "\n";
        std::cout << "save_outputs = \t"    << save_outputs << "\n";
        std::cout << "output_dir = \t"       << output_dir << "\n";

    }



    std::cout << "------------------\n";
    //-----------------------------------------
    //Make grid, b-splines & Diffs
    std::cout << "Making Splines...\n";
    std::vector<double> rgrid = make_rgrid(rmin, rmax, ngrid);
    double dr = (rmax-rmin) / (ngrid-1);

    list_of_vecs bsplines       =  generate_splines(nsplines, rgrid);
    list_of_vecs bsplines_diff  =  generate_spline_diffs(nsplines, rgrid);

    std::cout << "Done.\n";

    //-----------------------------------------
    //Make voltage vector
    std::cout << "Generating Voltage potentials \n";

    // Make potential func using greens function
    // Values are hardcoded for this question as we only care about lithium
    double d = 0.2;
    double h = 1.0;
    int Z = 3;

    //Get potentials for both s and l orbitals
    auto Vnuc_s_f = [Z](double r) { return V_hydrogen(r, Z, 0);};
    auto Vnuc_l_f = [Z](double r) { return V_hydrogen(r, Z, 1);};
    auto Vgr_f = [Z, d, h](double r) { return V_green(r, Z, d, h);};

    std::vector<double> Vnuc_s  = vec_from_func(Vnuc_s_f, rgrid);
    std::vector<double> Vnuc_l  = vec_from_func(Vnuc_l_f, rgrid);
    std::vector<double> Vdir     = vec_from_func(Vgr_f, rgrid); //e-e interaction, begin as greens function

    std::vector<double> Vs      = Vnuc_s + Vdir;
    std::vector<double> Vl      = Vnuc_l + Vdir;

    std::cout << "Done.\n";
    std::cout << "----------------------------------\n";

    //-----------------------------------------
    //Calculate first itteration
    std::cout << "Calculating initial energy eigenstates...\n";

    energy_and_waves solutions_s = solve_energies(Vs, bsplines, bsplines_diff, dr);
    energy_and_waves solutions_l = solve_energies(Vl, bsplines, bsplines_diff, dr);

    std::cout << "Done.\n";

    //-----------------------------------------
    //Perform hartree itterations
    std::cout << "----------------------------------\n";
    std::cout << "Doing Hartree Procedure.\n";

    std::vector<energy_and_waves> hartree_sols = hartree_fast(rgrid, Vnuc_s, Vnuc_l,
                                                              solutions_s, solutions_l,
                                                              maxits, tol*10,ens_to_check);

    solutions_s = hartree_sols[0];
    solutions_l = hartree_sols[1];

    std::cout << "Done.\n";
    std::cout << "----------------------------------\n";

    //-----------------------------------------
    //Perform hartree_dock itterations
    std::cout << "----------------------------------\n";
    std::cout << "Doing Hartree-Fock Procedure.\n";

    hartree_sols = hartree_fock_fast(rgrid, Vnuc_s, Vnuc_l,
                                solutions_s, solutions_l,
                                maxits, tol, ens_to_check);

    solutions_s = hartree_sols[0];
    solutions_l = hartree_sols[1];

    std::cout << "Done.\n";
    std::cout << "----------------------------------\n";

    //-----------------------------------------

    std::cout << "\nCalculated energies: \t s orbital \t l orbital \n";
    for (int i=0; i<5; i++){
        std::cout << "Energy Level "<< i+1 << ":\t" << solutions_s.energies[i] << "\t" << solutions_l.energies[i] << "\n";
    }

    std::cout << "\nAverage Positions: \t s orbital \t l orbital \n";
    for (int i=0; i<5; i++){
        std::cout << "Energy Level "<< i+1 << ":\t";
        std::cout << vint(solutions_s.waves[i]*solutions_s.waves[i]*rgrid,dr)<< "\t";
        std::cout << vint(solutions_l.waves[i]*solutions_l.waves[i]*rgrid,dr)<< "\t";
        std::cout << "\n";
    }

    //Get relaxation time
    double R = vint(solutions_l.waves[0] * rgrid * solutions_s.waves[1],dr);
    double omega = fabs(solutions_l.energies[0] - solutions_s.energies[1]);
    double gamma = 2.0/3.0 * pow(R,2) * pow(omega,3) * 1.071 * 10;
    double T = 1/gamma;
    std::cout << "Relaxation time is ~ \t " << T <<" ns \n";
    //-----------------------------------------
    //Outputs and saving

    //Output solutions to .txt files
    if(save_outputs){

        //Save Potential
        std::vector<double> y0_1s1s = YK::ykab(0, solutions_s.waves[0],solutions_s.waves[0],rgrid);
        std::ofstream potfile(output_dir + "potential.txt");
        for (int i=0; i<ngrid; i++){
            potfile << rgrid[i] << "\t ";
            potfile << Vnuc_s[i] << "\t" << Vnuc_l[i] << "\t";
            potfile << Vs[i] << "\t" << Vl[i] << "\t";
            potfile << Vnuc_s[i]+2.0*y0_1s1s[i]<< "\t" << Vnuc_l[i]+2.0*y0_1s1s[i];
            potfile << "\n";
        }
        potfile.close();

        //Save Solutions (Waves & Energies)
        std::ofstream wavefile_s(output_dir + "waves_s.txt");
        std::ofstream wavefile_l(output_dir + "waves_l.txt");
        std::ofstream energyfile(output_dir + "energies.txt");
        std::ofstream avposfile(output_dir + "avpos.txt");
        for (int i = 0; i < ngrid; i++) {

            wavefile_s << rgrid[i] << "\t ";
            wavefile_l << rgrid[i] << "\t ";

            for (int j = 0; j < nsplines; j++) {
                wavefile_s << solutions_s.waves[j][i] << "\t ";
                wavefile_l << solutions_l.waves[j][i] << "\t ";
            }
            wavefile_s << "\n";
            wavefile_l << "\n";
        }
        for (int j = 0; j < nsplines; j++) {
            energyfile << j+1 << "\t";
            energyfile << solutions_s.energies[j] << "\t" << solutions_l.energies[j] << "\n";

            avposfile << j+1 << "\t";
            avposfile << vint(solutions_s.waves[j]*solutions_s.waves[j]*rgrid,dr) << "\t" << vint(solutions_l.waves[j]*solutions_l.waves[j]*rgrid,dr) << "\n";
        }
        energyfile.close();
        wavefile_s.close();
        wavefile_l.close();
        avposfile.close();
    }

    std::cout << "Done.\n";
    //-----------------------------------------
    std::cout << "\t\t Program Done.\n";


    return 0;
}
