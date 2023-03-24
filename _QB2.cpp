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
    /// Order of inputs : nsplines,   ngrid,   rmin,   rmax,   save_outputs,   output_dir
    /// Default inputs  : 60,         5001,   0.001,  100,     False,           "./outputs_B2"

    std::cout << "Program Starting\n";
    int nsplines, ngrid;
    double rmin, rmax;
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
            ngrid = 5001;
        }

        if (argc > 3){
            rmin = std::stof(argv[3]);
        } else{
            rmin = 0.001;
        }

        if (argc > 4){
            rmax = std::stof(argv[4]);
        } else{
            rmax  = 100;
        }

        if (argc > 5){
            save_outputs = std::stoi(argv[5]);
        } else{
            save_outputs  = false;
        }

        if (argc > 6){
            output_dir = argv[6];
        } else{
            output_dir  = "./outputs/B2/";
        }


    }


    //-----------------------------------------
    //Make grid, b-splines & Diffs
    std::cout << "Parameters:\n";
    std::cout << "nsplines = \t"    << nsplines << "\n";
    std::cout << "ngrid = \t"       << ngrid << "\n";
    std::cout << "rmin = \t"    << rmin << "\n";
    std::cout << "rmax = \t"       << rmax << "\n";
    std::cout << "save_outputs = \t"    << save_outputs << "\n";
    std::cout << "output_dir = \t"       << output_dir << "\n";

    std::cout << "------------------\n";

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
    std::vector<double> Vgr     = vec_from_func(Vgr_f, rgrid);

    std::vector<double> Vs      = Vnuc_s + Vgr;
    std::vector<double> Vl      = Vnuc_l + Vgr;

    std::cout << "Done.\n";
    std::cout << "----------------------------------\n";

    //-----------------------------------------
    //Calculate first itteration
    std::cout << "Calculating initial energy eigenstates...\n";

    energy_and_waves solutions_s = solve_energies(Vs, bsplines, bsplines_diff, dr);
    energy_and_waves solutions_l = solve_energies(Vl, bsplines, bsplines_diff, dr);

    std::cout << "Done.\n";

    std::cout << "Calculated energies: \t s orbital \t l orbital \n";
    for (int i=0; i<5; i++){
        std::cout << "Energy Level "<< i+1 << ":\t" << solutions_s.energies[i] << "\t" << solutions_l.energies[i] << "\n";
    }

    std::cout << "Average Positions: \t s orbital \t l orbital \n";
    for (int i=0; i<5; i++){
        std::cout << "Energy Level "<< i+1 << ":\t";
        std::cout << vint(solutions_s.waves[i]*solutions_s.waves[i]*rgrid,dr)<< "\t";
        std::cout << vint(solutions_l.waves[i]*solutions_l.waves[i]*rgrid,dr)<< "\t";
        std::cout << "\n";
    }

    //-----------------------------------------
    //Calculate first order adjustments
    std::cout << "----------------------------------\n";
    std::cout << "Calculating first order energy perturbations...\n";
    //Get y^{0}_{1s1s}
    std::vector<double> y0_1s1s = YK::ykab(0, solutions_s.waves[0],solutions_s.waves[0],rgrid);
    std::vector<double> adjust = 2 * y0_1s1s - Vgr; //No need to calculate this twice

    for (int i=0; i<nsplines; i++){
        solutions_s.energies[i]+=vint(solutions_s.waves[i]*solutions_s.waves[i]*adjust, dr);
        solutions_l.energies[i]+=vint(solutions_l.waves[i]*solutions_l.waves[i]*adjust, dr);
    }

    std::cout << "Done.\n";

    std::cout << "Calculated energies: \t s orbital \t l orbital \n";
    for (int i=0; i<5; i++){
        std::cout << "Energy Level "<< i+1 << ":\t" << solutions_s.energies[i] << "\t" << solutions_l.energies[i] << "\n";
    }

    //-----------------------------------------
    //Outputs and saving

    //Output solutions to .txt files
    if(save_outputs){

        //Save Potential
        std::ofstream potfile(output_dir + "potential.txt");
        for (int i=0; i<ngrid; i++){
            potfile << rgrid[i] << "\t ";
            potfile << Vs[i] << "\t" << Vl[i] << "\t" << Vs[i]+adjust[i]<< "\t" << Vl[i]+adjust[i];
            potfile << "\n";
        }
        potfile.close();

        //Save Solutions (Waves & Energies)
        std::ofstream wavefile_s(output_dir + "waves_s.txt");
        std::ofstream wavefile_l(output_dir + "waves_l.txt");
        std::ofstream energyfile(output_dir + "energies.txt");
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
        }
        energyfile.close();
        wavefile_s.close();
        wavefile_l.close();
    }

    std::cout << "Done.\n";
    //-----------------------------------------
    std::cout << "\t\t Program Done.\n";


    return 0;
}
