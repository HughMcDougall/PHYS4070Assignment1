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
#include "vector_utils.hpp"
#include "spline_eval.hpp"
#include "Potentials_And_Solvers.hpp"

//=======================================================
//Aliases & Holder Structs
using function_1D  = std::function<double(double)>;
using list_of_vecs = std::vector<std::vector<double>>;

//=======================================================
//Main runtime
int main(int argc, char *argv[]) {
    /// Read Inputs
    /// Order of inputs : Z,    l,  nsplines,   ngrid,  rmin,   rmax,   save_outputs,   output_dir
    /// Default inputs  : 1,    0,  60,         1001,   0.001,  20,     True,           "./outputs_B1"

    std::cout << "Program Starting\n";
    int Z, l;
    int nsplines, ngrid;
    double rmin, rmax;
    bool save_outputs;
    std::string output_dir;

    // Read inputs. Scoped to allow to collapse easily in editor
    {
        if (argc > 1){
            Z = std::stoi(argv[1]);
        } else{
            Z = 1;
        }

        if (argc > 2){
            l = std::stoi(argv[2]);
        } else{
            l = 0;
        }

        if (argc > 3){
            nsplines = std::stoi(argv[3]);
        } else{
            nsplines = 60;
        }

        if (argc > 4){
            ngrid = std::stoi(argv[4]);
        } else{
            ngrid = 1001;
        }

        if (argc > 5){
            rmin = std::stof(argv[5]);
        } else{
            rmin = 0.001;
        }

        if (argc > 6){
            rmax = std::stof(argv[6]);
        } else{
            rmax  = 20;
        }

        if (argc > 7){
            save_outputs = std::stoi(argv[7]);
        } else{
            save_outputs  = false;
        }

        if (argc > 8){
            output_dir = argv[8];
        } else{
            output_dir  = "./outputs/B1/";
        }
    }


    //-----------------------------------------
    //Make grid, b-splines & Diffs
    std::cout << "Parameters:\n";
    std::cout << "Z = \t"           << Z << "\n";
    std::cout << "l = \t"           << l << "\n";
    std::cout << "nsplines = \t"    << nsplines << "\n";
    std::cout << "ngrid = \t"       << ngrid << "\n";
    std::cout << "rmin = \t"    << rmin << "\n";
    std::cout << "rmax = \t"       << rmax << "\n";
    std::cout << "save_outputs = \t"    << save_outputs << "\n";
    std::cout << "output_dir = \t"       << output_dir << "\n";

    std::vector<double> rgrid = make_rgrid(rmin, rmax, ngrid);

    list_of_vecs bsplines       =  generate_splines(nsplines, rgrid);
    list_of_vecs bsplines_diff  =  generate_spline_diffs(nsplines, rgrid);

    std::cout << "Done.\n";

    //-----------------------------------------
    //Make voltage vector
    std::cout << "Generating Voltage potential for Z= " << Z << " and l= " << l << "\n";

    // Make potential func
    double d = 0.2;
    double h = 0.2;
    function_1D  Vfunc;
    Vfunc = [Z, l](double r) { return V_hydrogen(r, Z, l); };
    std::vector<double> V = vec_from_func(Vfunc, rgrid);

    std::cout << "Done.\n";

    //-----------------------------------------
    //Actually perform calculations
    std::cout << "Calculating energy eigenstates...\n";

    energy_and_waves energy_sols = solve_energies(V, bsplines, bsplines_diff);

    //-----------------------------------------
    //Outputs and saving

    std::cout << "Calculated energies:\n";
    for (int i=0; i<5; i++){
        std::cout << "Energy Level "<< i+1 << ":\t" << energy_sols.energies[i] << "\n";
    }

    std::cout << "Average Positions:\n";
    for (int i=0; i<5; i++){
        std::cout << "Energy Level "<< i+1 << ":\t" << vint(energy_sols.waves[i]*energy_sols.waves[i]*rgrid) /  vint(energy_sols.waves[i]*energy_sols.waves[i]) << "\n";
    }

    //Output solutions to .txt files
    if(save_outputs){

            //Save Splines
            std::ofstream spline_file(output_dir+"splines.txt");
            std::ofstream spline_diffile(output_dir+"spline_diffs.txt");
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

            //Save Potential
            std::ofstream potfile(output_dir + "potential.txt");
            for (int i=0; i<ngrid; i++){
                potfile << rgrid[i] << "\t ";
                potfile << V[i];
                potfile << "\n";
            }
            potfile.close();

            //Save Solutions (Waves & Energies)
            std::ofstream wavefile(output_dir + "waves.txt");
            std::ofstream energyfile(output_dir + "energies.txt");
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


    return 0;
}
