//
// Created by hughm on 25/03/2023.
//

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

int main(){
    double rmin=1;
    double rmax=101;
    int ngrid = 10;

    std::vector<double> rgrid = make_rgrid(rmin, rmax, ngrid);
    double dr = (rmax-rmin) / (ngrid-1);

    std::vector<double> y = rgrid*rgrid;
    std::vector<double> ydiff1 = vdiff(y);
    std::vector<double> ydiff2 = vdiff(y, dr);

    printv(rgrid);
    printv(y);
    printv(ydiff1);
    printv(ydiff2);

    return 0;
}