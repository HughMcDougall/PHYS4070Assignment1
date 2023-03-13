#include "bspline.hpp"
#include <iostream>

int main() { //

  // double dr = (rmax - r0) / (nsteps - 1);

  double r0 = 1.0e-4;
  double rmax = 50.0;
  int k = 7;
  int n = 30;

  // All splines go to zero at r=0, excapt the 0th spline (index 0)
  // The 0th spline (index 0) is only non-zero below r0
  // The (k)th spline (and above) [index k] are non-zero only above r0
  // All splines go to zero at rmax, except for the last one (n-1)th spline
  BSpline bspl(k, n, r0, rmax, false);

  std::cout
      << "B-spline example:\n\n"
         "bspl.b(i,r); //gives ith spline at r=r (i in range i=[0,N-1])\n"
         "Note: Each B-spline is zero over most of the domain\n"
         "b(i,r) is non-zero only if t_i <= x < t_{i+n}, where t_i are knots\n"
         "Knots are splaced logarithmically (by default) over [r0,rmax]\n\n";

  // evaluate all the B-splines at 10 points between (r0, rmax):
  int num_steps = 10;
  double dr = (rmax - r0) / (num_steps - 1);

  for (int i_spl = 0; i_spl < n; ++i_spl) {
    std::cout << i_spl << ", ";
    for (int i_r = 0; i_r < num_steps; ++i_r) {
      double r = r0 + dr * i_r;
      // std::cout << bspl.b(i_spl, r) << ", ";
      printf("%.2e, ", bspl.b(i_spl, r));
    }
    std::cout << "\n";
  }

  std::cout << "done\n";

  return 0;
}
