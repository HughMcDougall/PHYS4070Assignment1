echo off
echo "Beginning compile"
echo "	---------------------------------"

wsl g++ -o B1 B1.cpp matrix.cpp LP_solvers.cpp vector_utils.cpp spline_eval.cpp Potentials_and_Solvers.cpp -llapack -lblas -O3

echo "Compile complete. Running: \n"
echo "	---------------------------------"

wsl ./B1

echo "	---------------------------------"
pause