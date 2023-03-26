echo off
echo "Beginning compile"
echo "	---------------------------------"

wsl g++ -o 0QB1 _QB1.cpp matrix.cpp LP_solvers.cpp vector_utils.cpp spline_eval.cpp Potentials_and_Solvers.cpp -llapack -lblas -O3

echo "Compile complete. Running: \n"
echo "	---------------------------------"

wsl ./0QB1 60 10001 0.00001 120 40 0.000001 2

echo "	---------------------------------"
pause