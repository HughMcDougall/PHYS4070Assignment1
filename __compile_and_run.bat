echo off
echo "Beginning compile"
echo "	---------------------------------"

wsl g++ -o 0QB1 _QB1.cpp matrix.cpp LP_solvers.cpp vector_utils.cpp spline_eval.cpp Potentials_and_Solvers.cpp -llapack -lblas -O3
wsl g++ -o 0QB2 _QB2.cpp matrix.cpp LP_solvers.cpp vector_utils.cpp spline_eval.cpp Potentials_and_Solvers.cpp -llapack -lblas -O3

echo "Compile complete. Running: \n"
echo "	---------------------------------"

wsl ./0QB2

echo "	---------------------------------"
pause