echo off
echo "Beginning compile"
echo "	---------------------------------"

wsl g++ -o 0QB1 _QB1.cpp matrix.cpp LP_solvers.cpp vector_utils.cpp spline_eval.cpp Potentials_and_Solvers.cpp -llapack -lblas -O3
wsl g++ -o 0QB2 _QB2.cpp matrix.cpp LP_solvers.cpp vector_utils.cpp spline_eval.cpp Potentials_and_Solvers.cpp -llapack -lblas -O3
wsl g++ -o 0QB3 _QB3.cpp matrix.cpp LP_solvers.cpp vector_utils.cpp spline_eval.cpp Potentials_and_Solvers.cpp -llapack -lblas -O3

echo "Compile complete. Running: \n"
echo "	---------------------------------"

wsl ./0QB1 1 0 60 12001 0.0001 120

echo "	---------------------------------"
pause