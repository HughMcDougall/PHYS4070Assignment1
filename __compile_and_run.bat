echo off
echo "Beginning compile"
echo "	---------------------------------"

:: wsl g++ -o 0QB1 _QB1.cpp matrix.cpp LP_solvers.cpp vector_utils.cpp spline_eval.cpp Potentials_and_Solvers.cpp -llapack -lblas -O3
::  wsl g++ -o 0QB2 _QB2.cpp matrix.cpp LP_solvers.cpp vector_utils.cpp spline_eval.cpp Potentials_and_Solvers.cpp -llapack -lblas -O3
:: wsl g++ -o 0QB3 _QB3.cpp matrix.cpp LP_solvers.cpp vector_utils.cpp spline_eval.cpp Potentials_and_Solvers.cpp -llapack -lblas -O3
:: wsl g++ -o 0QB4 _QB4.cpp matrix.cpp LP_solvers.cpp vector_utils.cpp spline_eval.cpp Potentials_and_Solvers.cpp -llapack -lblas -O3
echo "	---------------------------------"

echo "Compile complete. Running: \n"
echo "	---------------------------------"

:: wsl ./0QB1 60 10001 0.0001 120 1
wsl ./0QB2 60 10001 0.0001 120 1
:: wsl ./0QB3 60 50001 0.00001 120 40 0.000001 5 1
:: wsl ./0QB4 60 50001 0.00001 120 40 0.000001 5 1

echo "	---------------------------------"
pause