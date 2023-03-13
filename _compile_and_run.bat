echo off
echo "Beginning compile"
echo "	---------------------------------"

wsl g++ -o out main.cpp matrix.cpp LP_solvers.cpp -llapack -lblas

echo "Compile complete. Running: \n"
echo "	---------------------------------"

wsl ./out

echo "	---------------------------------"
pause