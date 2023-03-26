echo "Beginning compile"
echo "	---------------------------------"

g++ -o _QB1 _QB1.cpp matrix.cpp LP_solvers.cpp vector_utils.cpp spline_eval.cpp Potentials_and_Solvers.cpp -llapack -lblas -O3
g++ -o _QB2 _QB2.cpp matrix.cpp LP_solvers.cpp vector_utils.cpp spline_eval.cpp Potentials_and_Solvers.cpp -llapack -lblas -O3
g++ -o _QB3 _QB3.cpp matrix.cpp LP_solvers.cpp vector_utils.cpp spline_eval.cpp Potentials_and_Solvers.cpp -llapack -lblas -O3
g++ -o _QB4 _QB4.cpp matrix.cpp LP_solvers.cpp vector_utils.cpp spline_eval.cpp Potentials_and_Solvers.cpp -llapack -lblas -O3

g++ -o _SingleHydrogen _SingleHydrogen.cpp matrix.cpp LP_solvers.cpp vector_utils.cpp spline_eval.cpp Potentials_and_Solvers.cpp -llapack -lblas -O3
g++ -o _QB3-alt _QB3-alt.cpp matrix.cpp LP_solvers.cpp vector_utils.cpp spline_eval.cpp Potentials_and_Solvers.cpp -llapack -lblas -O3
g++ -o _QB4-alt _QB4-alt.cpp matrix.cpp LP_solvers.cpp vector_utils.cpp spline_eval.cpp Potentials_and_Solvers.cpp -llapack -lblas -O3
echo "	---------------------------------"

echo "Compile complete. \n"
echo "	---------------------------------"

pause