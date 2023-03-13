cmake_minimum_required(VERSION 3.24)
project(Assignment1)

set(CMAKE_CXX_STANDARD 17)

add_executable(Assignment1 main.cpp matrix.cpp matrix.hpp diff_int.cpp diff_int.hpp LP_solvers.cpp LP_solvers.hpp vector_utils.cpp vector_utils.h)
