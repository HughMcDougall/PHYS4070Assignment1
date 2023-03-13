//
// Created by hughm on 11/03/2023.
// Differential and integral functions
//

#include <iostream>
#include <functional>

#ifndef ASSIGNMENT1_DIFF_INT_H
#define ASSIGNMENT1_DIFF_INT_H

double diff(std::function<double(double)> f, double dx);
double integ(std::function<double(double)> f, float a, float b, double dx);
#endif //ASSIGNMENT1_DIFF_INT_H

namespace diff_int{
    class testclass{
    public:
        void go(){std::cout << "hello world \n";}
        int a = 1;
    };

    testclass operator+(testclass a, testclass b);
}