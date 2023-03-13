//
// Created by hughm on 11/03/2023.
// Differential and integral functions
//

#include "diff_int.h"
#include <iostream>
#include <functional>

namespace diff_int{
    testclass operator+(testclass a, testclass b){
        testclass out;
        out.a=a.a+b.a;
    };
}