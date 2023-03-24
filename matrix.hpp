//
// header file for matrix class. Adapted from workshop 2 examples
//

#ifndef ASSIGNMENT1_MATRIX_H
#define ASSIGNMENT1_MATRIX_H

#pragma once
#include <cassert>
#include <vector>
#include <iostream>

using unt = std::size_t;
namespace matrix{

    class sqmatrix{

        private:
            std::vector<double> _datavec;
            unt _N;
            unt _size;

        public:
            //Init
            sqmatrix(int N): _N(N), _size(N*N){_datavec.resize(_size);}

            //Gets
            double& at(unt i, unt j)        {return _datavec[i * _N + j];}            //For assigning
            double  at(unt i, unt j) const  {return _datavec[i * _N + j];}            //Const version
            unt N() const {return _N;}                                                 //Shield for _N
            unt size() const {return _size;}                                           //Shield for _size=_N*_N
            double *data() {return _datavec.data();}                                    //Shortcut to data array
            double &operator()(std::size_t i, std::size_t j) { return at(i, j); }       //Quick format of .at
            double operator()(std::size_t i, std::size_t j) const { return at(i, j); }  //Const Version

            //Misc Utilities
            void print(){
                for (unt i=0; i < _N; i++){
                    for (unt j=0; j < _N; j++){
                        std::cout << _datavec[i*_N+j]<< "\t";
                        }
                    std::cout << "\n";
                    }
                }

            //Trace
            double trace(){
                /// Returns the matrix's trace, the sum of its diagonal elements
                double out=0;
                for (int i=0;i<_N;i++){out+=_datavec[i * _N + i];}
                return(out);
            }

            //Diagonal
            std::vector<double> diag(){
                /// Returns the diagonal elements of the matrix in a vector
                std::vector<double> out(_N);
                for (int i=0;i<_N;i++){out[i]=_datavec[i * _N + i];}
                return(out);
            }

        };//sqmatrix
    // -----------------------------
    //Declare overload operators. Defined properly in cpp file
    sqmatrix operator+=(sqmatrix &a, const sqmatrix &b);
    sqmatrix operator+=(sqmatrix &a, const double &b);
    sqmatrix operator-=(sqmatrix &a, const sqmatrix &b);
    sqmatrix operator-=(sqmatrix &a, const double &b);
    sqmatrix operator*=(sqmatrix &a, const sqmatrix &b);
    sqmatrix operator*=(sqmatrix &a, const double &b);
    sqmatrix operator/=(sqmatrix &a, const sqmatrix &b);
    sqmatrix operator/=(sqmatrix &a, const double &b);

    sqmatrix operator+(sqmatrix a, const sqmatrix &b);
    sqmatrix operator+(sqmatrix a, const double &b);
    sqmatrix operator+(const double &b, sqmatrix a);
    sqmatrix operator-(sqmatrix a, const sqmatrix &b);
    sqmatrix operator-(sqmatrix a, const double &b);
    sqmatrix operator*(sqmatrix a, const sqmatrix &b);
    sqmatrix operator*(sqmatrix a, const double &b);
    sqmatrix operator*(const double &b, sqmatrix a);
    sqmatrix operator/(sqmatrix a, const sqmatrix &b);
    sqmatrix operator/(sqmatrix a, const double &b);

    } //Namespace

#endif //ASSIGNMENT1_MATRIX_H
