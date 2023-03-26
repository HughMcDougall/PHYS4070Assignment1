//
// Grab-bag of general utility functions for working with std::vectors
// Over-loads basic arithmetic (+, -, *, /) for vectors and doubles / ints
// HM Mar 23
//

#include <iostream>
#include <vector>
#include <cassert>
#include <functional>

using function_1D  = std::function<double(double)>;
using list_of_vecs = std::vector<std::vector<double>>;

//=======================================================
//Overload vector operations to make direct products easier
//V-V Multiplication
std::vector<double> operator*=(std::vector<double> & a, const std::vector<double> & b){
    assert(a.size()==b.size() && "Tried to multiply two vectors with different lengths");
    for (int i=0; i<a.size(); i++){
        a[i]*=b[i];
    }
    return(a);
}
std::vector<double> operator*(std::vector<double> a, const std::vector<double> & b){return a*=b;}

//V-V Division
std::vector<double> operator/=(std::vector<double> & a, const std::vector<double> & b){
    assert(a.size()==b.size() && "Tried to divide two vectors with different lengths");
    for (int i=0; i<a.size(); i++){
        a[i]/=b[i];
    }
    return(a);
}
std::vector<double> operator/(std::vector<double> a, const std::vector<double> & b){return a/=b;}

//V-V Addition
std::vector<double> operator+=(std::vector<double> & a, const std::vector<double> & b){
    assert(a.size()==b.size() && "Tried to add two vectors with different lengths");
    for (int i=0; i<a.size(); i++){
        a[i]+=b[i];
    }
    return(a);
}
std::vector<double> operator+(std::vector<double> a, const std::vector<double> & b){return a+=b;}

//V-V Subtraction
std::vector<double> operator-=(std::vector<double> & a, const std::vector<double> & b){
    assert(a.size()==b.size() && "Tried to subtract two vectors with different lengths");
    for (int i=0; i<a.size(); i++){
        a[i]-=b[i];
    }
    return(a);
}
std::vector<double> operator-(std::vector<double> a, const std::vector<double> & b){return a-=b;}

//V-D Multiplication
std::vector<double> operator*=(std::vector<double> & v, const double & a){
    for (int i=0; i<v.size(); i++){
        v[i]*=a;
    }
    return(v);
}
std::vector<double> operator*(std::vector<double> v, const double & a){return(v*=a);}
std::vector<double> operator*(const double & a, std::vector<double> v){return(v*a);}

//V-D Division
std::vector<double> operator/=(std::vector<double> & v, const double & a){
    for (int i=0; i<v.size(); i++){
        v[i]/=a;
    }
    return(v);
}
std::vector<double> operator/(std::vector<double> v, const double & a){return(v/=a);}
std::vector<double> operator/(const double & a, std::vector<double> v){return(v/a);}

//V-D Addition
std::vector<double> operator+=(std::vector<double> & v, const double & a){
    for (int i=0; i<v.size(); i++){
        v[i]+=a;
    }
    return(v);
}
std::vector<double> operator+(std::vector<double> v, const double & a){return(v+=a);}
std::vector<double> operator+(const double & a, std::vector<double> v){return(v+a);}

//V-D Subtraction
std::vector<double> operator-=(std::vector<double> & v, const double & a){return(v+=a*-1);}
std::vector<double> operator-(std::vector<double> v, const double & a){return(v-=a);}
std::vector<double> operator-(const double & a, std::vector<double> v){return(v-a);}

//=======================================================
//Conversions & Integrations
//Function to integrate over a vector
double vint(const std::vector<double>& a, double dx){
    /// Integrates over a vector with constant spacing. Uses simspsons rule
    int n = a.size();
    if (n%2==0){std::cerr<<"Warning! Integrating vector with even number of points: \t" << n;}
    double out = 0;

    //Peform weighted summation in keeping /w simpsons rule
    out += a.front();
    for (int i=1; i<n; i+=2){ //Odd Indices
        out+=a[i]*4;
    }
    for (int i=2; i<n; i+=2){ //Even Indices
        out+=a[i]*2;
    }
    out += a.back();

    out /= 3;

    // If differental element has been provided, use it
    if (dx!=0){out*=dx;}

    return out;
}

//Function to integrate over a vector
std::vector<double> vdiff(const std::vector<double>& a, double dx){
    /// Differentiates a vector with constant spacing

    int n = a.size();
    std::vector<double> out(n);

    out[0] = a[1]-a[0];             //forward difference
    for (int i=1; i<n-1; i++){
        out[i]=(a[i+1]-a[i-1])/2;       //central difference
    }
    out[n-1] = a[n-1]-a[n-2];   //backward difference

    // If differental element has been provided, use it
    if (dx!=0){out/=dx;}

    return out;
}

std::vector<double> make_rgrid(double rmin = 0.001, double rmax = 100, int n_grid = 101){
    /// Creates a linspace of radial grid points. Saves time on passing the same args to every function
    double dr = (rmax-rmin) / (n_grid-1);
    double r = rmin;
    std::vector<double> out(n_grid);

    for (int i=0; i<n_grid; i++){
        out[i]=r;
        r+=dr;
    }

    return out;
}

std::vector<double> vec_from_func(const function_1D& V, const std::vector<double>& rgrid){
    ///Converts a 1D potential function into a std::vector for quick-swapping of potential functions
    int n_grid = rgrid.size();
    std::vector<double> out (n_grid);

    for (int i=0; i<n_grid; i++){
        out[i]=V(rgrid[i]);
    }

    return out;
}

//=======================================================
//Utility Functions
void printv(const std::vector<double>& a){
    /// Prints all elements in a 1D vector of doubles to std::cout. Good for debugging
    for (int i=0; i<a.size(); i++){
        std::cout << a[i] << ", ";
    }
    std::cout<<"\n";
}

std::vector<double> vcopy(std::vector<double> a){
    /// Returns a copy of vector to assign a second version of it.
    return(a);
}

std::vector<double> zeros_like(std::vector<double> a){
    /// Returns a copy of vector to assign a second version of it.
    return(a*0);
}
