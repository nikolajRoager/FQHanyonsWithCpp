#include "approx.hpp"

#include<vector>
#include<iostream>
#include<string>
#include<complex>
#include<exception>

//The actual math library
#include <armadillo>



using namespace std;



bool approx(std::complex<double> a, std::complex<double> b,double tau, double epsilon)
{
    return approx(a.real(),b.real(),tau,epsilon) && approx(a.imag(),b.imag(),tau,epsilon);
}

bool approx(double a, double b,double tau, double epsilon)
{

    if (abs(a-b)<=tau)
        return true;
    else if (abs(a-b)/(abs(a)+abs(b))<=epsilon)//No risk of divide by 0, if the user asks for a=b=0, then the above always returns true (unless some idiot asks for a negative tolerance and then asks if 0=0)
        return true;
    else
        return false;
}

