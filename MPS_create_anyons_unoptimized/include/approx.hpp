#pragma once

#include<complex>

bool approx(double a, double b,double tau = 1e-9, double epsilon =1e-9);


bool approx(std::complex<double> a, std::complex<double> b,double tau = 1e-9, double epsilon =1e-9);
