#pragma once
#include<cstdint>
#include<iostream>
#include<string>
#include<armadillo>
#include<functional>

using state = uint64_t;
using namespace std; //Namespace containing in and out stream, std::vector (dynamically sized list) and complex numbers
using namespace arma;//Namespace containing matrices and vectors (actual vectors)


vec qnewton( function<double(vec)> func , vec x0, size_t& steps,size_t max_steps=128, double acc=1e-2, bool verbose=false);

//Mirror function to maximize instead
vec max_qnewton( function<double(vec)> func , vec x0, size_t& steps, double acc=1e-2, bool verbose=false);

//Versions without counting number of steps
vec qnewton( function<double(vec)> func , vec x0, double acc=1e-2, bool verbose=false);

//Mirror function to maximize instead
vec max_qnewton( function<double(vec)> func , vec x0, double acc=1e-2, bool verbose=false);


vec downhill_simplex( function<double(vec)> func , vec x0, size_t& steps, double acc=1e-2, bool verbose=false);

//Mirror function to maximize instead
vec max_downhill_simplex( function<double(vec)> func , vec x0, size_t& steps, double acc=1e-2, bool verbose=false);

//Versions without counting number of steps
vec downhill_simplex( function<double(vec)> func , vec x0, double acc=1e-2, bool verbose=false);

//Mirror function to maximize instead
vec max_downhill_simplex( function<double(vec)> func , vec x0, double acc=1e-2, bool verbose=false);


//Taking same input as before, use alglib because I care more about having a well tested and optimized optimization algorithm, than makign everything from scratch. the limited memory Broyden-Fletcher-Goldfarb-Shanno is a quasi newton method which should hopefully make many fewer calls.
//vec libalg_lbfgs( function<double(vec)> func , vec x0 , uint64_t steps, double acc=1e-2, bool verbose=false);//I will keep this third party library contained wihtin minimize.cpp alone, and convert the output back to armadillo style vectors, just to avoid confusion, armadillo is generally superior when it comes to matrix operations, hence why it is still my main library

