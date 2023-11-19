#pragma once

#include<cstdint>

//NOTE! std::vector is a dynamic sized list, it is not a vector in the linear algebra sense of the word, that will be the class vec from the armadillo library
#include<vector>
#include<iostream>
#include<string>
#include<complex>
#include<exception>
#include<iomanip>
#include<fstream>

//The actual math library
#include <armadillo>
#include <chrono>

using namespace std; //Namespace containing in and out stream, std::vector (dynamically sized list) and complex numbers
using namespace arma;//Namespace containing matrices and vectors (actual vectors)


#include "mps.hpp"

//A matrix product operator,
class mpo
{
private:

    vector<sp_cx_mat> matrices;//There are L*(d**2) matrices


/*
W have both sig and sig1 and are stored

W^{l=0,sig=0,sig1=0}
W^{l=0,sig=0,sig1=1}
W^{l=0,sig=1,sig1=0}
W^{l=0,sig=1,sig1=1}
*/


    uint Dw;//matrix dimension, I assume the first matrix is 1*Dw, the last is Dw*1, the rest are Dw

    uint d = 2;
    uint64_t L;
    uint dTwo=4;





    bool is_n=false;//This is actually the density operator n, in that case we will let n[l] "jump" between sites to test them all, without making 64 or more different operators


    //Get the L[l] tensor, calculated given the tensor L[l-1]. L[l] has indices a0_l, a1_l and b_l, where  a_l are the indices of the left matrices of the MPO
    cx_cube getL_recursive(const mps& mps1, uint l,const cx_cube& L_prev) const;
    //Same for R[l]
    cx_cube getR_recursive(const mps& mps1, uint l,const cx_cube& R_next) const;

    //For all rows (except the last one) what singular column do we have an entry in
    vector<uint> column_lookup;


    //Get submatrix indices for converting reduced H indices into what submatrix (sig,k) and index (a,ap1) we put this entry into
    struct submatrix_index
    {
        uint sig;
        uint k  ;
        uint a  ;//Specifically a submatrix coordinate
        uint ap1;
        uint a_min;
        uint ap1_min;
submatrix_index(uint s, uint _k, uint _a, uint _ap1, uint _a_min, uint _ap1_min)
{sig = s; k = _k; a=_a; ap1=_ap1;a_min=_a_min;ap1_min=_ap1_min;}
    };


    static void make_H_thread(uint I_min, uint I_max, uint reduced_size,uint mW, uint nW,uint n_c,const vector<submatrix_index>& submatrix_indices, const cx_cube& tensorL, const cx_cube& tensorR, sp_cx_mat* grand_H,uint d, uint dTwo, const vector<uint>& column_lookup, const vector<sp_cx_mat>& matrices);


public:
    //Generate A matrices of a matrix product state
    mpo(uint64_t L,uint d,vector< vector<complex<double> > > Matrices);


    //n operator
    mpo(
        uint64_t Lx,
        uint64_t Ly);


    //1D hopping
    mpo(
        uint64_t _L,
        const vector<double>& V,
        complex<double> t);

    mpo(
        uint64_t Lx,
        uint64_t Ly,
        vector<double>& V,
        vector<complex<double> >& tx,
        vector<complex<double> >& ty
        );


    void set_potential(const vector<double>& V);

    mps find_ground(const mps& that, double& energy, uint max_sweeps=4);


    complex<double> get_expectation_value(const mps& mps1, uint j = -1/*ONLY FOR expectation value of n at a particular j*/, bool verbose=false);

    void DEBUG_PRINT_ALL() const;
};
