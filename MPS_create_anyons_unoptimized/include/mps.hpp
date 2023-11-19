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

//A matrix product state,
class mps
{
private:

    vector< vector<cx_mat> > sub_matrices;

    vector<uint > D_list;//A list with L numbers with the height the matrices, the height is always the width of the previous

    vector< vector<uint> > subdimension_height;//List of height of submatrices
    vector< vector<uint> > row_start;//Where, in the full matrix, does this subrow start
    vector< uint > n_subrows;//How many subrows are there at this site

    vector< uint > firstsubcol_totalparticles;//Total particles after step j, in the top matrix


    //The top submatrix subrow/subcol, the rest follow diagonally after that (d*L)
    vector<pair<uint,uint> > top_submatrix_location;

    uint d = 2;
    uint64_t L;

    uint N_particles;

    //Is this right, left or mixed canonical, can_center is the center of mixed canonization, if -1 it is all right, if L it is all left, if it is -2 it does not have a canonization center
    int can_center=L;



    //Better, if you know k
    complex<double> get(uint sig,uint j, uint k,uint a,uint ap1) const;
    void set(uint sig,uint j, uint k,uint a,uint ap1,complex<double> data) ;


    //At this j,sig get submatrix k's row and column range (return false if there is none)
    bool sub_matrix_range(uint j, uint sig, uint k, uint& a0, uint& a1, uint& a0p1/*p1 = plus 1, so the next index*/, uint& a1p1) const;



    friend class mpo;
public:


    mps(uint64_t L,uint N_particles,
    const vector<pair<vector<bool> ,complex<double> >
        >& states,uint min_D=1, bool left=true);

    mps(uint64_t L,uint N_particles, uint local_D, bool left=true);


    void print_canonization() const;

    ~mps();

    uint get_d() const {return d;}
    uint64_t get_L() const {return L;}
    int get_can_center() const {return can_center;}

    void DEBUG_PRINT_ALL() const;

    //Step canonization center one place left, from can_center to can_center-1 or can_center+1, also works to canonize mixed states using QR with no compression
    bool step_left_QR();//return true if all right
    bool step_right_QR();//return true of all left


    complex<double> get_coefficient(vector<uint> state) const;
    complex<double> get_coefficient(vector<bool> state) const;

    complex<double> overlap(const mps& that) const;

};
