#include"generate_states.hpp"
#include"get_state.hpp"

//Has already been done in my header, but I like to include everything explicitly in every file (the compiler will ignore it)
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
#include"minimize.hpp"

using namespace std; //Namespace containing in and out stream, std::vector (dynamically sized list) and complex numbers
using namespace arma;//Namespace containing matrices and vectors (actual vectors)


#include "approx.hpp"

int main(int argc, char* argv[])
{
    //Read arguments
    if ( argc!=6)
    {
        cout<<"Invalid, argument, need "<<argv[0]<<" width height particle_number n0_out n0_out_RAW"<<endl;
        return 1;
    }




    //If this fails, we get 0
    uint64_t nw = atoi(argv[1]);
    uint64_t nh = atoi(argv[2]);

    uint64_t N_sites = nh*nw;


    uint64_t N_particles = atoi(argv[3]);


    if (nw==0)
    {
        cout << "Invalid argument "<<argv[1]<<", must be integer"<<endl;
        return 1;
    }
    if (nh==0)
    {
        cout << "Invalid argument "<<argv[2]<<", must be integer"<<endl;
        return 1;
    }
    if (N_sites >sizeof(state_t)*8)
    {
        cout <<"This program only allows up to "<<(sizeof(state_t)*8)<<" sites, got  "<<nw<<"x"<<nh<<"="<<N_sites<<endl;
        return 1;
    }

    if (N_sites < N_particles)
    {
        cout <<"More particles "<<N_particles<<" than lattice sites  "<<nw<<"x"<<nh<<"="<<N_sites<<endl;
        return 1;
    }


    //Basis states with 3 or 2 particles
    uint64_t N_states0=0;
    vector<state_t> states0(N_states0);

    generate(states0,N_states0,N_sites,N_particles);


    double empty_eigval;
    cx_vec empty_eigvec;
    sp_cx_mat H0;




    cout<<N_sites<<" "<<N_particles<<" "<<N_states0<<" "<<flush;

    {
    get_state(
    N_states0,
    N_sites,
    N_particles,
    0,//N_anyons,
    nw,
    nh,
    states0,
    empty_eigval,
    empty_eigvec,
    H0);

        cout<<"DONE: Energy= "<<empty_eigval<<endl;
    }
    cout<<endl;



    //To calculate the density of particles in the lattice, we need the probability of being in each basis state
    vector<double> prop0_list(N_states0);
    for (uint64_t j = 0; j < N_states0; ++j)
    {
        prop0_list[j]= norm(empty_eigvec[j]);
    }

    vec lattice_dens0= get_density(states0,prop0_list, N_states0, nw, nh);
    //print_density_data(lattice_dens0,  nw, nh);//Ok now print it so gnuplot can plot it
    //cerr << endl;


    ofstream n_file(argv[4]);

    print_density_data(n_file,lattice_dens0 , nw, nh);
    n_file.close();

    ofstream raw_file(argv[5]);

    raw_file<<empty_eigval<<endl;

    for (uint64_t j = 0; j < N_sites; ++j)
    {
        raw_file<<lattice_dens0(j)<<endl;
    }

    raw_file.close();

    return 0;
}
