#pragma once
//The function for generating the states, and some of the math required to do that

#include<cstdint>
#include<vector>
#include<armadillo>

//It is critically important that I use the same integer type for all binary operations, and seeing that the state is 64 bit, then almost everything else need to be 64 bit as well.
using state_t = uint64_t;

using namespace std;
using namespace arma;

uint64_t binSearch(const vector<state_t>& states, state_t target, uint64_t low=0 ,uint64_t high=-1 /*Intentional underflow, is truncated to the end of the list automatically */);

uint64_t bin_coef(uint64_t n, uint64_t m);


//Print this state (human readable), just as a sanity check
void print_state(state_t S,  int w, int h);//Assumes little Endian states

//Print this state in a format gnuplot can plot
void print_state_data(state_t S,  int w, int h);//Assumes little Endian states

//Print this state in a format gnuplot can plot
void print_superpositions_data(const vector<state_t>& states,vector<double> superpositions , uint64_t& n_states,  int w, int h);//Assumes little Endian states

void print_density_data(ofstream& Out,vec lattice_density,  uint w, uint h);//Assumes little Endian states

//Generate all states and stores them here, assumes that size has been allocated correctly
void generate(vector<state_t>& states,
              uint64_t& N_states,
              uint64_t N_sites,
              uint64_t N_particles);



//Get the density at all sites,
vec get_density(const vector<state_t>& states,vector<double> superpositions , const uint64_t& n_states,  int w, int h);

