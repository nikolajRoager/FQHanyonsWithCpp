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
    if (argc<8)
    {
        cout<<"Invalid, argument, need "<<argv[0]<<" width height q_output_file q_all_file N_particles position_file (repeat with as many V as there are key-frames)"<<endl;
        return 1;
    }

    //If this fails, we get 0
    uint64_t nw = atoi(argv[1]);
    uint64_t nh = atoi(argv[2]);

    uint64_t N_sites = nh*nw;


    uint64_t N_particles = atoi(argv[5]);


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

    if (0 > N_particles)
    {
        cout<<"Negative particle number "<< N_particles<<" not allowed"<<endl;
        return 1;
    }



    //I am very unlikely to get a number of states which needs uint64_t here, but I do need to use the same type as I use for the states (written as a 64 bit number), some of the binary operations may fail otherwise

    //Basis states with 3 or 2 particles
    uint64_t N_states0=0;
    vector<state_t> states0(N_states0);

    uint64_t N_states1=0;
    vector<state_t> states1(N_states1);

    generate(states0,N_states0,N_sites,N_particles);
    generate(states1,N_states1,N_sites,N_particles-1);


    double empty_eigval;
    cx_vec empty_eigvec;
    sp_cx_mat H0;

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

    //To calculate the density of particles in the lattice, we need the probability of being in each basis state
    vector<double> prop0_list(N_states0);
    for (uint64_t j = 0; j < N_states0; ++j)
    {
        prop0_list[j]= norm(empty_eigvec[j]);
    }

    vec lattice_dens0= get_density(states0,prop0_list, N_states0, nw, nh);
    //print_density_data(lattice_dens0,  nw, nh);//Ok now print it so gnuplot can plot it
    //cerr << endl;


    //Now add in two anyons, removing one particle (so keeping phi the same) and adding a potential

    //But first. Set up all the variables we would like to save, and need to use

    //Traping potentials AND auxillary potential on all other sites put together
    vec V(N_sites);

    vector<vec> potentials(argc-6,vec(N_sites));
    vector<vector<uint> > anyon_sites(argc-6,vector<uint>(N_sites));

    for (uint j = 0; j < argc-6; ++j)
    {
        ifstream pos_file(argv[6+j]);

        if (!pos_file.is_open())
        {
            cout<<"Could not open potential file "<<argv[6+j]<<endl;
            return 1;
        }

        string potential_name;
        string operation;
        pos_file>>potential_name;
        pos_file>>operation;

        ifstream potential_file(potential_name);

        for (uint64_t i = 0; i < N_sites; ++i)
        {
            if (! (potential_file>>potentials[j](i)))
            {
                cout<<potential_name<<" did not match number of sites (site "<<i<<" could not be read)" <<endl;
                return 1;
            }
            uint64_t I =0;
            if (! (pos_file>>I))
            {
                cout<<argv[3]<<" did not match number of sites (site "<<i<<" could not be read)" <<endl;
                return 1;
            }

            if (I>2)
                cout<<"Site "<<i<<" could not be read in "<<argv[3]<<" is outside range, should be 0 (no anyons), 1 or 2 (either anyon), got "<<I<<endl;

            anyon_sites[j][i] = I;

        }
        potential_file.close();


        if (operation.compare("mirror_x")==0)
        {
            for (uint y = 0; y < nh; ++y)
            {
                for (uint x = 0; x < nw/2; ++x)
                {
                    double temp = potentials[j](x+y*nw);
                    potentials[j](x+y*nw)=potentials[j](nw-1-x+(y)*nw);
                    potentials[j](nw-1-x+(y)*nw)=temp;
                }
                cout<<endl;
            }


        }
        else if (operation.compare("rotate")==0)
        {
            vec new_potential (N_sites);
            for (uint x = 0; x < nw; ++x)
                for (uint y = 0; y < nh; ++y)
                {
                    double temp = potentials[j](x+y*nw);
                    new_potential (x+y*nw)=potentials[j](nw-1-y+x*nw);
                    new_potential (nw-1-y+x*nw)=temp;
                }
            potentials[j]=new_potential;
        }
        else if (operation.compare("mirror_rotate")==0)
        {
            vec new_potential (N_sites);
            for (uint x = 0; x < nw; ++x)
                for (uint y = 0; y < nh; ++y)
                {
                    double temp = potentials[j](x+y*nw);
                    new_potential (x+y*nw)=potentials[j](y+x*nw);
                    new_potential (y+x*nw)=temp;
                }
            potentials[j]=new_potential;
        }


        pos_file.close();

    }

    cout<<"Potentials"<<endl;
    for (uint j = 0; j < argc-6; ++j)
    {
        cout<<j<<endl;
            for (uint y = 0; y < nh; ++y)
            {
                for (uint x = 0; x < nw; ++x)
                {
                    cout<<potentials[j][x+y*nw]<<' ';
                }
                cout<<endl;
            }

    }

/*
    for (uint j = 0; j < argc-6; ++j)
    {
    cout<<j<<endl;
        for (uint y = 0; y < nh; ++y)
        {
            for (uint x = 0; x < nw; ++x)
            {
                cout<<anyon_sites[j][x+y*nw]<<' ';
            }
            cout<<endl;
        }

    cout<<endl;
    }
*/

    vec lattice_dens1;//Will be written to, if we want to plot the actual density

    //The actual density difference we found
    vec q(N_sites);

    get_H0(
        N_states1,
        N_sites,
        N_particles-1,
        2,//N_anyons,
        nw,
        nh,
        states1,
        H0
    );

    //I will use a lambda function, to make getting the fitness, densities and charge difference easier, though C++ allows me to modify the captures, the things actually written to are all given as arguments
    auto get_q=[&lattice_dens0,N_states1, &states1,N_sites,N_particles,&H0,nw,nh](vec V, vec& lattice_dens1, vec& q, cx_vec& eigvec, double& Energy) -> void
    {


        get_state(
            N_states1,
            N_sites,
            N_particles-1,
            2,//N_anyons,
            nw,
            nh,
            states1,
            Energy,  //written to
            eigvec,   //written to
            H0,
            V
        );

        //Get the propability of finding the system in either basis state, this is one step before finding the expected number of particles in each site

        vector<double> prop1_list(N_states1);
        //I strictly don't need to use uint64 here, but the compiler just gets real mad if I use different types
        for (uint64_t j = 0; j < N_states1; ++j)
        {
            prop1_list[j]= norm(eigvec[j]);
        }

        lattice_dens1 = get_density(states1,prop1_list, N_states1, nw, nh);


        //Check
        q=vec(N_sites);
        double fitness=0;

        double anyon0_charge=0;
        double anyon1_charge=0;


        for (uint64_t i = 0; i < N_sites; ++i)
        {
                double n0 =lattice_dens0(i);
                double n1 =lattice_dens1(i);
                q(i)=n0-n1;
        }
    };



    size_t steps;
    cerr<<std::setprecision(16);
    cout<<std::setprecision(4);

    //Use input potential
    ofstream q_file(argv[3]);

    if (!q_file.is_open())
    {
        cout<<"Could not open q output "<<argv[3]<<endl;
        return 1;
    }

    ofstream q_all_file(argv[4]);

    if (!q_all_file.is_open())
    {
        cout<<"Could not open q output "<<argv[4]<<endl;
        return 1;
    }

    cx_vec first_eigvec;
    cx_vec p_eigvec;
    cx_vec   eigvec;
    for (uint i = 0; i < argc-6-1; ++i)
    {


        cout<<"Saving "<<i<<endl;
        float default_step_size = 0.02;


        float stepsize = default_step_size ;
        uint zoom_level =1;

        bool emergency_redo =false;
        float prev_fac=0;

        uint steps_since_emergency=0;

        bool fin_final = false;
        uint steps_here=0;

        for (float fac = 0; fac< 1.0 || (!fin_final && i==argc-8 ); fac+= stepsize)
        {
            ++steps_here;
            if (fac>=1.0)
            {
                fac=1.0;
                fin_final=true;//An extra final step for the very last one
                cout<<"FINAL STEP"<<endl;
            }


            //float fac = float(j)/Max;
            cout<<"Substep "<<fac<<endl;
            double Energy=0.0;
            get_q(potentials[i]*(1.f-fac)+potentials[i+1]*fac,lattice_dens1,q,eigvec,Energy);


            double fitness=0;

            double anyon0_charge=0;
            double anyon1_charge=0;


            for (uint64_t k = 0; k < N_sites; ++k)
            {
                if (anyon_sites[i][k]==1 || anyon_sites[i+1][k]==1)//But I don't care where the charge of the anyons are, as long as the charge in total adds up to 0.5
                    anyon0_charge+=std::abs(q(k));
                else if (anyon_sites[i][k]==2 || anyon_sites[i+1][k]==2)
                    anyon1_charge+=std::abs(q(k));
                else
                    fitness+=std::abs(q(k));//Whatever charge is here IS WRONG
            }

/*
            if (j%(Max/5)==0)
            {
                cout<<"ON SITE "<<i<<" : "<<j<<endl;
                for (uint y = 0; y < nh; ++y)
                {
                    for (uint x = 0; x < nw; ++x)
                    {
                        uint k = x+y*nw;
                if (anyon_sites[i][k]==1 || anyon_sites[i+1][k]==1)//But I don't care where the charge of the anyons are, as long as the charge in total adds up to 0.5
                    cout<<":"<<std::abs(q(k));
                else if (anyon_sites[i][k]==2 || anyon_sites[i+1][k]==2)
                    cout<<"|"<<std::abs(q(k));
                else
                    cout<<" "<<std::abs(q(k));//Whatever charge is here IS WRONG
                    }
                    cout<<endl;
                }

            cout<<endl;

            }*/


            double q_fitness =fitness;
            fitness+=std::abs(anyon0_charge-0.5)+std::abs(anyon1_charge-0.5);


            emergency_redo =false;

            if (i==0 && fac==0)
            {
                first_eigvec=eigvec;
                cerr<<(i+fac)<<" "<<fitness<<" "<<q_fitness<<" "<<std::abs(anyon0_charge-0.5)<<" "<<std::abs(anyon1_charge-0.5)<<" 1.0 "<<Energy<<endl;
            }
            else
            {
                complex<double> overlap= dot(p_eigvec.t(),eigvec);

                double phase = std::arg(overlap);
                complex<double> phase_offset = {cos(phase),-sin(phase)};
                //Fix phase

                for (complex<double>& D : eigvec)
                    D*=phase_offset;

                overlap= dot(p_eigvec.t(),eigvec);

                cout<<"Overlap "<<overlap<<endl;


                if (overlap.real()<0.995 && !fin_final )
                {
                    cout<<"EMERGENCY STEPSIZE DECRESE "<<endl;
                    emergency_redo =true;
                }
                else
                    cerr<<(i+fac)<<" "<<fitness<<" "<<q_fitness<<" "<<std::abs(anyon0_charge-0.5)<<" "<<std::abs(anyon1_charge-0.5)<<' '<<overlap.real()<<" "<<Energy<<endl;
            }

            if (!emergency_redo)
            {

                q_all_file<<fac<<endl;
                for (uint y = 0; y < nh; ++y)
                {
                    for (uint x = 0; x < nw; ++x)
                    {
                        q_all_file<<q(x+y*nw)<<' ';
                    }
                    q_all_file<<endl;
                }
                q_all_file<<endl;
                if (uint (fac*20)>uint(prev_fac*20) || fac==0 || fin_final )
                {
                    print_density_data(q_file,q , nw, nh);
                    q_file<<endl;

                    cout<<"Saved "<<i<<' '<<fac<<endl;
                    cout<<"With fitness "<<fitness<<endl;
                }



                p_eigvec=eigvec;




                prev_fac=fac;
                ++steps_since_emergency;

                if (zoom_level>0)
                    if (steps_since_emergency>= pow(2,zoom_level))//Hopefully the emergency is over
                    {
                        zoom_level =0;
                        stepsize=default_step_size;
                        cout<<"RETURNING TO DEFAULT ZOOM"<<endl;
                    }


            }
            else
            {
                if (zoom_level==20)
                {
                    cout<<"ZOOM 20 EXCEEDED, THIS IS NOT POSSIBLE !"<<endl;
                    return 1;
                }
                stepsize*=0.5;
                ++zoom_level;
                cout<<"GO TO ZOOM LEVEL "<<zoom_level<<" AND REDO STEP WITH NEW STEP-SIZE "<<stepsize<<endl;
                fac=prev_fac;
                steps_since_emergency=0;
            }



        }

        cout<<"DID "<<steps_here<<" DYNAMIC STEPS"<<endl;
    }


    cout<<"Final step "<<endl;

    complex<double> overlap= dot(p_eigvec.t(),first_eigvec);

    cout<<"Final overlap "<<overlap<<endl;
    q_file.close();
    q_all_file.close();

    return 0;
}
