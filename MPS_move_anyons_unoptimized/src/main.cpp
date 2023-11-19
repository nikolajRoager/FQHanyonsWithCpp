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

using namespace std; //Namespace containing in and out stream, std::vector (dynamically sized list) and complex numbers
using namespace arma;//Namespace containing matrices and vectors (actual vectors)



#include "approx.hpp"
#include "mps.hpp"
#include "mpo.hpp"

int main(int argc, char* argv[])
{
    //Read arguments
    if (argc<7 )
    {
        cout<<"Invalid, argument, need "<<argv[0]<<" width height N_particles D max_sweeps q_output_file anyon_location_file (Add more files, to get more states in the same output)"<<endl;
        return 1;
    }

    cout<<"Execute :";

    for (uint i = 0; i<argc; ++i)
        cout<<' '<<argv[i];
    cout<<endl;




    //If this fails, we get 0
    uint64_t nw = atoi(argv[1]);
    uint64_t nh = atoi(argv[2]);

    uint64_t N_sites = nh*nw;


    uint64_t N_particles = atoi(argv[3]);
    uint64_t N_anyons = 2;


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


    uint64_t D = atoi(argv[4]);

    uint64_t max_sweeps = atoi(argv[5]);

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


    vector<double> V(nw*nh);

    double phi = (N_particles+0.5*N_anyons)*2.0/((nw-1)*(nh-1));

    //cout<<"anyon sites "<<endl;
    //for (uint i = 0; i < nw*nh; ++i)
    //{
    //    V[i]=0;
    //    cout<<anyon_sites[i]<<' '<<endl;
    //}

    double t = 1.0;//Use t as basic unit

    cout<<"Creating MPO"<<endl;
    vector<complex<double> > tx(N_sites);
    vector<complex<double> > ty(N_sites);

    for (uint jy = 0; jy<nh; ++jy)
        for (uint jx = 0; jx<nw; ++jx)
        {
            uint j = jx*nh+jy;
            if (jx+1<nw)
            {
                tx[j]=t;
            }
            if (jy+1<nh)
            {
                double angle = M_PI*2.0*phi*double(jx);
                ty[j]= {t*cos(angle),t*sin(angle)};
            }
        }

    mpo hofstadter_mpo(nw,nh,V,tx,ty);




    cout<<"Creating MPS (no potentials)"<<endl;
    mps mps0 (N_sites,N_particles+1,D,true);
//    mps0.step_left_QR();
//    while (mps0.step_right_QR())
//    {cout<<"Step right"<<endl;}

    cout<<"Get MPO expectation value "<<endl;
    complex<double> expectation_value =hofstadter_mpo.get_expectation_value(mps0);
    cout<<"got expectation value "<<expectation_value <<endl;
    //mps0.DEBUG_PRINT_ALL();
    //hofstadter_mpo.DEBUG_PRINT_ALL();


    cout<<"Starting MPO ground search (no potential)"<<endl;
    double energy = 0.0;
    mps0 = hofstadter_mpo.find_ground(mps0,energy,max_sweeps);

    cout<<"Ground state energy "<<energy<<endl;


    mpo n_mpo(nw,nh);

    vector<double> n0(N_sites);

    ofstream q_out(argv[6]);




    for (uint y = 0 ; y < nh; ++y)
    {
        for (uint x = 0 ; x < nw; ++x)
        {
            uint i=x*nh+y;

            n0[i]=n_mpo.get_expectation_value(mps0,i).real();
        }
    }





    q_out<<N_sites<<' ';
    for (uint x = 0 ; x < nw; ++x)
        q_out<<x-0.5<<' '<<(x+0.5)<<' ';

    q_out<<endl;
    for (uint y = 0 ; y < nh; ++y)
    {
        //Top and bottom of the cell
        for (uint j = 0; j < 2; ++j)
        {

            q_out<<(y+j-0.5);
            for (uint x = 0 ; x < nw; ++x)
            {

                uint i=x*nh+(nh-y-1);//Plot upside down, for gnuplot to plot it the way I want


                q_out<<' '<<n0[i]<<' '<<n0[i];
            }
            q_out<<endl;
        }
    }





    for (uint arg = 7; arg<argc; ++arg)
    {
        vector<uint64_t> anyon_sites(N_sites);//The sites containing 0: no anyons, 1: first anyon 2: second anyon

        ifstream anyon_file(argv[arg]);

        if (!anyon_file.is_open())
        {
            cout<<"Could not open anyon file "<<argv[arg]<<endl;
            return 1;
        }

        for (uint y = 0 ; y < nh; ++y)
            for (uint x = 0 ; x < nw; ++x)
            {
                uint64_t I =0;
                uint i=x*nh+y;
                if (! (anyon_file>>I))
                {
                    cout<<argv[arg]<<" did not match number of sites (site "<<i<<" could not be read)" <<endl;
                    return 1;
                }

                if (I>2)
                    cout<<"Site "<<i<<" could not be read in "<<argv[arg]<<" is outside range, should be 0 (no anyons), 1 or 2 (either anyon), got "<<I<<endl;

                anyon_sites[i]= I;
            }


        anyon_file.close();





        //Now add in a potential
        {
            //Just make the Naïve guess that the potential at the location of the anyons is 1 and that it is 0 everywhere else.
            for (uint64_t i = 0; i < N_sites; ++i)
            {

                V[i]=(anyon_sites[i]==0) ? 0 : 1;

            }
        }

        cout<<"Creating MPS (with potentials)"<<endl;
        mps mps_anyons (N_sites,N_particles,D,true);

        hofstadter_mpo.set_potential(V);

        cout<<"Starting MPO ground search (with potential)"<<endl;
        energy = 0.0;
        mps_anyons = hofstadter_mpo.find_ground(mps_anyons,energy,max_sweeps);

        cout<<"This has ground state energy "<<energy<<endl;


        q_out<<endl;
        q_out<<N_sites<<' ';
        for (uint x = 0 ; x < nw; ++x)
            q_out<<x-0.5<<' '<<(x+0.5)<<' ';
        q_out<<endl;



        vector<double> n1(N_sites);


        double fitness=0;

        double anyon0_charge=0;
        double anyon1_charge=0;

        for (uint y = 0 ; y < nh; ++y)
        {
            for (uint x = 0 ; x < nw; ++x)
            {
                uint i=x*nh+(nh-y-1);//Plot upside down, for gnuplot to plot it the way I want
                n1[i]=n_mpo.get_expectation_value(mps_anyons ,i).real();

                double qi=n0[i]-n1[i];
                if (anyon_sites[i]==0)
                    fitness+=std::abs(qi);//Whatever charge is here IS WRONG
                else if (anyon_sites[i]==1)//But I don't care where the charge of the anyons are, as long as thje charge in total adds up to 0.5
                    anyon0_charge+=std::abs(qi);
                else if (anyon_sites[i]==2)
                    anyon1_charge+=std::abs(qi);

            }
        }




        fitness+=std::abs(anyon0_charge-0.5)+std::abs(anyon1_charge-0.5);




        for (uint y = 0 ; y < nh; ++y)
        {
            //Top and bottom of the cell
            for (uint j = 0; j < 2; ++j)
            {
                q_out<<(y+j-0.5);
                for (uint x = 0 ; x < nw; ++x)
                {
                    uint i=x*nh+(nh-y-1);//Plot upside down, for gnuplot to plot it the way I want

                    q_out<<' '<<(n0[i]-n1[i])<<' '<<(n0[i]-n1[i]);
                }
                q_out<<endl;
            }
        }

    }



    q_out.close();





    return 0;
}
