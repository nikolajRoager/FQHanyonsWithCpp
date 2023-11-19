#include<generate_states.hpp>
#include "approx.hpp"

//Has already been done in my header, but I like to include everything explicitly in every file (the compiler will ignore it)
#include<cstdint>
#include<vector>
#include<iostream>
#include<armadillo>

using namespace std;
using namespace arma;



uint64_t binSearch(const vector<state_t>& states, state_t target, uint64_t low ,uint64_t high)
{
    high = std::min(high,states.size()-1);//Truncate (default guess is always end of list)

    while (low!=high)
    {
        //Take average
        uint64_t midpoint =(low + high)>>1;//Bitshift right once = divide unsigned integer by 2 and round down
            if (target == states[midpoint])
            {
               return midpoint;
            }
            else if (target > states[midpoint]) // x is on the right side
                low = midpoint + 1;//THIS IS SAFE as low +high<= 2 high-1, so midpoint <high<end of array
            else//x is on the left side
                high = midpoint - 1;//THIS IS ALSO SAFE as low+high>=low>=0, but if midpoint=0 then x can never be on the left side, so there is no risk of overflow

    }

    return low;
}

//We need the largest integer we can get when working with factorials
uint64_t factorial(uint64_t n)
{
    uint64_t out = 1;
    for (;n>0;--n)
        out*=n;
    return out;
}

//There are likely smarter ways of doing this, but this method does at least not try to store all of n! which is likely to overflow. We do still need to calculate and store whichever is smaller of m! or (n-m)! which is up to (n/2)! overflow is still likely so I will just use uint64_t
uint64_t bin_coef(uint64_t n, uint64_t m)
{
    if (n>=m && m>=0)
    {
        //Overflow is extremely likely if we use our formula
        //factorial(n)/(factorial(m)*factorial(n-m));
        //But in practice we are taking (n)*(n-1)*(n-2)*... (m+1) or (n-m+1), whichever is larger, divided by (n-m)! or m! whichever is smaller, lets label the larger m

        uint64_t M;//Whichever is larger of m and n-m; (this is really just bin_coef(n,m)=bin_coef(n,n-m) )

        if (2*m > n)
        {
            M=m;
            m=n-m;//Use m to store whichever is smaller
        }
        else
            M=n-m;

        uint64_t out = 1;
        for (uint64_t i = n; i>M; --i)
            out*=i;//Calculate (n)*(n-1)*(n-2)*... (m+1) or (n-m+1), whichever is larger
        out/=factorial(m);//Divide with whatever is smaller

        return out;
    }
    else
        return 0;
}

//Print this state, just as a sanity check
void print_state(state_t S,  int w, int h)//Little Endian
{
    int N_sites = w*h;
    cout << "|"<<S<<">=|";
    for (int i = 0 ; i < N_sites; ++i)
        cout <<  ((S >> i) & 1);

    cout<<">"<<endl;

    //Print as a table

    int i=0;//Running index, start at 0
    for (int y = 0 ; y < h; ++y)
    {
        cout << " |";
        for (int x = 0 ; x < w; ++x)
        {
            if((S >> i) & 1)
                cout << " x";
            else
                cout << " o";
            ++i;
        }
        cout << " |"<<endl;
    }

}

void print_state_data(state_t S,  int w, int h)
{
    int N_sites = w*h;

    //Print as a table


    //Print this in the format, which makes gnuplot plot this as a surface plot
    cerr << N_sites<<' ';
    for (int x = 0 ; x < w; ++x)
        cerr << x <<' '<<(x+1)<<' ';

    cerr <<endl;
    int i=0;//Running index, start at 0
    for (int y = 0 ; y < h; ++y)
    {
        //Top and bottom of the cell
        for (int j = 0; j < 2; ++j)
        {
            cerr <<(y+j);
            for (int x = 0 ; x < w; ++x)
            {
                if((S >> (i+x)) & 1)
                    cerr <<" 1.0 1.0";
                else
                    cerr <<" 0.0 0.0";
            }

            cerr<<endl;
        }

        i+=w;
    }

}

        //From smallest to larges the states (example 3 partilces 6 states)
        //000111 0
        //001011 1
        //001101 1
        //001110 1
        //010011 2
        //010101 2
        //011001 2
        //011010 2
        //010110 2
        //011010 2
        //011100 2
        //...

        //The Way I generate this, in ascending order, is I start by moving the front particle left:
        // 0001xx
        // 001xxx
        // 01xxxx
        // 1xxxxx
        //Where xx... are all the ways of aranging M-1 particles on however many sites are left



void generate_recursive(vector<state_t>& states,//The list of states we are writing to
              state_t Template,//What shall we set the bits outside range to?
              uint64_t& N,           //The current index we have gotten to
              uint64_t N_sites,      //How many sites do we have
              uint64_t N_particles)  //How many particles should we place
{
    if (N_particles>1)//More particles, require recursion
    {
        //Loop through:
        // 0001xx
        // 001xxx
        // 01xxxx
        // 1xxxxx

        // Now we want to fill in x x x x x as the ways we can arrange Np-1 particles on i-1 sites, we will do this recursively
        for (uint64_t i = N_particles-1; i<N_sites; ++i)
        {
            //Now, we want to keep the 0001 or whatever in front of the xx, and if this is not the first level of recursion, keep anything before this also
            generate_recursive(states,Template | ((state_t)1)<<i,N,i,N_particles-1);
        }

    }
    else if (N_particles==1)
    {
        //Write down all single-particle states

        for (uint64_t i = 0; i < N_sites;++i)
            states[N+i]=Template | (((state_t)1)<<i);

        (N)+=N_sites;//Advance
    }
    else//The only way we get here, is if 0 particles was inserted to begin with, but that can still happen soo...
    {//Huh, no particles, well there is only one option then
        states[N]=Template;
        ++(N);
    }
}

//Default version
void generate(vector<state_t>& states,
              uint64_t& N_states,
              uint64_t N_sites,
              uint64_t N_particles)
{
    N_states = bin_coef(N_sites,N_particles);
    states=vector<state_t>(N_states,0);
    //Set everything to 0
    uint64_t N=0;
    generate_recursive(states,0,N,N_sites,N_particles);
}

void print_density_data(ofstream& Out,vec lattice_density,  uint w, uint h)//Assumes little Endian states
{
    uint N_sites = w*h;


    //Print this in the format, which makes gnuplot plot this as a surface plot
    Out<<N_sites<<' ';
    for (uint x = 0 ; x < w; ++x)
        Out<<x-0.5<<' '<<(x+0.5)<<' ';

    Out<<endl;
    for (uint y = 0 ; y < h; ++y)
    {
        //Top and bottom of the cell
        for (uint j = 0; j < 2; ++j)
        {

            Out<<(y+j-0.5);
            for (uint x = 0 ; x < w; ++x)
            {

                uint i=x+(h-y-1)*w;//Plot in reverse order, for gnuplot to plot it the way I want
                Out<<' '<<lattice_density(i)<<' '<<lattice_density(i);

            }
            Out<<endl;
        }

    }
}



//Print this state in a format gnuplot can plot
void print_superpositions_data(const vector<state_t>& states,vector<double> superpositions , uint64_t& n_states,  int w, int h)//Assumes little Endian states
{
    int N_sites = w*h;

    //Print as a table


    //Print this in the format, which makes gnuplot plot this as a surface plot
    cerr<<N_sites<<' ';
    for (int x = 0 ; x < w; ++x)
        cerr<<x-0.5<<' '<<(x+0.5)<<' ';

    cerr<<endl;
    for (int y = 0 ; y < h; ++y)
    {
        //Top and bottom of the cell
        for (int j = 0; j < 2; ++j)
        {

            cerr<<(y+j-0.5);
            for (int x = 0 ; x < w; ++x)
            {
                int i=x+(h-y-1)*w;//Plot in reverse order, for gnuplot to plot it the way I want
                double sum = 0;
                //Sum up all possible states
                for (int k = 0; k < n_states; ++k)
                {
                    if( (states[k] >> (i)) & ((state_t)1))
                        sum+=superpositions[k];
                }
                cerr<<' '<<sum<<' '<<sum;;

            }

            cerr<<endl;
        }

    }
}


//Get the density at all sites, I prefer to get this data as a matrix
vec get_density(const vector<state_t>& states,vector<double> superpositions , const uint64_t& n_states,  int w, int h)//Assumes little Endian states
{
    vec Out(w*h);//Starts as 0

    //This time I try the approach of using a running index, arguably *moderately* better
    uint64_t i=0;
    for (uint64_t y = 0 ; y < h; ++y)
    {
        for (uint64_t x = 0 ; x < w; ++x)
        {
            //Sum up all possible states
            for (uint64_t k = 0; k < n_states; ++k)
            {
                if( (states[k] >> (i+x)) & ((state_t)1))
                    Out(i+x)+=superpositions[k];
            }
        }
        i+=w;
    }

    return Out;
}
