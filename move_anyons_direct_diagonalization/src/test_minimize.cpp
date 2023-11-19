#include<cstdint>
#include<iostream>
#include<string>
#include<armadillo>

#include"minimize.hpp"

using state = uint64_t;
using namespace std; //Namespace containing in and out stream, std::vector (dynamically sized list) and complex numbers
using namespace arma;//Namespace containing matrices and vectors (actual vectors)


//double precision approximation
bool approx(double a,double b,double tau=1e-9,double eps=1e-9)
{
    if (std::abs(a-b)<2*tau)
        return true;
    if (std::abs(a-b)/(std::abs(a)+std::abs(b))<2*eps)
        return true;
    return false;
}


double simpleGauss(vec X)
{
    return exp( -(X[0]-2)*(X[0]-2));
}

int main(int argc, char* argv[])
{

    cout<< "----------------------------------------------"<<endl;
    cout<< "Demonstrating Quasi newton minimization method"<<endl;
    cout<< "----------------------------------------------"<<endl;


    bool verbose = false;


    for (int i = 0; i < argc; ++i)
        if (string(argv[i]).compare("-v")==0)
            verbose = true;

    cout <<"Demonstration one, find max of f(x)=exp(-(x-2)^2) (true solution x=2)"<<endl;

    cout<<"Starting at x0 = (0.0) with precision 10^-5. Now Running ..."<<endl;

    size_t steps=0;


    vec out = max_qnewton(&simpleGauss,vec("0"),steps,1e-5,verbose );

    cout<<" Got root "<<out[0]<<" in "<<steps<<"steps "<<endl;

    if (approx(out[0],2,1e-5,1e-5))//vector double is understood  as a vector with this element only
    {
        cout<<"PASS this is within 10^-5 of 2.0"<<endl;
    }
    else
    {
        cout<<"FAIL this is not within 10^-5 of 2.0"<<endl;
        return 0;
    }


    cout<<"----------------------------------------------"<<endl;

    cout<<"Demonstration Rosenbrock, find minimum of f(x,y)=((1-x)^2+100*(y-x^2)^2) (minimum is at 1,1)"<<endl;
    cout<<"Starting at x0 = (-1,2) with precision 10^-5. Now Running ..."<<endl;

    out = qnewton([](vec X) -> double {return ((1-X[0])*(1-X[0])+100*(X[1]-X[0]*X[0])*(X[1]-X[0]*X[0]) );},vec("-1 2"),steps,1e-5,verbose);

    cout<<endl;
    cout<<"Got predicted root x="<<out<<endl;

    cout<<endl;
    if (approx(out[0],1,1e-5,1e-5) && approx(out[1],1,1e-5,1e-5))
    {
       cout<<"PASS this is within 10^-5 of (1,1)"<<endl;
    }
    else
    {
        cout<<"FAIL this is not within 10^-5 of (1,1)"<<endl;
        return 0;
    }

    cout<<"----------------------------------------------"<<endl;

    cout<<"Demonstration Himmelblau, find minimum of f(x,y)=(x^2+y-11)^2+(x+y^2-7)^2\nminima at: (3.0,2.0),(-2.805118, 3.131312),(-3.779310, -3.283186) and (3.584428, -1.848126)"<<endl;
    cout<<"Starting at x0 = (0,0) with precision 10^-5. Now Running ..."<<endl;

    out = qnewton([](vec X) -> double {return ((X[0]*X[0]+X[1]-11)*(X[0]*X[0]+X[1]-11)+(X[0]+X[1]*X[1]-7)*(X[0]+X[1]*X[1]-7));},vec("0 0"),steps,1e-5,verbose);

    cout<<endl;
    cout<<"Got predicted root x="<<out<<endl;

    cout<<endl;
    return 0;
}
/*


    WriteLine($"----------------------------------------------");
    WriteLine($);
"minimize.hpp"       Func<vector,double> f2 = (X)  => ((X[0]*X[0]+X[1]-11)*(X[0]*X[0]+X[1]-11)+(X[0]+X[1]*X[1]-7)*(X[0]+X[1]*X[1]-7));
    WriteLine($"Starting at x0 = (0,0) with precision 10^-5. Now Running ...");
    (vector root2, int steps2) = qnewton(f2,new vector(0,0),1e-5,verbose , new System.IO.StreamWriter("himmelblau.tsv"));

    WriteLine("");
    WriteLine(root2.getString("Got predicted root x="));

    WriteLine("");
    WriteLine($"In {steps2} steps: Has f(x)={f2(root2)}");
    WriteLine("");
    if ( approx(f2(root2),0.0,1e-5,1e-5))
    {
        WriteLine("PASS the function is within 10^-5 of 0, here; which is the known minimum");
    }
    else
    {
        PASS=false;
        WriteLine("FAIL the function is not within 10^-5 of 0, here; which is the known minimum");
    }
    WriteLine($"----------------------------------------------");
    if (PASS)
        WriteLine("ALL TESTS PASSED");
    else
        WriteLine("SOME TESTS FAILED");
*/

/*


        return 0;
    }
*/
