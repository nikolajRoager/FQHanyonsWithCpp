#include "minimize.hpp"

#include<cstdint>
#include<iostream>
#include<string>
#include<armadillo>
#include<functional>
#include<thread>

using state = uint64_t;
using namespace std; //Namespace containing in and out stream, std::vector (dynamically sized list) and complex numbers
using namespace arma;//Namespace containing matrices and vectors (actual vectors)




vec qnewton( function<double(vec)> F, vec x0, size_t& step,size_t max_steps, double acc, bool verbose)
{

    //I want ALL THE THREADS
    uint64_t max_threads = thread::hardware_concurrency();

    if (max_threads==0)
    {
        cout<<"WARNING, could not get number of available cores, running single threaded"<<endl;
        max_threads=1;
    }
    else
        max_threads =1;
    uint64_t calls=0;

    auto func = [&F,&calls](vec V) -> double
    {
        ++calls;
        return F(V);
    };


    double little = pow(2,-26);
    double alpha = 1e-4;//Used in acceptance condition for backtracking
    double epsilon = 1e-9;//Used to see if we should reset hessian

    vec x = x0;//Calls copy constructor by default

    //Read the number of parameters the function outputs, and the number of parameters in the input
    uint64_t n = x.size();


    bool redo=true;//I have an emergency redo function, in case deltax is too small, in some rare cases, I have seen deltax be treated as a 0, breaking the algorithm uppon division, if that happen, we will retry here

    uint64_t total_calls =1;//We will start with one call, before starting optimization

    vec best_x = x;
    double best_fx=-1;
    double fx =-1;
    double delta_fx =-1;
    uint64_t break_counter=0;//If the break condition is only fulfilled for one singular step, don't break
    uint64_t break_counter_max=16;//But if it keeps happening, do break

    while (redo)
    {
        redo=false;//In 99% of cases, this will only need to be done once, but maybe we are on a computer where machine epsilon is smaller than expected.

        step = -1;
        fx=func(x);
        if (best_fx==-1)
        {
            //Sometimes we end up going UP the hill, save the best we ever found
            best_fx=fx;
            best_x=x;
        }
        mat B = mat(n,n,fill::eye);//Inverse hessian, start guess: identity

        //We will save the gradient between steps, as we need the gradient both at the start and at the end
        vec GradFx = vec(n);


        if (verbose)
        {

            cout<<"Start f(x)="<<fx<<endl;
        }


        do
        {//Each step
            calls = 0;
            ++step;
            //PART 1: Calculate the Gradient
            // ∇ f(x)
            //
            // PART 2: backtracking:
            // Let Delta x = - B ∇ f(x)
            // Let s = lambda Delta x
            // Perform backtracking lambda -> lambda/2 and update s until f(x+s)<f(x)+alpha dot(s,∇ f(x)) (or until lambda too small)
            //
            // Symmetric Rank 1 update:
            //
            // CALCULATE ∇f(x + s) again !!?
            //
            // Let y = ∇f(x + s)-∇f(x) and u = s- By
            //
            // if dot(u,y)>epsilon
            //     let B = outer_product(u,u)/dot(u,y)
            // else
            //     let B = 1
            // endif


            //If this is the first step, we need to get gradient here matrix
            if  (step==0)
                for (uint64_t i = 0; i < n; ++i)
                {
                    //Get the input, offset on the relevant element
                    vec offsetx = x;
                    double deltax = std::max(std::abs(x[i])*little,0.25*little);//The tiny step to use when calculating the Jacobi matrix, note I do not enjoy 0 divisions, and it is possible that the user guesses that x[k]=0, so I add a floor to deltax

                    offsetx[i]+=deltax;

                    double foffsetxk=func(offsetx);
                    offsetx[i]+=deltax;
                    GradFx[i] = (foffsetxk-fx)/deltax;

                    if(std::isinf(GradFx[i]))//This is very very rare
                    {
                        //CRAP, deltax was too small, and we got a 0 division, we must retry with a larger little number... unless we have done that too many times already

                        if (little>0.25)
                            cout<<("ABORT: Got 0 division in Quasi Newton methods gradient too many times")<<endl;

                        cout<<("WARNING: Got 0 division in Quasi Newton method gradient, this may be caused by your computer having less precision for double point numbers than expected. Will retry with a larger Delta x")<<endl;
                        redo = true;
                        x = x0;//Reset our guess, and try again with a larger little number
                        little*=2;


                        break;
                    }
                }



            try
            {

            if (!redo)//If finding the gradient worked
            {
                if (verbose)
                {
                    //cout<<"Step "<<step<<endl;
                    //cout<<"x= "<<x<<endl<<endl;
                    //cout<<"numeric = ∇ f(x)="<<GradFx<<endl<<endl<<endl;
                }

                //Now need to do backtracking

                // Let Delta x = - B ∇ f(x)
                vec Dx = -B*GradFx;


                // Let s = lambda Delta x
                double lambda = 1;
                vec S = Dx;
                vec x_new = x+S;


                bool linesearch_fail = true;

                double fx_new = fx;


                for (uint64_t i = 0; i < 32; ++i)
                {

                    // Perform backtracking lambda -> lambda/2 and update s until f(x+s)<f(x)+alpha dot(s,∇ f(x)) (or until lambda too small)
                    fx_new=func(x_new);//Hopefully we won't calls this 32 times
                    if (fx_new < fx+ alpha*dot(S,GradFx)  )
                    {
                        linesearch_fail =false;
                        break;
                    }

                    lambda*=0.5;
                    S = lambda*Dx;
                    x_new = x+S;
                }

                x = x_new;
                delta_fx=fx_new-fx;
                fx=fx_new;

                if (fx<best_fx || best_fx==-1)
                {
                    //Sometimes we end up going UP the hill, save the best we ever found
                    best_fx=fx;
                    best_x=x;
                }

                //With enough steps, this will actually work without resetting B. But we can do better now entering SR1 update


                // CALCULATE ∇f(x + s) for the next step

                vec GradFx_new = vec(n);

                uint64_t i=0;//Current index


                auto get_gradi =[&func,&fx](double* out, vec offsetx, double deltax) -> void

                {
                    double foffsetxk=func(offsetx);
                    *out = (foffsetxk-fx)/deltax;
                };




                while (i<n)
                {




                    std::thread myThreads[max_threads];
                    uint64_t active_threads =0;
                    for (uint64_t j=0; (j<max_threads ); j++){


                        //Get the input, offset on the relevant element
                        vec offsetx = x;
                        double deltax = std::max(std::abs(x[i])*little,0.25*little);//The tiny step to use when calculating the Jacobi matrix, note I do not enjoy 0 divisions, and it is possible that the user guesses that x[k]=0, so I add a floor to deltax

                        offsetx[i]+=deltax;

 //                       AllThreads.push_back(thread(f1,1/*, func,fx,GradFx_new[i],offsetx,deltax*/));

                        //double foffsetxk=func(offsetx);
                        //GradFx_new[i]= (foffsetxk-fx)/deltax;


                        myThreads[j]=thread(get_gradi,&(GradFx_new[i]),offsetx,deltax);
                        ++active_threads ;
                        ++i;

                        if (i>=n)
                            break;
                    }
                    for (uint64_t j=0; j<active_threads ; j++){
                        myThreads[j].join();
                    }





                }

                //Do sanity check
                for (uint64_t i = 0;   i<n; ++i)
                {

                    if(std::isinf(GradFx_new[i]))//This is very very rare
                    {
                        //CRAP, deltax was too small, and we got a 0 division, we must retry with a larger little number... unless we have done that too many times already

                        if (little>0.25)
                            cout<<"ABORT: Got 0 division in Quasi Newton methods gradient too many times"<<endl;

                        cout<<"WARNING: Got 0 division in Quasi Newton method gradient, this may be caused by your computer having less precision for double point numbers than expected. Will retry with a larger Delta x"<<endl;
                        redo = true;
                        x = x0;//Reset our guess, and try again with a larger little number
                        little*=2;


                        break;
                    }
                }

                if (linesearch_fail)
                {
                    B = mat(n,n,fill::eye);
                    if (verbose)
                    {
                        cout<<"RESET B AS LINESEARCH FAILED"<<endl;
                        //cout<<"B="<<B<<endl;
                    }

                }
                else
                {

                    // Let y = ∇f(x + s)-∇f(x) and u = s- By
                    vec y = GradFx_new-GradFx;
                    vec u = S-B*y;


                    if (verbose)
                    {
                        //cout<<endl;
                        //cout<<"∇f(x + s)="<< GradFx_new<<endl;
                        //cout<<"y="<<y<<endl;
                        //cout<<endl;
                        //cout<<"u="<<u<<endl;
                        //cout<<"So dot(y,u)={u.dot(y)}"<<endl;
                    }

                    if (std::abs(dot(u,y))>epsilon)
                    {

                        B = B+ mat(u)*trans(mat(u))/dot(u,y);
                        if (verbose)
                        {
                            //cout<<"UPDATED B"<<endl;
                            //cout<<B<<endl;
                        }
                    }
                    else
                    {
                        if (verbose)
                        {
                            //cout<<"KEPT B"<<endl;
                            //cout<<B<<endl;
                        }
                    }

                }

                //Now, we already have whatever the gradient should be next step.
                GradFx=GradFx_new;

            }



            }
            catch(...)//Crap, lets just hope this does not happen again
            {

                if (verbose)
                {
                    cout<<"PREMATURE BREAK (matrix could not be diagonalized), the optimization can not continue."<<endl;
                }
                redo=false;
                break;
            }

            if (verbose)
            {
                cout<<"Step "<<step<<"\t f(x)="<<fx<<" Function calls: "<<calls<<endl;

            }

            if (step>=max_steps)//I am not sure how you would have > max steps ... ok actually you can, if we redo at the second to last step
            {
                if (verbose)
                    cout<<"ALGORITHM FAILED: too many steps"<<endl;
                break;
            }


            if (norm(GradFx)<acc)
            {
                if (verbose)
                    cout<<"BREAK DUE TO: gradient reached target"<<endl;
                break;
            }
            if (std::abs(delta_fx)<acc)
            {
                ++break_counter;


                if (break_counter>break_counter_max)
                {
                    if (verbose)
                        cout<<"BREAK DUE TO: Change too small"<<endl;
                    break;
                }

                if (verbose)
                    cout<<"Change too small, if this continues we will break ("<<break_counter<<'/'<<break_counter_max<<")"<<endl;

            }
            else
                break_counter=0;


            total_calls+=calls;

        }//check if we should stop, from too many steps, or reaching the goal
        while (true);

    }

    if (verbose)
    {
        cout<<"Finished in "<<total_calls<<" calls"<<endl;
        if (best_fx!=fx)
            cout<<"The final step ("<<fx<<") is worse than the best ever step ("<<best_fx<<"), returning to the best ever step"<<endl;
    }


    return best_x;
}



//Mirror function to maximize instead
vec max_qnewton( function<double(vec)> func , vec x0, size_t& steps, double acc, bool verbose)
{
    return  qnewton( [func](vec x) -> double {return -func(x);},x0,steps,acc,verbose);
}

//Versions without counting number of steps
vec qnewton( function<double(vec)> func , vec x0, double acc, bool verbose)
{
    size_t steps;
    return  qnewton( func,x0,steps,acc,verbose);
}

//Mirror function to maximize instead
vec max_qnewton( function<double(vec)> func , vec x0, double acc, bool verbose)
{
    size_t steps;
    return  qnewton( [func](vec x) -> double {return -func(x);},x0,steps,acc,verbose);
}


vec downhill_simplex(function<double(vec)> F, vec x0, size_t& step, double acc, bool verbose)
{
    step =0;
    uint64_t calls=0;
    auto func = [&F,&calls](vec V) -> double
    {
        ++calls;
        return F(V);
    };


    //Read the number of parameters the function outputs, and the number of parameters in the input
    uint64_t n = x0.size()+1;
    uint64_t max_steps=1024;//should never need that much


    vector<vec> vertices(n+1);//The simplex vertices + the centroid
    vector<double> fs = vector<double>(n);

    for (uint64_t i = 0; i<n ; ++i)
    {
        vertices[i] = x0;//Calls copy constructor by default


        if (i!=0)
        {
            vertices[i][i-1]+=1.0;//Start with an offset in each dimension, except for the first point
        }
        fs[i] = func(vertices[i]);

    }


   uint64_t lowest  = 0;//Which of the points are lowest and highest
   uint64_t highest = 0;

    //Check if any of the other actual points are higher/lower, and get the first centroid
    vertices[n] = vertices[0]/(n-1);//The centroid is (vertices[0]+...vertices[n-1])/(n-1) over all points which are NOT the highest

    for (uint64_t i = 1; i < n; ++i)
    {
        if (fs[i] < fs[lowest])
        {
            lowest = i;
        }
        if (fs[i] >= fs[highest])//The = here guarantees that the highest and lowest will not be the same, even if all points are identical
        {
            highest = i;
        }
        vertices[n] += vertices[i]/(n-1);//The centroid is (vertices[0]+...vertices[n-1])/(n-1) over all points which are NOT the highest
    }

    //Remove the highest point from the centroid sum
    vertices[n] -= vertices[highest]/(n-1);

    uint64_t total_calls=calls;

    bool do_end = false;
    bool reduced = false;//Was the last step reduction
    do
    {//Each step
        ++step;

        calls=0;

        //Try reflection
        vec new_vertex = vertices[n]-(vertices[highest]-vertices[n]);//This takes highest -> centroid - (vector from centroid to highest)
        double new_vertex_f = func(new_vertex);

        //Did it work
        if (new_vertex_f<fs[highest])
        {

            //If reflection works perfectly, try expansion
            if (new_vertex_f<fs[lowest])
            {//Oh we are going the right way, try expansion
                lowest = highest;//Regardless of what happens, this is the lowest point now

                //Try reflection
                vec expanded= vertices[n]-2*(vertices[highest]-vertices[n]);//This takes highest -> centroid - 2*(vector from centroid to highest)
                double expanded_f = func(expanded);

                if (expanded_f<new_vertex_f)
                {//accept expansion
                    new_vertex = expanded;
                    new_vertex_f = expanded_f;
                    if (verbose)
                        cout<<"Accept Expansion"<<endl;

                }
                else
                    if (verbose)
                            cout<<"Reject Expansion, Accept Reflection"<<endl;


            }
        }
        else//REJECT reflection
        {//Try contraction

            //Try contraction
            new_vertex = vertices[n]+0.5*(vertices[highest]-vertices[n]);//This takes highest -> centroid + 1/2 (vector from centroid to highest)
            new_vertex_f = func(new_vertex);

            if (new_vertex_f <fs[highest])
            {//accept expansion
                if (verbose)
                    cout<< "Accept Contraction"<<endl;

                if (new_vertex_f <fs[lowest])
                {
                    lowest = highest;//If, by change this new point is actually the lowest, reset the lowest marker to that
                }

            }
            else
            {
                if (verbose)
                    cout<< "Do reduction"<<endl;

                //Shrink everything towards the lowest point
                for (uint64_t i = 0; i < n; ++i)
                {
                    if (i != lowest)//The = here guarantees that the highest and lowest will not be the same, even if all points are identical
                    {
                        //This becomes the midpoint between this and the lowest
                        vertices[i] = 0.5*(vertices[i]+vertices[lowest]);



                    }
                }

                //also reset the centroid and the highest and lowest point, using a full recalculating as EVERYTHING changed
                lowest  = 0;
                highest = 0;

                //Check if any of the other actual points are higher/lower, and get the first centroid
                vertices[n] = vertices[0]/(n-1);//The centroid is (vertices[0]+...vertices[n-1])/(n-1) over all points which are NOT the highest

                if (norm(vertices[0]-vertices[n-1])<acc)
                {
                    if (verbose)
                        cout<<"Distance to neighbour within limit"<<endl;
                    do_end=true;
                }

                for (uint64_t i = 1; i < n; ++i)
                {
                    if (fs[i] < fs[lowest])
                    {
                        lowest = i;
                    }
                    if (fs[i] >= fs[highest])//The = here guarantees that the highest and lowest will not be the same, even if all points are identical
                    {
                        highest = i;
                    }
                    vertices[n] += vertices[i]/(n-1);//The centroid is (vertices[0]+...vertices[n-1])/(n-1) over all points which are NOT the highest

                    if (norm(vertices[i]-vertices[i-1])<acc)
                    {
                        if (verbose)
                            cout<<"Distance to neighbour within limit"<<endl;
                        do_end=true;
                    }

                }
                //Remove the highest point from the centroid sum
                vertices[n] -= vertices[highest]/(n-1);


                reduced = true;



            }







        }

        if (!do_end && !reduced)//If we rand reduction, we already reset the lowest and highest, if we should end, we did not find a new vertex
        {//Reset the centroid and highest and lowest each step

            //If the new point is not the the lowest, the lowest is guaranteed to be unchanged

            vertices[highest] = new_vertex;
            //Reset highest and centroid
            fs[highest] =  new_vertex_f;

            //The highest may have changed, our first guess is that the new_vertex poitn is still the highest
            uint64_t new_vertex_point = highest;
            for (uint64_t i = 0; i < n; ++i)
            {
                if (fs[i] >= fs[highest])//The = here guarantees that the highest and lowest will not be the same, even if all points are identical
                {
                    highest = i;
                }
            }

            if (new_vertex_point != highest)//Only change the centroid if the new point we added is not still the highest
            {
                //Remove the point which now is the highest from the centroid sum
                vertices[n] -= vertices[highest]/(n-1);
                vertices[n] += vertices[new_vertex_point]/(n-1);
            }

            //Check if the point we just moved is too close to one neighbor (it is the only distance which changed, and it always has the same distance to both neighbour )

            int neighbour =(new_vertex_point == 0) ? n-1 : new_vertex_point -1;


            if (norm(vertices[new_vertex_point]-vertices[neighbour ])<acc)
            {
                if (verbose)
                    cout<<"Distance to neighbour within limit"<<endl;
                do_end=true;
            }
        }







        if (step>=max_steps)//I am not sure how you would have > max steps
        {
            if (verbose)
                cout<<"ALGORITHM FAILED: too many steps"<<endl;
            do_end=true;
        }

        if (verbose)
        {
            cout<<"Step "<<step<<" ; calls : "<<calls<<" ; highest point : "<<fs[highest]<<" lowest point : "<<fs[lowest]<<endl;
        }
        total_calls+=calls;
    }
    while (!do_end);


    if (verbose)
    {
        cout<<"Finished in "<<total_calls<<" calls"<<endl;
    }

    return vertices[lowest];
}



//Mirror function to maximize instead
vec max_downhill_simplex( function<double(vec)> func , vec x0, size_t& steps, double acc, bool verbose)
{
    return  downhill_simplex( [func](vec x) -> double {return -func(x);},x0,steps,acc,verbose);
}

//Versions without counting number of steps
vec downhill_simplex( function<double(vec)> func , vec x0, double acc, bool verbose)
{
    size_t steps;
    return  downhill_simplex( func,x0,steps,acc,verbose);
}

//Mirror function to maximize instead
vec max_downhill_simplex( function<double(vec)> func , vec x0, double acc, bool verbose)
{
    size_t steps;
    return  downhill_simplex( [func](vec x) -> double {return -func(x);},x0,steps,acc,verbose);
}



/*
//Taking same input as before, use alglib because I care more about having a well tested and optimized optimization algorithm, than makign everything from scratch. the limited memory Broyden-Fletcher-Goldfarb-Shanno is a quasi newton method which should hopefully make many fewer calls.

#include <libalglib/stdafx.h>
#include <libalglib/optimization.h>

using namespace alglib;


//EVIL HACK, we NEED to convert a lambda function, which has CAPTURED all the setup parameters, into a simple void pointer, we will just define the things the function NEED to know in global space then
function<double(vec)> F_evil_hack;
uint64_t n =0;
uint64_t al_calls =0;



void libalg_func(const real_1d_array &x, double &func, void *ptr)
{
    ++al_calls;
    //We have to do a data copy, which absolutely sucks!... ok this is just the allocation of a single vector, it is not going to break anything... but still
    const vec V(x.getcontent(),n);
    func = F_evil_hack(V);

    cout<<"f="<<func<<" : calls="<<al_calls<<endl;
}


vec libalg_lbfgs(function<double(vec)> F, vec x0, uint64_t steps, double acc,bool verbose)//I will keep this third party library contained wihtin minimize.cpp alone, and convert the output back to armadillo style vectors, just to avoid confusion, armadillo is generally superior when it comes to matrix operations, hence why it is still my main library
{
    al_calls =0;
    steps=0;
    F_evil_hack = F;//copy the function over, so that I libalg doesn't get its fealings hurt.

    al_calls=0;
    n = x0.size();

    double epsg = acc;//Gradient break condition
    double epsf = 0;//function difference condition 0= nope
    double epsx = 0;//stepsize break condition 0=nope
    double diffstep = 1.0e-6;
    minlbfgsstate state;
    minlbfgsreport rep;

    real_1d_array x;
    x.setcontent(n,x0.memptr());//Does not borrow pointer, but creates an independent object

    minlbfgscreatef(1, x, diffstep, state);
    minlbfgssetcond(state, epsg, epsf, epsx, 256);

    if (verbose)
        cout<<"Calling alglib optimization, this might take quite a while ... maybe even hours"<<endl;


    alglib::minlbfgsoptimize(state, libalg_func);
    minlbfgsresults(state, x, rep);
    steps = (uint64_t)rep.iterationscount;
    //Again, copies data


    vec arma_out(x.getcontent(),n);

    if (verbose)

    {
        cout<<"alglib finished optimization in "<<steps<<" steps , with "<<al_calls <<" function calls with output status: "<<endl;
        switch(int(rep.terminationtype))
        {
            case 1:
                cout<<"-(SUCCESS) Function change reached target "<<endl;
            break;
            case 2:
                cout<<"-(INTERRUPT) stepsize too small"<<endl;
            break;
            case 4:
                cout<<"-(SUCCESS) Gradient norm reached target "<<endl;
            break;
            case 5:
            default:
                cout<<"-(INTERRUPT) Max steps taken"<<endl;
            break;
            case 7:
                cout<<"-(INTERRUPT) The library decided it has reached a minimum"<<endl;
            break;
            case 8:
                cout<<"-(INTERRUPT) Received manual interrupt command."<<endl;
            break;
            case -8:
                cout<<"-(FAIL) Zero division or infinity encountered in function"<<endl;
            break;
            case -2:
                cout<<"-(INTERRUPT) Serious rounding errors encountered"<<endl;
            break;
            case -1:
                cout<<"-(FAIL) Parameters of function did not match input"<<endl;
            break;
        }
    }

    return arma_out;
}*/
