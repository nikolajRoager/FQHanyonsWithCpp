/*
LICENCE GOES HERE
*/

#define DEBUG_BUILD

#include "IO.hpp"
#include "get_assets.hpp"

#include<string>
#include<iostream>
#include<fstream>
#include<optional>
#include<filesystem>
#include<exception>

#include<cstdint>


//I use these types so much that these aliases are well worth it
using ushort = uint16_t;
using uint = uint32_t;
using ulong = uint64_t;//I need to be sure that this is large enough because I use it for keeping time, and 16 bit will break in 1 minute and 5 seconds, this should last long enough for even the most dedicated of people

using namespace std;

namespace fs = std::filesystem;

//The main loop, I try to keep anything graphics related out of this.
int main(int argc, char* argv[])
{
    fs::path Assets = get_assets();
    //fs::path Fonts = Assets/"fonts";
    fs::path Keymaps  = Assets/"keymap.txt";
    fs::path Shaders  = Assets/"shaders";
    fs::path Textures = Assets/"textures";
    fs::path Audio    = Assets/"audio";
    fs::path Materials= Assets/"materials.json";



    cout << "Running with "<<argc<<" arguments, command executed: " << endl;

    for (int i = 0; i < argc; ++i)
    {
        cout << argv[i] << " " << flush;
    }

    if (argc!=6)
    {
        cout<<"Need arguments "<<argv[0]<<" w h min max file" << endl;
        return 1;
    }

    cout << endl;
    ushort cursor =0;

    //If this fails, we get 0
    uint64_t nw        = atoi(argv[1]);
    uint64_t nh        = atoi(argv[2]);
    float min_q = atof(argv[3]);
    float max_q = atof(argv[4]);

    uint64_t steps=0;
    cout<<"Starting main display"<<endl;
    try
    {
        IO::init("Animated charge",nw,nh,Textures,Audio,Shaders,Keymaps);
    }
    catch(std::exception&  error)
    {
        IO::end();
        cout<<"couldn't set up graphics engine: "<<error.what()<<endl;
        return 0;
    }





    ulong millis=0;
    ulong pmillis=IO::input_devices::get_millis();

    int mx=0;
    int my=0;


    cout<<"Start loop"<<endl;


    double fixed_dt = 1.0/60.0;

    ulong fixed_dt_millis = fixed_dt*1000;

    vector <vector<float> > q;
    vector <float> f;
    vector <uint> k;

    ifstream q_file(argv[5]);

    if (!q_file.is_open())
    {
        cout<<"Could not open "<<argv[5]<<endl;
    }

    bool do_break=false;
    uint K=0;
    double pf=0;
    while (true)
    {

        vector<float> this_q(nw*nh);

        double F=0;
        if(!(q_file>>F))
            break;



        for (uint j = 0; j < nw*nh; ++j)
            if (!(q_file>>this_q[j]))
            {
                do_break=true;
                break;
            }
        if (do_break)
            break;

        q.push_back(this_q);
        f.push_back(F);
        if (F<pf)
        {
            ++K;
        }
        else
        {
        }

        k.push_back(K);
        pf=F;
        ++steps;
    }
    cout<<"Playing "<<steps<<" Frames "<<endl;

    q_file.close();
    uint frame=0;
    uint millis_last_move_frame=0;

    uint ms_per_frame = 1/24.f;

    do
    {
        millis = IO::input_devices::get_millis();
        double delta_t = fixed_dt;


        if (IO::input_devices::up_key())
        {

            if (millis>millis_last_move_frame+ms_per_frame)
            {
                ++frame;
                millis_last_move_frame=millis;
                frame%=steps;
            }
        }
        else if (IO::input_devices::down_key())
        {

            if (millis>millis_last_move_frame+ms_per_frame)
            {
                if (frame>0)
                    --frame;
                else
                    frame=steps-1;
                millis_last_move_frame=millis;
            }
        }
        if (IO::input_devices::right_key_click())
        {
            ++frame;
            frame%=steps;
        }
        if (IO::input_devices::left_key_click())
        {
            if (frame>0)
                --frame;
            else
                frame=steps-1;
        }

        cout<<"Frame : "<<frame<<" | key state "<<k[frame] <<" f="<<f[frame]<<endl;



        IO::graphics::draw_data(q[frame],min_q,max_q);

        try
        {
            IO::post_loop(millis+fixed_dt_millis );
        }
        catch(std::exception& Oof)
        {
            cout<<"Error in post loop: "<<Oof.what()<<endl;
            break;
        }
        catch(...)
        {
            cout<<"Error in post loop"<<endl;
            break;
        }

        //Do not exceed max refresh rate



        pmillis=millis;
    }while (!IO::input_devices::should_quit());
    cout<<"SHUTDOWN"<<endl;
    IO::end();


    return 0;
    //All is gone, the program has quit.
}
