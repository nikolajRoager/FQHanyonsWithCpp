/*THIS FILE
Definition to the master IO functions (those which secretly call the other namespaces behind the scenes)
These are put in the same header, so that we only need to include one header, the definitions are separate files though
*/


//openGL Extension Wrangler library
#include <GL/glew.h>
//Main sdl
#include <SDL2/SDL.h>
//loading of png and jpg
#include <SDL2/SDL_image.h>
//openGL functionalities
#include <SDL2/SDL_opengl.h>
//sound
#include <SDL2/SDL_mixer.h>


//std string
#include <string>
//string-stream, for turning a string into a "stream" object, useful for loading plain text files
#include <sstream>
//text output, should never be used in looping functions
#include <iostream>
//operating independent file-system functions
#include<filesystem>
//file streams, we use filesystem for paths, and fstream to actually stream the data from files in.
#include <fstream>


//Dynamic sized list
#include <vector>

//Not so dynamic sized lists
#include <array>


//I use these types so much that these aliases are well worth it
using tex_index = uint32_t;
using uint = uint32_t;
using ulong = uint64_t;

using namespace std;

namespace fs = std::filesystem;

//The separate init or end functions should not be visible for anyone who includes IO.hpp, but I need to see them here, so declare them here, the linker will link them to their definitions
//I don't think this is great style, but I do not see any usable alternative


namespace IO::graphics
{
    // ---- one time functions ----
    void init(string window_header,uint w, uint h, fs::path tex,  fs::path shaders);
    void end();

        // ---- Functions which must be called once per loop ----
        void clear();//Updates window  and resets the image, must be called each loop
        void flush();//Actually display all the things we made ready to display



        //---Functions for getting information about the graphical setup ar currently not public---
}
namespace IO::audio
{
    // ---- one time functions ----
    void init(fs::path audio);
    void end();
    // ---- Functions which must be called once per loop ----
    void update_loops();//Update looping audio
}

namespace IO::input_devices
{
    // ---- one time functions ----
    void init(fs::path keymaps);
    void end();
    // ---- Functions which must be called once per loop ----
    //Polling and updating before and after ensures we can keep track of not only what is being done, but whether something has just changed
    void poll_events();
    //NOTE all inputs are updated at poll_events. Hence calling any of these input functions returns the same result regardless of when and how many times you call it
}

#include "IO.hpp"


namespace IO
{
    // Some internal variables
    //All the paths we need to know
    fs::path texture_path;
    fs::path sound_path;
    fs::path shader_path;



    double avg_mspf=0;
    uint frames_this_second = 0;
    uint last_second=0;


    void init( string window_header,uint w, uint h, fs::path tex, fs::path audio, fs::path shaders,fs::path keymaps)
    {
        texture_path = tex;
        sound_path = audio;
        shader_path = shaders;

        //If this is a restart command, actually just stop and start again
        end();


        //Get the error ready, just in case
        string error = "";

        //Initialize SDL globally, it does not make sense to include in any specific subspace
        if (SDL_Init(SDL_INIT_EVERYTHING) < 0)
        {
            error.append("Could not initialize SDL:\n\t\t");
            error.append(SDL_GetError());
            throw std::runtime_error(error);
        }

        graphics::init(window_header,w,h, tex, shaders);


        audio::init(audio);

        input_devices::init(keymaps);

        last_second = input_devices::get_millis();

        avg_mspf=0;
        frames_this_second = 0;
        last_second=0;

        //Do one false loop, this is just to make the screen all black while we load
        graphics::clear();
        post_loop();
    }

    double get_mspf()
    {
        return avg_mspf;
    }



    bool should_quit()
    {
        //The little x on the top of the screen has been pressed, the program can also be quit if this is not set
        return input_devices::should_quit();//The poor program, it does not know that returning true will have it murdered instantly
    }

    //the part where all we have done this frame actually ends up on the screen, this function both posts the now finished frame to the screen, and makes ready for the next frame
    void post_loop(ulong wait_until)
    {
        ++frames_this_second;

        ulong millis = input_devices::get_millis();
        if (millis >last_second+1000)
        {
            avg_mspf=double(millis-last_second)/frames_this_second;

            //cout <<"Cycle time: τ="<<avg_mspf<< " ms; ν = "<<(1000/avg_mspf)<<" Hz"<<endl;
            last_second = millis;
            frames_this_second = 0;
        }

        //Now update looping soun, it really does not matter if we do this before or after the graphics is flushed.
        audio::update_loops();
        //Any drawing until now have been done to the internal framebuffer, now flush this onto the screen, after this point and until we call clear the drawing functions do not work.
        graphics::flush();

        //TBD, Opengl VSync is not reliable on all machines, implement a custom waiting function to not exceed 60 Hz. ... right here is where we need it


        graphics::clear();//Alright, now we have officially started next frame, doing this in one function ensures we never have any draw functions between flush and clear
        input_devices::freeze_until(wait_until);
        input_devices::poll_events();

    }



    void print_ram()
    {
        cout<<"===Texture and Audio memory report==="<<endl;
        audio::print_ram();
        cout<<"====================================="<<endl;
    }

    void end()
    {
        input_devices::end();
        graphics::end();
        audio::end();

        //Close SDL
        SDL_Quit();
    }

}//ends namespace IO
