/*THIS FILE
Definition for IO:input_devices functions
//This includes: mouse, keyboard and the text-input from the developer commandline (While these things sound separate, SDL treats them the same way, so separating them would not be practical)
*/


//Main sdl
#include <SDL2/SDL.h>


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

#include "IO.hpp"

//I use these types so much that these aliases are well worth it
using uchar= uint8_t;
using uint = uint32_t;
using ulong = uint64_t;

using namespace std;

namespace fs = std::filesystem;



//The SDL keyboard codes include things like up, down left right, etc. there are too many of these to fit into the normal 8-bit char, so SDL uses int, conveniently all regular letters (in the latin alphabet, not extended unicode) have the same number as the numerical value of them as chars. I like to rename this to Char so that I remember this is a keyboard character
using Char = int;

namespace IO::input_devices
{
    //--Internal variables and flags--
    bool quit=true;

    //milli seconds since start of SDL
    ulong millis=0;

    //Mouse down? l/r = left or right,
    bool l_mouse_down = false;
    bool r_mouse_down = false;
    bool l_p_mouse_down = false;//p indicates previous loop, this is because in some cases we don't want to spam certain actions just because the key is held down.
    bool r_p_mouse_down = false;
    int mouse_x, mouse_y;//mouse location, must be signed integers as this is what SDL sends them out as
    int scroll = 0;


    //Whether or not the selected keys are pressed
    bool up_press     = false;
    bool down_press   = false;
    bool left_press   = false;
    bool right_press  = false;
    bool A_press   = false;
    bool B_press  = false;
    bool C_press   = false;
    bool D_press  = false;
    bool E_press   = false;
    bool F_press  = false;
    bool tab_press = false;
    bool ctrl_press   = false;
    bool shift_press  = false;
    bool space_press  = false;
    bool esc_press  = false;

    //whether or not the keys were pressed last loop
    bool pup_press    = false;
    bool pdown_press  = false;
    bool pleft_press  = false;
    bool pright_press = false;
    bool pctrl_press  = false;
    bool pshift_press = false;
    bool pspace_press = false;
    bool pA_press   = false;
    bool pB_press  = false;
    bool pC_press   = false;
    bool pD_press  = false;
    bool pE_press   = false;
    bool pF_press  = false;
    bool ptab_press = false;
    bool pesc_press  = false;


    //Default keybindings, can (should) be changed
    Char UP_KEY         = SDLK_UP;
    Char LEFT_KEY       = SDLK_LEFT;
    Char RIGHT_KEY      = SDLK_RIGHT;
    Char DOWN_KEY       = SDLK_DOWN;
    Char TAB_KEY        = SDLK_LSHIFT;
    Char ESC_KEY        = SDLK_ESCAPE;
    Char SPACE_KEY      = SDLK_SPACE;
    Char ENTER_KEY      = SDLK_RETURN;
    Char DELETE_KEY     = SDLK_BACKSPACE ;

    Char A_KEY          = 'q';
    Char B_KEY          = 'w';
    Char C_KEY          = 'e';
    Char D_KEY          = 'a';
    Char E_KEY          = 's';
    Char F_KEY          = 'd';
    Char inventory_KEY  = 'i';
    fs::path keymaps;//The file containing the keyboard map

    void reload_keys();


    // ---- one time functions ----
    void updateMouse();//pre declare these things
    void init(fs::path _keymaps)
    {
        //We are not stopped
        quit = false;

        keymaps=_keymaps;
        reload_keys();


        updateMouse();
    }
    void end()
    {
    }

    void reload_keys()
    {
        ifstream keys(keymaps);

        if (!keys.is_open())
        {
            throw std::runtime_error("File "+keymaps.string()+" could not be opened");
        }

        //TBD, really, intead of having a number of pre-defined allowed keys, could we define more keys without recompiling?

        string input;
        //The keys and their SDL key code are stored in one long list, this is so int that plain-text loading is fine
        while (getline(keys,input))//Read one line at the time
        {
            //Skip comments
            if (input[0]=='#')
                continue;

            //Now stream the line content word for word
            stringstream ss(input);
            string key="";
            ss >> key;

            if (key.compare("LEFT")==0)
                ss>>LEFT_KEY;
            if (key.compare("RIGHT")==0)
                ss>>RIGHT_KEY;
            if (key.compare("DOWN")==0)
                ss>>DOWN_KEY;
            if (key.compare("UP")==0)
                ss>>UP_KEY;
            if (key.compare("TAB")==0)
                ss>>TAB_KEY;
            if (key.compare("A")==0)
                ss>>A_KEY;
            if (key.compare("B")==0)
                ss>>B_KEY;
            if (key.compare("C")==0)
                ss>>C_KEY;
            if (key.compare("D")==0)
                ss>>D_KEY;
            if (key.compare("E")==0)
                ss>>E_KEY;
            if (key.compare("F")==0)
                ss>>F_KEY;
            if (key.compare("ESC")==0)
                ss>>ESC_KEY;
            if (key.compare("SPACE")==0)
                ss>>SPACE_KEY;
            if (key.compare("ENTER")==0)
                ss>>ENTER_KEY;
            if (key.compare("DELETE")==0)
                ss>>DELETE_KEY;

        }
        keys.close();
    }


    //USE WAIT EXTREME CAUTION WILL FREEZE THE THREAD!
    void freeze_until(ulong time)
    {
        while (SDL_GetTicks()<time)
        {/*Just looping around in circles*/}
    }
    // ---- Functions which must be called once per loop ----
    //Polling and updating before and after ensures we can keep track of not only what is being done, but whether something has just changed
    void poll_events()
    {
        millis = SDL_GetTicks();

        //Save our old events
        pup_press    = up_press;
        pdown_press  = down_press;
        pleft_press  = left_press;
        pright_press = right_press;
        pA_press = A_press;
        pB_press = B_press;
        pC_press = C_press;
        pD_press = D_press;
        pE_press = E_press;
        pF_press = F_press;
        ptab_press = tab_press;
        pctrl_press  = ctrl_press;
        pshift_press = shift_press;
        pspace_press = space_press;
        pesc_press = esc_press ;

        r_p_mouse_down = r_mouse_down;
        l_p_mouse_down = l_mouse_down;

        scroll = 0;


        //Some things are gotten as keyboard mod states, i.e. is the effect of shift present, this gets both caps-lock, left and right shift
        shift_press = (SDL_GetModState() & KMOD_SHIFT);
        ctrl_press = (SDL_GetModState() & KMOD_CTRL);

        SDL_Event e;
        //Look at every single event, SDL queues events, so we can read more things happening at the same time
        while (SDL_PollEvent(&e) != 0)
        {
            //Match against those we actually care about

            switch (e.type)
            {
            //Death to the application!
            case SDL_QUIT:
                quit = true;
                break;

            //Keydown is only the act of pressing down this thing, releasing it is another act
            case  SDL_KEYDOWN:
                {
                    //Switch don't work with non-const statements for some reason, so use if-else instead... it is likely going to be interpreted the same way by the compiler anyway
                    if (e.key.keysym.sym == UP_KEY)
                    {
                        up_press = true;
                    }//once again, we can detect multiple keys at the same time, they just give two events in the while loop, each event can only be one thing
                    if (e.key.keysym.sym == DOWN_KEY)
                    {
                        down_press = true;
                    }
                    if (e.key.keysym.sym == LEFT_KEY)
                    {
                        left_press = true;
                    }
                    if (e.key.keysym.sym == RIGHT_KEY)
                    {
                        right_press = true;
                    }
                    if (e.key.keysym.sym == A_KEY)
                    {
                        A_press = true;
                    }
                    if (e.key.keysym.sym == B_KEY)
                    {
                        B_press = true;
                    }
                    if (e.key.keysym.sym == C_KEY)
                    {
                        C_press = true;
                    }
                    if (e.key.keysym.sym == D_KEY)
                    {
                        D_press = true;
                    }
                    if (e.key.keysym.sym == E_KEY)
                    {
                        E_press = true;
                    }
                    if (e.key.keysym.sym == F_KEY)
                    {
                        F_press = true;
                    }
                    if (e.key.keysym.sym == TAB_KEY)
                    {
                        tab_press = true;
                    }
                    if (e.key.keysym.sym == ESC_KEY)
                    {
                        esc_press = true;
                    }
                    if (e.key.keysym.sym == SPACE_KEY)
                    {
                        space_press = true;
                    }
                }
                break;

            case SDL_MOUSEBUTTONDOWN:
                if (e.button.button == SDL_BUTTON_LEFT)
                    l_mouse_down = true;
                else if (e.button.button == SDL_BUTTON_RIGHT)
                    r_mouse_down = true;
                //Also make sure mouse position is current
                updateMouse();

                break;

            case SDL_MOUSEBUTTONUP:
                if (e.button.button == SDL_BUTTON_LEFT)
                    l_mouse_down = false;
                else if (e.button.button == SDL_BUTTON_RIGHT)
                    r_mouse_down = false;
                //Also make sure mouse position is current
                updateMouse();
                break;

            case SDL_MOUSEMOTION:
                //Get mouse position
                updateMouse();

                break;

            //key just released
            case  SDL_KEYUP:
                {
                    if (e.key.keysym.sym == UP_KEY)
                    {
                        up_press = false;
                    }
                    if (e.key.keysym.sym == DOWN_KEY)
                    {
                        down_press = false;
                    }
                    if (e.key.keysym.sym == LEFT_KEY)
                    {
                        left_press = false;
                    }
                    if (e.key.keysym.sym == RIGHT_KEY)
                    {
                        right_press = false;
                    }
                    if (e.key.keysym.sym == A_KEY)
                    {
                        A_press = false;
                    }
                    if (e.key.keysym.sym == B_KEY)
                    {
                        B_press = false;
                    }
                    if (e.key.keysym.sym == C_KEY)
                    {
                        C_press = false;
                    }
                    if (e.key.keysym.sym == D_KEY)
                    {
                        D_press = false;
                    }
                    if (e.key.keysym.sym == E_KEY)
                    {
                        E_press = false;
                    }
                    if (e.key.keysym.sym == F_KEY)
                    {
                        F_press = false;
                    }
                    if (e.key.keysym.sym == TAB_KEY)
                    {
                        tab_press = false;
                    }
                    if (e.key.keysym.sym == ESC_KEY)
                    {
                        esc_press = false;
                    }
                    if (e.key.keysym.sym == SPACE_KEY)
                    {
                        space_press = false;
                    }

                }
                break;

            case SDL_MOUSEWHEEL://How scrolled is it
                if (e.wheel.y != 0) // scroll up or down
                {
                    scroll = e.wheel.y;
                }
                break;
            }//end switch (e.type)
        }//end while (SDL_PollEvent(&e) != 0)
    }



    void updateMouse()
    {
        int w_mouse_x,w_mouse_y;
        //Get mouse position, in a 2D world that is fairly easy
        SDL_GetMouseState(&w_mouse_x, &w_mouse_y);
        //These are pixels on the window, but 0,0 on the screen is not 0,0 internally, and the pixels are way smaller

        //Just read wherever the mouse is now
        mouse_x = w_mouse_x ;
        mouse_y = w_mouse_y ;
    }

    // ---- input functions ----
    ulong get_millis()//Get current milli seconds since the graphical display started (Measured at poll events)
    {
        return millis;
    }

    //--- Mouse ---
    void get_mouse_px(int& x, int& y)
    {
        x = mouse_x;
        y = mouse_y;
    }


    bool mouse_click(bool right)
    {
        if (right)
            return r_mouse_down && !r_p_mouse_down;
        else
            return l_mouse_down && !l_p_mouse_down;
    }
    bool mouse_press(bool right)
    {
        if (right)
            return r_mouse_down;
        else
            return l_mouse_down;

    }
    //Click is only true for one loop (useful for menu buttons), mouse press is true as long as the mouse is held down. Here is an example:
    /*
    Mouse button held down? no no no yes yes yes no yes yes
    mouse_click             0  0  0  1   0   0   0  1   0
    mouse_press             0  0  0  1   1   1   0  1   1
    */
    //SDL only allows integer scroll "steps" may be positive or negative
    int get_scroll()
    {
        return scroll;//I hope this is the scroll you wanted
    }


    //--- Keyboard keys held down ---

    //These keys are pre-defined keys we might want to look out for.
    //TBD, as the functions are defined, remapping is possible, but this should be possible without re-compiling

    //TBD, it seems -- wasteful, to define unique functions for each key, would it be possible to define a "key" class where this key then could map to whatever keys we need
    // Is a key currently held down?
    bool up_key()   { return up_press; }
    bool down_key() { return down_press; }
    bool left_key() { return left_press; }
    bool right_key(){ return right_press; }
    bool ctrl_key() { return ctrl_press; }
    bool shift_key(){ return shift_press; }
    bool space_key(){ return space_press; }
    bool esc_key()  { return esc_press; }
    bool A_key()    { return A_press; }
    bool B_key()    { return B_press; }
    bool C_key()    { return C_press; }
    bool D_key()    { return D_press; }
    bool E_key()    { return E_press; }
    bool F_key()    { return F_press; }
    bool tab_key()  { return tab_press; }
    bool should_quit(){return quit;}//Has the little x at the top of the window been pressed? Technically a button so it is here
    void do_quit(){quit=true;}

    //Was a key just clicked (true) or just released (false)
    bool A_key_click(bool click)
    {
        return (click && A_press && !pA_press) || (!click && !A_press && pA_press);
    }
    bool B_key_click(bool click)
    {
        return (click && B_press && !pB_press) || (!click && !B_press && pB_press);
    }
    bool C_key_click(bool click)
    {
        return (click && C_press && !pC_press) || (!click && !C_press && pC_press);
    }
    bool D_key_click(bool click)
    {
        return (click && D_press && !pD_press) || (!click && !D_press && pD_press);
    }
    bool E_key_click(bool click)
    {
        return (click && E_press && !pE_press) || (!click && !E_press && pE_press);
    }
    bool F_key_click(bool click)
    {
        return (click && F_press && !pF_press) || (!click && !F_press && pF_press);
    }
    bool tab_key_click(bool click)
    {
        return (click && tab_press && !ptab_press) || (!click && !tab_press && ptab_press);
    }
    bool down_key_click(bool click)
    {
        return (click && down_press && !pdown_press) || (!click && !down_press && pdown_press);
    }
    bool right_key_click(bool click)
    {
        return (click && right_press && !pright_press) || (!click && !right_press && pright_press);
    }
    bool left_key_click(bool click)
    {
        return (click && left_press && !pleft_press) || (!click && !left_press && pleft_press);
    }
    bool up_key_click(bool click)
    {
        return (click && up_press && !pup_press) || (!click && !up_press && pup_press);
    }
    bool ctrl_key_click(bool click)
    {
        return (click && ctrl_press && !pctrl_press) || (!click && !ctrl_press && pctrl_press);
    }
    bool shift_key_click(bool click)
    {
        return (click && shift_press && !pshift_press) || (!click && !shift_press && pshift_press);
    }
    bool space_key_click(bool click)
    {
        return (click && space_press && !pspace_press) || (!click && !space_press && pspace_press);
    }
    bool esc_key_click(bool click)
    {
        return (click && esc_press && !pesc_press) || (!click && !esc_press && pesc_press);
    }

}
