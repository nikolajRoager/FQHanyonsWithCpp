//Only include once, you can write pragma once, or you can use
#pragma once

/*THIS FILE
The public graphics, sound and control functions
These are put in the same header, so that we only need to include one header, the definitions are separate files though
*/

#define TWO_PI 6.283185307179586

#include<string>
#include<vector>
#include<cstdint>

//operating independent file-system functions
#include<filesystem>

//Vector math
#include <glm/glm.hpp>

#include <SDL2/SDL.h>

//TBD HIDE THIS
#include <GL/glew.h>

//I use these types so much that these aliases are well worth it
using tex_index = uint32_t;
using uint = uint32_t;
using ulong = uint64_t;

using namespace std;
using namespace glm;

namespace fs = std::filesystem;

//To capture taps and arrowkeys we need more than 8 bit, and it feels wrong to refer to it as an int
using Char = int;


namespace IO
{
    //Watch me struggle to put these functions in any sort of logical order:  ---- indicate a new chapter ----  ==== indicates a larger chapter =====

    //==== Global functions ====

    // ---- One time functions, which change some properties ----
    //This calls all the associated init functions IN THE RIGHT ORDER, which is important
    void init(string window_header,uint w, uint h, fs::path tex, fs::path audio, fs::path shaders,fs::path keymaps);
    void end();

    // ---- Functions which must be called once per loop ----

    void post_loop(ulong wait_until=0/*0=will not wait*/);//Actually display all the things we made ready to display, and update events. This ALSO Updates window, polls input (keyboard mouse) and resets the image for the next frame, must be called each loop


    //---- Debug functions ----
    double get_mspf();
    void print_ram();

    //==== Graphics =====
    namespace graphics
    {

        // ----  Display functions ----

        // ---- Functions for moving the camera around ----

        void add_rotation(float angle);
        void move_camera_origin(vec2 origin_offset);
        void multiply_zoom(float factor);


        //=== Get functions, yeah get_windows should not be there, it breaks encapsulation, sorry about that
        SDL_Window * get_window();
        uint get_w();
        uint get_h();

        uint load_hex_tex(const string& name);


        vec2 pixel_to_worldspace(uint x, uint y);

        //Get edges of all within the camera

        //--Mesh rendering functions--

        //Prepare to render to the lightmap or default display

        //Single hex, with planet to worldspace and local hex coordinates already given
        void draw_hex( const mat3& planetspace_to_worldspace, vec2 hex_coord, uint tex);

        void draw_circle( const vec2& offset,const vec4& colour);

        void draw_data(const vector<float>& Data,float min, float max);


        void draw_hex( const vec2& offset, uint tex);
        void draw_lines(const vector<glm::vec2>& verticies, uint size,glm::vec4 color,vec2 offset= vec2(0));
        void draw_triangles_simple(const vector<glm::vec2>& vertices, uint size,glm::vec4 color,vec2 offset= vec2(0));
        void draw_triangles(const vector<glm::vec2>& vertices, uint size,glm::vec4 color,glm::vec2 origin,float range=-1,vec2 offset= vec2(0));
        void draw_segments(const vector<glm::vec2>& vertices, uint size,glm::vec4 color,vec2 offset= vec2(0));

        void move_camera_origin(vec2 origin_offset);
        void add_rotation(float angle);
    }

    namespace audio
    {


        // --- Sound loading ---
        tex_index load_sound(string name);//Load this sound effect, and return its ID
        void delete_sound(tex_index ID);//Same procedures as with textures
        void add_sound_user(tex_index ID);

        // ---- Play sound functions ----
        void play_sound(tex_index ID);//"Display" this sound effect
        void loop_sound(tex_index ID);//Keep looping this sound as long as this is called every frame, for practical reasons we can keep at most 16 sounds looping at once

        //---Print ram ----
        void print_ram();//Print currently loaded sound effects to the external terminal (usually doesn't fit in the developer commandline)


    }




    //This includes: mouse, keyboard and the text-input from the developer commandline (While these things sound separate, SDL treats them the same way, so separating them would not be practical)
    namespace input_devices
    {


        // ---- input functions ----
        ulong get_millis();//Get current milli seconds since the graphical display started (Measured at poll events)

        //USE WAIT EXTREME CAUTION WILL FREEZE THE THREAD!
        void freeze_until(ulong time);

        //--- Mouse ---
        //Get position in pixels
        //Really this is what function SHOULD have been used
        //void get_mouse(int& x, int& y);
        void get_mouse_px(int& x, int& y);

        bool mouse_click(bool right = false);//Was the mouse just clicked
        bool mouse_press(bool right = false);//Is the mouse being held down

        //Click is only true for one loop (useful for menu buttons), mouse press is true as long as the mouse is held down. Here is an example:
        /*
        Mouse button held down? no no no yes yes yes no yes yes
        mouse_click             0  0  0  1   0   0   0  1   0
        mouse_press             0  0  0  1   1   1   0  1   1
        */
        //SDL only allows integer scroll "steps" may be positive or negative
        int get_scroll();


        //--- Keyboard keys held down ---

        //These keys are pre-defined keys we might want to look out for.

        //TBD, it seems -- wasteful, to define unique functions for each key, would it be possible to define a "key" class where this key then could map to whatever keys we need
        // Is a key currently held down?
        bool up_key();
        bool down_key();
        bool left_key();
        bool right_key();
        bool A_key();
        bool B_key();
        bool C_key();
        bool D_key();
        bool E_key();
        bool F_key();
        bool tap_key();
        bool ctrl_key();
        bool shift_key();
        bool space_key();
        bool esc_key();
        bool should_quit();//Has the little x at the top of the window been pressed? Technically a button so it is here
        void do_quit();//Set the above flag to true (maybe this current point in time is an inopportune moment to quite, and maybe it is easier just to wait until next time we check the flag)

        //Was a key just clicked (true) or just released (false)
        bool up_key_click(bool click=true);
        bool down_key_click(bool click=true);
        bool left_key_click(bool click=true);
        bool right_key_click(bool click=true);
        bool ctrl_key_click(bool click=true);
        bool shift_key_click(bool click=true);
        bool space_key_click(bool click=true);
        bool A_key_click(bool click=true);
        bool B_key_click(bool click=true);
        bool C_key_click(bool click=true);
        bool D_key_click(bool click=true);
        bool E_key_click(bool click=true);
        bool F_key_click(bool click=true);
        bool tab_key_click(bool click=true);
        bool esc_key_click(bool click=true);

    }


}
