/*THIS FILE

The definition for the IO::graphics namespace, which handles the 2D pixelated texture rendering. This namespace also handles text-box rendering, this might sound like something completely different which need to be given its own namespace, however text boxes are also textures, and openGL functionalities are used to generate them, so all attempts at separating this into its own namespace have not succeeded

General note, OpenGL is NOT object oriented, the OpenGL functions I use here are actually regular c functions.
*/

//Needed to unlock some more glm features

//openGL Extension Wrangler library
#include <GL/glew.h>
//Main SDL library
#include <SDL2/SDL.h>
//loading of png and jpg
#include <SDL2/SDL_image.h>
//openGL functionalities
#include <SDL2/SDL_opengl.h>

//Vector math
#include <glm/glm.hpp>
//Matrix math: gtx require experimental flag


using namespace glm;


#include "IO.hpp"
#include "load_shader.hpp"
//A texture wrapper, which includes not only the graphical information but also the name and path required to reload this texture (which may also be a text-box)

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


namespace IO::graphics
{
    //Watch me continue to struggle defining the internal graphics variable in any sort of sensible order

    // ---- OpenGL setup and renderbuffers ----
    //Just one window for now ... actually ever, turns out SDL does not handle more than one window well.
    SDL_Window* window = nullptr;
    //Another thing which we only need to define one of ... it does important things behind the scene
    SDL_GLContext context;

    GLuint VertexArrayObject = -1;//I will create this, but I won't use it for anything actually, I just need to create it

    //In openGL we render to "framebuffers" which include both the rendered image and more buffers (depth buffers), framebuffer 0 always refers to the main screen

    //The main display is rendered to a texture for post processing and only rendered to the actual screen at the end, hence why I call this "display"
    GLuint display_Framebuffer = -1;
    GLuint display_texture = -1;


    //window width/height in external window pixels
    uint win_h;
    uint win_w;

    //lattice width height
    uint w;
    uint h;

    vector<float> data;
    float Min=0;
    float Max=1.0;



    // ---- Other internal flags and data which we need to keep around ----
    //That is, things the graphics needs to remember between functions.
    bool rsz_scr = false;//Flag to tell us that the size of the window just changed, we may need to do a lot of updates later down the line, but these updates are relatively slow, so only do them if we need

    //Paths we load things from
    fs::path texture_path;
    fs::path shader_path;

    const mat2 hexcoord_to_cartesian=mat2(
        sqrt(3), 0,
        -sqrt(3)/2 , 1.5
    );


    const mat3 unit_matrix_ref = mat3(1);//An identity matrix, this will be the default matrix used for non-animated textures:
    /*
    1 0 0
    0 1 0
    0 0 1
    */

    //---- The textured rectangle shader ----
    //A simple textured rectangle
    GLuint surf_ProgramID = -1;

    //The singular attribute, the position of the rectangle
    GLuint surf_VertexPosAttribID = -1;

    GLuint surf_PosBuffer = -1;

    //A colored circle
    GLuint circle_ProgramID = -1;

    //The singular attribute, the position of the rectangle
    GLuint circle_VertexPosAttribID = -1;

    GLuint circle_baseColour_ID = -1;//The base colour
    GLuint circle_r_ID = -1;//The radius
    GLuint circle_matrix_ID=-1;//The tranformation matrix from whatever coordinates are given to screenspace coordinates in the surf program



    GLuint hex_PosBuffer = -1;


    //These IDs are where the uniforms are stored on the graphics card, so that we can upload data to them
    GLuint surf_colorTex_ID = -1;//The identity of the texture in the surf program
    GLuint surf_W_ID= -1;
    GLuint surf_H_ID= -1;

    GLuint surf_min_ID= -1;
    GLuint surf_max_ID= -1;


    GLuint surf_Q_ID= -1;

    GLuint surf_matrix_ID=-1;//The tranformation matrix from whatever coordinates are given to screenspace coordinates in the surf program
    GLuint surf_tex_matrix_ID=-1;//Transformation of the texture, for instance for animation


    // ----- The Line shader -------
    //A line which is the basis of the mesh-rendering engine
    GLuint Line_ProgramID = -1;

    //The vertex position of the lines drawn, in world-space
    GLuint Line_VertexPosAttribID = -1;

    GLuint Line_matrix_ID=-1;//The tranformation matrix from line-worldspace coordinates to device coordinates in the Line program
    GLuint Line_color_ID=-1; //The color of the line to be drawn


    vector<GLuint> hexagon_textures;


    float px_per_m = 50;
    vec2 origin(0,0);
    float camera_rotation = 0;


    void reset_screenspace_matrix();

    void move_camera_origin(vec2 origin_offset)
    {

        float cosR = cos(camera_rotation);
        float sinR = sin(camera_rotation);

        mat2 inv_rotate( cosR,sinR,
                        -sinR,cosR);



        origin+=inv_rotate*origin_offset;



        reset_screenspace_matrix();
    }

    void add_rotation(float angle)
    {
        camera_rotation+=angle;
        reset_screenspace_matrix();
    }
    void multiply_zoom(float factor)
    {
        px_per_m/=factor;
        if (px_per_m<10)
            px_per_m=10;
        else if (px_per_m>200)
            px_per_m = 200;
        reset_screenspace_matrix();
    }



    float screen_width_per_meter= px_per_m/win_w;
    float screen_height_per_meter= px_per_m/win_h;

    mat3 worldspace_to_screenspace =
            mat3(vec3(screen_width_per_meter,0,0),
                 vec3(0,screen_height_per_meter,0),
                 vec3(0,0,1)
                );


    mat3 screenspace_to_worldspace =
            mat3(vec3(1.0/screen_width_per_meter,0,0),
                 vec3(0,1.0/screen_height_per_meter,0),
                 vec3(0,0,1)
                );


        mat3 hextex_matrix =
        mat3(
            vec3(1.0/sqrt(3.f), 0, 0),
            vec3(0,  1/3.f, 0),
            vec3(0.0, 1/3.f+0.5, 1)
        );


    //Pre-declare some things we need in init
    void set_output_size(uint _w, uint _h);


    GLuint load_tex(const fs::path& path);

    // ---- one time functions ----
    void init(string window_header,uint W, uint H, fs::path tex,  fs::path shaders)

    {
        w=W;
        h=H;
        texture_path = tex;
        shader_path = shaders;

        //Get the error ready, just in case
        string error = "";

        cout<<"Init openGL"<<endl;

        //Initialize openGL version, this version should be good enough, it is not the newest, so some crazy experimental features are not there.
        SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 4);
        SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);
        SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);

        //Note this is a c style binary flag, the int here is just a list of 0's and 1's, each 1 corresponds to one thing we want to turn on, for instance this could be 0011 or something like that, the | operator is a binary or, (i.e. 0010 | 0001 = 0011)
        int imgFlags = IMG_INIT_JPG | IMG_INIT_PNG;
        cout<<"Init images"<<endl;
        if (!(IMG_Init(imgFlags) & imgFlags))//the & is a binary and, if this thing is all 0, then something was not initialized correctly
        {
            error.append(" Could not initialize SDL_image with JPEG and PNG loading:\n\t\t");
            error.append(IMG_GetError());
            throw std::runtime_error(error);
        }


        win_w = 1920;
        win_h = 1000;


        data=vector<float>(w*h);
        Min=0.0;
        Max=1.0;


        //The transformation for the surface used to display the rescaled screen texture, will be manually calculated later
        //Create window and its surface
        SDL_DisplayMode DM;

        cout<<"Update screen width height"<<endl;
        SDL_GetCurrentDisplayMode(0, &DM);
        win_w = DM.h;
        win_h = DM.h;

        cout << "Width " << win_w << " Height " << win_h << endl;

        uint FLAGS = SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN;

        cout<<"Open Window"<<endl;
        window = SDL_CreateWindow(window_header.c_str(), 0, 0, win_w, win_h, FLAGS);


        if (window == nullptr)
        {
            end();
            throw std::runtime_error(SDL_GetError());
        }

        cout<<"Init cursor"<<endl;
        SDL_ShowCursor(SDL_DISABLE);

        cout<<"Create OpenGL"<<endl;
        context = SDL_GL_CreateContext(window);

        if (context == nullptr)
        {
            error.append("OpenGL context failed to be creted\n\t\t");
            error.append(SDL_GetError());
            throw std::runtime_error(error);
        }

        glewExperimental = GL_TRUE;
        cout<<"Init Glew"<<endl;
        GLenum glewError = glewInit();

        if (glewError != GLEW_OK)
        {
            error.append("Glew not initialized OK\n\t\t");
            error.append(string((char*)glewGetErrorString(glewError)));
        }

        cout<<"Set up OpenGL settings"<<endl;
        //I have never actually needed to use this thing explicitly, and I am sure at least 90% of applications don't either, and until OpenGL 3.2 this would be defined behind the scenes but now we need to define it to get our programs to work... I would have prefered it it would just be done by OpenGL behind the scenes still, OpenGL just seems to require a little to many lines of setup already.
        glGenVertexArrays(1, &VertexArrayObject);
        glBindVertexArray(VertexArrayObject);


        //Use Vsync, this seems a little questionable on some computers
        if (SDL_GL_SetSwapInterval(1) < 0)
        {
            cerr<<"-----------------"<<endl;
            cerr<<"WARNING VSYNC NOT SUPPORTED (SDL2 RETURNED THE FOLLOWING ERROR "<<SDL_GetError()<<")\nThe game will try to run, but will not perform optimally.\nVsync is known to be unavailable on older graphics card, but it is much more likely that SOMETHING IS WRONG WITH THE DRIVERS FOR YOUR GRAPHICS CARD Check if they are installed correctly."<<endl;
            #if defined( _WIN32)
            cout<<"You are using windows (I can't see which version though), if you use Windows 10 or 11 you can check this in the \"device manager\" find your graphics card (may be named something like \"Display adapters\", right click and select properties, ALTERNATIVELY write the following command in the Windows powershell as administrator \"Get-WmiObject Win32_PnPSignedDriver| select DeviceName, Manufacturer, DriverVersion\""<<endl;//Ha ha good luck
            #elif defined(__APPLE__) && defined(__MACH__)
            cout<<"You use an Apple product, which really should come with drivers installed correctly from the factory."<<endl;
            #elif defined(__FreeBSD__) || defined(BSD)
            cout<<"You use a Unix like system (FreeBSD), to check what driver you use on, try running: lspci -k | grep -EA3 \"VGA|3D\""<<endl;
            #elif defined(__ANDROID__)
            cout<<"You use an android device, which really should come with drivers installed correctly from the factory."<<endl;
            #elif defined(__linux__)//else means that ANDROID, which is Linux, will not get the regular linux error message
            cout<<"You use Linux (Good) (I can't see what distro though), to check what driver you use on, try running: lspci -k | grep -EA3 \"VGA|3D\""<<endl;
            cout<<"Proprietary drivers can usually be installed through your package manager (pacman, apt-get etc.);\nNOTE FOR ARCH LINUX, at the time of writing there are some known issues with Nvidia graphics card and newer intel CPUs (gen 11), switching kernel from linux to linux-lts (run sudo pacman -S linux-lts nvidia-lts) might fix this, at the cost of downgrading the kernel; check wiki.archlinux.org/title/NVIDIA for up to date instructions."<<endl;
            #else
            cout<<"I am not able to recognize your OS, so I can't help you further"<<endl;
            #endif
            cerr<<"-----------------"<<endl;
        }

        //We have an open window, but nothing in it yet, this looks bad (it may include a random chunk of whatever was on the screen before the window was opened and it will look as if the program has crashed), it might take a few seconds to load so lets change the window to black

        //A most wonderfull background color
        glClearColor(0.0f, 0.0f, 0.0f, 1.f);
        //behind the scenes OpenGL has two buffers for the things to be shown to the window, we set one of them to black here
        glClear(GL_COLOR_BUFFER_BIT);
        //Then we swap the buffers, putting the black buffer into the window and we will keep it there until we have rendered the first frame to the other buffer
        SDL_GL_SwapWindow(window);

        //Alpha blending, aka transparant textures
        glEnable(GL_BLEND);
        glBlendEquation(GL_FUNC_ADD);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);//Good for texture blending (front texture obscures back)

        //In case this was on by default!
        glDisable(GL_DEPTH_TEST);
        //We are not working in 3D, so this is not needed (likely not on by default) this would only draw triangles facing the camera, 2D makes this irrelevant
        glDisable(GL_CULL_FACE);

        cout<<"\nLoading Shaders"<<endl;
        //Loading the rectangular surface shader

        string log;
        surf_ProgramID = load_program(shader_path, "surface", log);

        surf_VertexPosAttribID = glGetAttribLocation(surf_ProgramID, "vertex_pos");

        surf_matrix_ID = glGetUniformLocation(surf_ProgramID, "Input_to_DC");
        surf_tex_matrix_ID = glGetUniformLocation(surf_ProgramID, "UV_transform");
        surf_colorTex_ID = glGetUniformLocation(surf_ProgramID, "colorSampler");

    surf_W_ID= glGetUniformLocation(surf_ProgramID, "w");
    surf_H_ID= glGetUniformLocation(surf_ProgramID, "h");
    surf_min_ID= glGetUniformLocation(surf_ProgramID, "Min");
    surf_max_ID= glGetUniformLocation(surf_ProgramID, "Max");
    surf_Q_ID= glGetUniformLocation(surf_ProgramID, "Q");



        cout << "\nLoaded surface program" <<endl;

        cout << log << endl;

        circle_ProgramID = load_program(shader_path, "circle", log);

        circle_VertexPosAttribID = glGetAttribLocation(circle_ProgramID , "vertex_pos");

        circle_matrix_ID= glGetUniformLocation(circle_ProgramID , "Input_to_DC");
        circle_baseColour_ID = glGetUniformLocation(circle_ProgramID , "base_colour");
        circle_r_ID = glGetUniformLocation(circle_ProgramID , "radius");

        cout << "\nLoaded Circle program" <<endl;


        Line_ProgramID = load_program(shader_path, "line", log);

        Line_VertexPosAttribID = glGetAttribLocation(Line_ProgramID , "vertex_worldspace");

        Line_matrix_ID = glGetUniformLocation(Line_ProgramID , "worldspace_to_DC");
        Line_color_ID= glGetUniformLocation(Line_ProgramID , "color");

        glLineWidth(1.f);

        cout << "Loaded Line program " <<endl;
        cout << log << endl;





        static const GLfloat surface_vertex_data[] = {
            -1.0f, -1.0f,
             1.0f, -1.0f,
             1.0f,  1.0f,
            -1.0f,  1.0f
        };


        cout<<"Set up OpenGL buffer objects"<<endl;


        //Create VBOs

        //Load up the basic surface vertices into the GPU
        glGenBuffers(1, &surf_PosBuffer);
        glBindBuffer(GL_ARRAY_BUFFER, surf_PosBuffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(surface_vertex_data), surface_vertex_data, GL_STATIC_DRAW);
        //Static draw since we will not need to update this at any point



        //Load up the basic surface vertices into the GPU
        float cos30 = sqrt(3.f)/2.f;
        float sin30 = 0.5f;
        static const GLfloat hex_vertex_data[] = {
            cos30, sin30,
              0.0f, 1.0f,
            -cos30, sin30,
            -cos30,-sin30,
              0.0f,-1.0f,
             cos30,-sin30
        };


        glGenBuffers(1, &hex_PosBuffer);
        glBindBuffer(GL_ARRAY_BUFFER, hex_PosBuffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(hex_vertex_data), hex_vertex_data, GL_STATIC_DRAW);

        //Just check that starting width and height is what we set it to
        int tw, th;//this width and height I first used just w and h and got confused when the later setup of the framebuffer gave me waaay to high resolution
        SDL_GetWindowSize(window, &tw, &th);
        win_h = th;
        win_w = tw;

        //set the internal output size
        set_output_size(win_w, win_h);

        hexagon_textures=vector<GLuint>();

        cout<<"Init graphics done"<<endl;

    }

    uint load_hex_tex(const string& name)
    {

        hexagon_textures.push_back(load_tex(texture_path/string(name+".png")));

        return hexagon_textures.size()-1;
    }
    GLuint load_tex(const fs::path& path)
    {
        GLuint textureID=-1;
        //Load to a simple SDL texture; SDL knows what to do with any file format; it can itself figure out what format we are using
        SDL_Surface* texture_surface = nullptr;

        {
            texture_surface = IMG_Load(path.string().c_str());


            if (texture_surface == NULL)
            {
                cout <<"THERE WAS AN ERROR loading a texture, if \""<< path.string().c_str() << "\" is different from \""<<path.string()<<"\", then something has gone wrong in the port, and it is all Microsoft's fault; if not look at the below error message"<<endl;
                throw std::runtime_error("Couldn't load texture " + path.string() + " to surface: SDL returned error: " + string(IMG_GetError()));
            }

            uint texture_w = texture_surface->w;
            uint texture_h = texture_surface->h;

            //Now let us try to figure out what on earth this image is; if this is jpg we know that the image format is 8 bit GL_RGB, but if this is PNG it is anybodies guess. 16 bit png won't work, sorry.

            GLenum image_format;//What data is used by the image file
            GLenum GL_format;//And what should our target data be

            //Does this do Alpha?
            if (texture_surface->format->BytesPerPixel == 4)
            {
                image_format = GL_RGBA;
                GL_format = GL_RGBA;//I assume that if the image has an alpha channel, it is meant to be there
            }
            else
            {
                image_format = GL_RGB;
                GL_format = GL_RGB;
            }//I am going to be using png and jpg, formats which do RGB not BGR so I won't even consider that


            //This is our destination
            glGenTextures(1, &(textureID));
            glBindTexture(GL_TEXTURE_2D, textureID);
            glTexImage2D
            (
                GL_TEXTURE_2D,//What texture type is this
                0,//Mipmap-Level, this project uses low-resolution textures intentionally
                GL_format,//internal format for graphics card
                texture_surface->w,//Width/height in pixels
                texture_surface->h,
                0,
                image_format,       //Format of input
                GL_UNSIGNED_BYTE,//Type of data
                texture_surface->pixels //A pointer to the actual data,
            );

            //This is the place to glGenerateMipmap(GL_TEXTURE_2D);
            //But for the purpose of this project I do not want that

            //Now set some more parameters,
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

            //Get best alignment for pixel storage, it is the largest element in {1,2,4,8} which divides the pitch of the texture
            uchar align = 8;
            while (texture_surface->pitch % align) align >>= 1;
            //FASTER binary notation for align/=2, (shift every bit one right).
            //For a texture with power of two width/height greater than 8, I sure hope that align=8.

            glPixelStorei(GL_UNPACK_ALIGNMENT, align);
            //Doing this may result in a very small performance boost, I am not sure it matters though ... no wait, just checked, it literally does not work without

            SDL_FreeSurface(texture_surface);

        }

        return textureID;
    }


    void end()
    {
        //Destroy everything
        SDL_DestroyWindow(window);
        window = nullptr;

        for (GLuint textureID : hexagon_textures)
        {
            if (textureID != (GLuint)-1)
            {
                glDeleteTextures(1, &textureID);
                textureID = -1;
            }
        }


        //Explicitly unload all buffers, framebuffers and textures provided they are not already set to default not set value (-1)
        if (VertexArrayObject != (GLuint)-1)
        {
            glDeleteVertexArrays(1, &VertexArrayObject);
            VertexArrayObject = -1;
        }

        if (hex_PosBuffer != (GLuint)-1)
        {
            glDeleteBuffers(1, &hex_PosBuffer);
            hex_PosBuffer = -1;
        }
        if (surf_PosBuffer != (GLuint)-1)
        {
            glDeleteBuffers(1, &surf_PosBuffer);
            surf_PosBuffer = -1;
        }


        if (surf_ProgramID != (GLuint)-1)
        {
            glDeleteProgram(surf_ProgramID);
            surf_ProgramID = -1;
        }
        if (Line_ProgramID != (GLuint)-1)
        {
            glDeleteProgram(Line_ProgramID);
            Line_ProgramID = -1;
        }

        if (display_Framebuffer != (GLuint)-1)
        {
            glDeleteFramebuffers(1, &display_Framebuffer);
            display_Framebuffer = -1;
        }
        if (display_texture != (GLuint)-1)
        {
            glDeleteTextures(1, &display_texture);
            display_texture = -1;
        }

        //Delet eopenGL context
        if (context != nullptr)
        {
            SDL_GL_DeleteContext(context);
            context = nullptr;
        }


        //Quit image module
        IMG_Quit();

    }



    // ---- Functions which must be called once per loop ----
    void clear()//Updates window  and resets the image, must be called each loop
    {

        //Clear all framebuffers
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glClearColor(0.2f, 0.2f, 0.2f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT);

        glBindFramebuffer(GL_FRAMEBUFFER, display_Framebuffer);
        glClearColor(0.0f, 0.0f, 0.0f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT);

        glViewport(0, 0, win_w, win_h);


        //Test if width/height of window has changed
        int nw = 0;
        int nh = 0;
        SDL_GetWindowSize(window, &nw, &nh);

        if ((uint)nh != win_h || (uint)nw != win_w)
        {
            win_w = nw;
            win_h = nh;
            rsz_scr = true;

            set_output_size(win_w, win_h);

        }
        else
        {
            rsz_scr = false;
        }



    }

    void flush()
    {
        //Simply draw a full square to the screen, which will include our rendered texture, stretching this to fit the entire screen with closest interpolation results in the pixelation effect
        //Now render the output to the screen
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);//Good for texture blending (front texture obscures back)
        //glBlendFunc(GL_SRC_ALPHA, GL_DST_ALPHA);//Needed for light blending (true additive)
        //Enable what we need
        glUseProgram(surf_ProgramID);

        //Screen viewport, the only one which can change, but we need to update it every frame because OpenGL only stores one set of viewport variables globally, and not individual ones on a per framebuffer basis
        glViewport(0, 0, win_w, win_h);


        //OpenGL has a number of different texture "slots" here I activate slot number 0
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D,  display_texture);//Here I bind the texture I want to display, it is bound to the currently active texture slot (0 here)
        //1i tells opengl that this is 1 integer, 2f would be a 2-d vector with floating point data. It is very important that the type matches the data in the glsl shader.
        glUniform1i(surf_colorTex_ID, 0);//Here I set the uniform texture sampler to integer 0, OpenGL understand that this means it shall look in slot 0


        glUniform1i(surf_W_ID, w);
        glUniform1i(surf_H_ID, h);

        glUniform1f(surf_min_ID, Min);
        glUniform1f(surf_max_ID, Max);



        glUniform1fv(surf_Q_ID, 36, data.data());

        mat3 Screen_matrix =
            mat3(
                vec3(2.0, 0, 0),
                vec3(0,  2.0, 0),
                vec3(-1.0, -1.0, 1)
            );



        //OpenGL uniform MatrixN sends an N by N matrix, f says it has floating values and v tells us that we want to send a pointer to the values in the CPU memory, rather than setting the values manually
        glUniformMatrix3fv(
        surf_matrix_ID,
        1,                       //How many (OpenGL does support array uniforms), in this case there is only 1
        GL_FALSE,                //Should openGL transpose this for us, (transposing a matrix swaps element j,i into i,j), I don't see why this is useful enough to be a parameter here, I can just transpose the matrix myself before sending it to OpenGL
        &(Screen_matrix [0][0])   //Due to the way the glm library works, the matrix elements is stored in a list with 0,0 being the first. The order is correct for OpenGL to copy the data correctly if we tell it to start at the adress (&) of element 0,0, the paranthesis around the element is not needed, but I think it makes it more clear what order things happen in
        );
        //This is not an animation, send a unit matrix to "transform" the UV coordinates nowhere
        glUniformMatrix3fv(surf_tex_matrix_ID, 1, GL_FALSE, &(unit_matrix_ref[0][0]));//Same procedure as before


        //We only have one vertex attribute, these are the values which change from vertex to vertex
        glEnableVertexAttribArray(surf_VertexPosAttribID);
        glBindBuffer(GL_ARRAY_BUFFER, surf_PosBuffer);
        glVertexAttribPointer
        (
            surf_VertexPosAttribID,//Attribute location, you can either locate this in the program (with glGetAttributeLocation), or you can force the shader to use a particular location from inside the shader, I do the former
            2,                   //Number of the below unit per element (this is a 2D vector, so 2)
            GL_FLOAT,            //Unit, this is single precition GL float
            GL_FALSE,             //Normalized? Nope
            0,                    //Stride and offset, I don't use these for anything
            (void*)0
        );


        glDrawArrays(GL_TRIANGLE_FAN,0,4);

        //Disable everything activated
        glDisableVertexAttribArray(surf_VertexPosAttribID);



        glUseProgram(0);

        //Now we want to see this
        SDL_GL_SwapWindow(window);//behind the scenes OpenGL has two buffers for the things to be shown to the window, now we put the buffer we just edited into the window so the image changes, and in the next loop we edit the other buffer.

    }



    // ---- Functions for getting information about the graphical setup ----

    //Get width and height of internal display in pixels
    uint get_w()
    {
        return win_w;
    }
    uint get_h()
    {

        return win_h;
    }

    //It is possible to run at different settings, but changing outside init is not recommended
    //Set up the render-targets which are not the screens, this looks very convoluted because OpenGL gives a lot of options which we won't use in this case: I just want to render to a texture, but I create a Framebuffer, which is the thing I actually renders too, which the texture is attached to; actually Framebuffers can contain many textures or buffers to encode more information, such as depth ... or whatever really, OpenGL is an extremely general tool.
    void set_output_size(uint _w, uint _h)
    {
        //First, make sure that the Framebuffers and textures are empty.
        if (display_Framebuffer != (GLuint)-1)
        {
            glDeleteFramebuffers(1, &display_Framebuffer);
            display_Framebuffer = -1;
        }
        if (display_texture != (GLuint)-1)
        {
            glDeleteTextures(1, &display_texture);
            display_texture = -1;
        }


        //Set up the display we will render to
        glGenFramebuffers(1, &display_Framebuffer);
        glBindFramebuffer(GL_FRAMEBUFFER, display_Framebuffer);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);//Good for texture blending (front texture obscures back)

        glGenTextures(1, &display_texture);
        glBindTexture(GL_TEXTURE_2D, display_texture);
        //Initialize empty, and at the size of the internal screen
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, win_w, win_h, 0, GL_RGB, GL_UNSIGNED_BYTE, 0);

        //No interpolation
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);


        //Now the display framebuffer renders to the texture we finally display
        glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, display_texture, 0);
        GLenum DrawBuffers[1] = { GL_COLOR_ATTACHMENT0 };
        glDrawBuffers(1, DrawBuffers);//Color attachment 0 as before


        //Now we can render to whatever we want by using:
        //glBindFramebuffer(GL_FRAMEBUFFER,0);//Final output i.e. what the user to see this, but not used internally because we want low resolution for pixelation effects
        //glBindFramebuffer(GL_FRAMEBUFFER,display_Framebuffer);//Default rendertarget used here, renders result to a texture which we will display to screen (with framebuffer 0)
        //by the way I DID NOT FORGET glViewport, but it does not go in the setup, because OpenGL only stores one viewport size globally, it can not be set on a per-framebuffer basis (As one would expect), 1and since we are using several different sized viewports per display we will need to manyally call it twice per loop anyway (in case you wonder: no, you do not need to set viewport before calling glClear, many code examples on the internet do put their viewport call very early because they only use one size framebuffer)




        reset_screenspace_matrix();

    }


    void reset_screenspace_matrix()
    {
        screen_width_per_meter= px_per_m/win_w;
        screen_height_per_meter= px_per_m/win_h;


        mat3 scale(screen_width_per_meter,0,0,
                   0,screen_height_per_meter,0,
                   0,0                      ,1);

        mat3 translate( 1       ,0        ,0,
                        0       ,1        ,0,
                       -origin.x,-origin.y,1);
        float cosR = cos(camera_rotation);
        float sinR = sin(camera_rotation);

        mat3 rotate( cosR ,-sinR ,0,
                     sinR , cosR ,0,
                     0    ,0     ,1);

        worldspace_to_screenspace = scale*rotate*translate;
        screenspace_to_worldspace =inverse(worldspace_to_screenspace );

    }


    //---- OpenGL rendering functions used internally----

    SDL_Window * get_window()
    {
        return window ;
    }

    //-- Mesh rendering functions --


    void draw_unicolor(const vector<vec2>& vertices, uint size, vec4 color,GLuint displaymode,vec2 offset= vec2(0));

    void draw_lines(const vector<vec2>& vertices, uint size, vec4 color,vec2 offset)
    {
        draw_unicolor(vertices,size,color,GL_LINE_STRIP,offset);
    }

    void draw_triangles_simple(const vector<vec2>& vertices, uint size, vec4 color,vec2 offset)
    {
        draw_unicolor(vertices,size,color,GL_TRIANGLE_FAN,offset);
    }


    void draw_segments(const vector<vec2>& vertices, uint size, vec4 color,vec2 offset)
    {
        draw_unicolor(vertices,size,color,GL_LINES,offset);
    }

    void draw_data(const vector<float>& Data,float _min, float _max)
    {

        data = Data;
        Min = _min;
        Max = _max;
    }


    void draw_circle( const vec2& offset,const vec4& colour)
    {

        {
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);//Good for texture blending (front texture obscures back)

            //Enable what we need
            glUseProgram(circle_ProgramID);


            mat3 localspace_to_worldspace=
                mat3(
                    vec3(1.0      , 0.0      , 0.0),
                    vec3(0.0      , 1.0      , 0.0),
                    vec3(offset.x, offset.y, 1.0)
                );


            glUniform4f(circle_baseColour_ID ,colour.x,colour.y,colour.z,colour.w);
            glUniform1f(circle_r_ID ,0.2);



            mat3 localspace_to_screenspace = worldspace_to_screenspace*localspace_to_worldspace;



            //OpenGL uniform MatrixN sends an N by N matrix, f says it has floating values and v tells us that we want to send a pointer to the values in the CPU memory, rather than setting the values manually
            glUniformMatrix3fv(circle_matrix_ID, 1, GL_FALSE, &localspace_to_screenspace [0][0]);

            //We only have one vertex attribute, these are the values which change from vertex to vertex
            glEnableVertexAttribArray(circle_VertexPosAttribID);
            glBindBuffer(GL_ARRAY_BUFFER, surf_PosBuffer);
            glVertexAttribPointer
            (
                circle_VertexPosAttribID,//Attribute location, you can either locate this in the program (with glGetAttributeLocation), or you can force the shader to use a particular location from inside the shader, I do the former
                2,                   //Number of the below unit per element (this is a 2D vector, so 2)
                GL_FLOAT,            //Unit, this is single precition GL float
                GL_FALSE,             //Normalized? Nope
                0,                    //Stride and offset, I don't use these for anything
                (void*)0
            );


            glDrawArrays(GL_TRIANGLE_FAN,0,4);

            //Disable everything activated
            glDisableVertexAttribArray(circle_VertexPosAttribID);
            glUseProgram(0);
        }
    }



    void draw_hex( const mat3& bodyspace_to_worldspace, vec2 hex_coord, uint tex)
    {

        if (hexagon_textures.size()!=0)
        {
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);//Good for texture blending (front texture obscures back)


            if (tex>=hexagon_textures.size())
                tex =hexagon_textures.size()-1;

            //Enable what we need
            glUseProgram(surf_ProgramID);

            //OpenGL has a number of different texture "slots" here I activate slot number 0
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D,  hexagon_textures[tex]);//Here I bind the texture I want to display, it is bound to the currently active texture slot (0 here)
            //1i tells opengl that this is 1 integer, 2f would be a 2-d vector with floating point data. It is very important that the type matches the data in the glsl shader.
            glUniform1i(surf_colorTex_ID, 0);//Here I set the uniform texture sampler to integer 0, OpenGL understand that this means it shall look in slot 0


            vec2 offset = hexcoord_to_cartesian*hex_coord;

            mat3 localspace_to_bodyspace =
                mat3(
                    vec3(1.0      , 0.0      , 0.0),
                    vec3(0.0      , 1.0      , 0.0),
                    vec3(offset.x, offset.y, 1.0)
                );



        mat3 localspace_to_screenspace = worldspace_to_screenspace*bodyspace_to_worldspace*localspace_to_bodyspace;



            //OpenGL uniform MatrixN sends an N by N matrix, f says it has floating values and v tells us that we want to send a pointer to the values in the CPU memory, rather than setting the values manually
            glUniformMatrix3fv(surf_matrix_ID, 1, GL_FALSE, &localspace_to_screenspace [0][0]);


            mat3 localsapce_to_texturespace=hextex_matrix * localspace_to_bodyspace;


            //This is not an animation, send a unit matrix to "transform" the UV coordinates nowhere
            glUniformMatrix3fv(surf_tex_matrix_ID, 1, GL_FALSE, &(localsapce_to_texturespace[0][0]));//Same procedure as before


            //We only have one vertex attribute, these are the values which change from vertex to vertex
            glEnableVertexAttribArray(surf_VertexPosAttribID);
            glBindBuffer(GL_ARRAY_BUFFER, hex_PosBuffer);
            glVertexAttribPointer
            (
                surf_VertexPosAttribID,//Attribute location, you can either locate this in the program (with glGetAttributeLocation), or you can force the shader to use a particular location from inside the shader, I do the former
                2,                   //Number of the below unit per element (this is a 2D vector, so 2)
                GL_FLOAT,            //Unit, this is single precition GL float
                GL_FALSE,             //Normalized? Nope
                0,                    //Stride and offset, I don't use these for anything
                (void*)0
            );


            glDrawArrays(GL_TRIANGLE_FAN,0,6);

            //Disable everything activated
            glDisableVertexAttribArray(Line_VertexPosAttribID);
            glUseProgram(0);
        }
    }




    void draw_unicolor(const vector<vec2>& vertices, uint size, vec4 color,GLuint displaymode, vec2 offset)
    {
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);//Good for texture blending (front texture obscures back)


        GLuint Buffer=-1;

        size =size > vertices.size() ? vertices.size() : size;
        glGenBuffers(1, &Buffer);
        glBindBuffer( GL_ARRAY_BUFFER, Buffer);
        glBufferData( GL_ARRAY_BUFFER,  sizeof(vec2)*(size), &(vertices[0]), GL_DYNAMIC_DRAW );

        //Enable what we need
        glUseProgram(Line_ProgramID);

        glUniformMatrix3fv(Line_matrix_ID, 1, GL_FALSE, &worldspace_to_screenspace[0][0]);
        glUniform4f(Line_color_ID,color.x,color.y,color.z,color.w);

        glEnableVertexAttribArray(Line_VertexPosAttribID);


        glBindBuffer( GL_ARRAY_BUFFER, Buffer);
        glVertexAttribPointer
        (
            Line_VertexPosAttribID,//Attribute location, you can either locate this in the program, or you can force the shader to use a particular location from inside the shader, I do the former
            2,                   //Number Number of the below unit per element (this is a 2D vector, so 2)
            GL_FLOAT,            //Unit, this is single precition GL float
            GL_FALSE,             //Normalized? Nope
            0,                    //Stride and offset, I don't use these for anything
            (void*)0
        );

        glDrawArrays(displaymode,0,size);

        //Disable everything activated
        glDisableVertexAttribArray(Line_VertexPosAttribID);
        glUseProgram(0);

        //Delete the buffer of this object
        if (Buffer != (GLuint)-1)
            glDeleteBuffers(1,&Buffer);
    }


    vec2 pixel_to_worldspace(uint x, uint y)
    {
        return screenspace_to_worldspace*vec3(2*float(x)/win_w-1.0,-2*float(y)/win_h+1.0,1);
    }

}

