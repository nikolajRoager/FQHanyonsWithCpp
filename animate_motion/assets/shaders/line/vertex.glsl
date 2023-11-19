#version 400 core

//Since this is 2D we can get away with 3D transformation matrices, much the same way 3D require 4D ones
in vec2 vertex_worldspace;

uniform mat3 worldspace_to_DC;

void main()
{
    vec2 fragment_pos =(worldspace_to_DC*vec3(vertex_worldspace,1)).xy;
    gl_Position = vec4(fragment_pos,0,1);//Should draw a square on the screen, to verrify all libraries and buffers are set up correctly
}

