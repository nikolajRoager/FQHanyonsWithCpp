#version 400 core

//Wait aren't transformation matrices 4D? yes, normally, but this is a 2D application, so I can chop of one dimension, alas, this does mean that I won't be using any of glm's build in transform functions, but nevermind, the matrices are pretty simple
uniform mat3 Input_to_DC;//Transformation from uv texture coordinates to device coordinates. The projection matrix is, in this case, just scaling from 0,1 by 0,1 ti -1,1 by -1,1
in vec2 vertex_pos;


out vec2 fragment_pos;

void main()
{
    fragment_pos=vertex_pos;
    gl_Position = vec4((Input_to_DC*vec3(vertex_pos,1)).xy,0,1);
}

