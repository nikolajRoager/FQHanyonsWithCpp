#version 400 core

uniform vec4 base_colour;
uniform float radius;

in vec2 fragment_pos;
out vec4 colour;


void main()
{

    float r = sqrt(fragment_pos.x*fragment_pos.x+fragment_pos.y*fragment_pos.y);


    colour = base_colour;
    if (r>radius)
        discard;

    colour = base_colour;

}

