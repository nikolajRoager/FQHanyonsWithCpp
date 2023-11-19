#version 400 core

uniform vec4 color;
void main()
{
    //It really doesn't get more basic than this.
    gl_FragColor =  color;
}

