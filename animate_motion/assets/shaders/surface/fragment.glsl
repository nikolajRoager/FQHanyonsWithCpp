#version 400 core

uniform sampler2D   colorSampler;

in vec2 fragment_uv;
in vec2 fragment_pos_uv;
out vec4 color;


uniform int w;
uniform int h;


uniform float Min;
uniform float Max;

uniform float Q[256];

void main()
{

    //color = texture(colorSampler,fragment_uv);

    int lat_x = int(fragment_uv.x*w);
    int lat_y = int(fragment_uv.y*h);

    float val =(Q[lat_x+(h-1-lat_y)*w]-Min)/(Max-Min);

    float fac = fract(val*4);
    switch (int (val*4))
    {
    case 0:
        color= vec4(vec3(0,1-fac,1),1.0);
    break;
    case 1:
        color= vec4(vec3(0,0,1-fac),1.0);
    break;
    case 2:
        color= vec4(vec3(fac,0,0),1.0);
    break;
    case 3:
        color= vec4(vec3(1,fac,0),1.0);
    break;
    default:
        if (val>1.0)
            color= vec4(vec3(1,1,0),1);
        else
            color= vec4(vec3(0,1,1),1);
    }

}

