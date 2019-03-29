#version 450

//in vec4 vertexToFragmentColor;
layout (location = 0) in vec4 evalColor;
layout (location = 0) out vec4 color;

in vec3 p;
in vec3 n;

void main(void)
{
    color = vec4(0.0, 0.0, 0.0, 0.0);///*vertexToFragmentColor;*/ evalColor;//vec4(0.0, 0.5, 0.9, 1.0);
}
