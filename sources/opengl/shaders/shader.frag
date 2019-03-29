#version 450

layout (location = 0) in vec4 vertexToFragmentColor;
layout (location = 0) out vec4 color;

void main(void)
{
    color = vertexToFragmentColor;//vec4(0.0, 0.5, 0.9, 1.0);
}
