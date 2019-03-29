#version 450

uniform vec3 singleColor;
layout (location = 0) out vec4 color;

void main(void)
{
    color = vec4(singleColor, 1.0);
}
