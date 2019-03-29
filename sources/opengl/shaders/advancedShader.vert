#version 450

layout (location = 0) in vec4 position;
layout (location = 1) in vec4 color;
layout (location = 2) in vec2 ab;
layout (location = 3) in vec3 ellipseOrientation;

layout (location = 0) out vec4 vertColor;
layout (location = 1) out vec2 vertAb;
layout (location = 2) out vec3 vertEllipseOrientation;

void main(void)
{
    vertColor = color;
    vertAb = ab;
    vertEllipseOrientation = ellipseOrientation;

    gl_Position = position;
}
