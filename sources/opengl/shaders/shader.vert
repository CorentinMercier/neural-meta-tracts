#version 450

uniform mat4 camera_modelview;
uniform mat4 camera_projection;

layout (location = 0) in vec4 position;
layout (location = 1) in vec4 color;

layout (location = 0) out vec4 vertexToFragmentColor;

void main(void)
{
    vertexToFragmentColor = color;
    gl_Position = camera_projection * camera_modelview * position;
}
