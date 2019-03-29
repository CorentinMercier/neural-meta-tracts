#version 450

uniform mat4 camera_modelview;
uniform mat4 camera_projection;

layout (location = 0) in vec4 position;

void main(void)
{
    gl_Position = camera_projection * camera_modelview * position;
}
