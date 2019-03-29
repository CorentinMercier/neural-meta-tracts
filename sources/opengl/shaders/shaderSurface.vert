#version 450

uniform mat4 camera_modelview;
uniform mat4 camera_projection;

layout (location = 0) in vec4 position;
layout (location = 1) in vec3 normal;

out vec3 n;
out vec3 p;

void main(void)
{
	n = normal;
	p= vec3(camera_modelview * position);
    gl_Position = camera_projection * camera_modelview * position;
}
