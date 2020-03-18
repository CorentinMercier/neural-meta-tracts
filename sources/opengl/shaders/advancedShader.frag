#version 450

layout (location = 0) in vec4 evalColor;
layout (location = 1) in vec3 n;
layout (location = 2) in vec3 p;
layout (location = 0) out vec4 color;

uniform vec3 light=vec3(0.5,0.3,5.0);

uniform bool wireframeMode;

void main(void)
{
	vec4 fiberColor = evalColor;
	if (wireframeMode)
		fiberColor = vec4(0.0, 0.0, 0.0, 0.0);
		
	vec3 normal=normalize(n);
	vec3 lightDirection=normalize(light-p);
	vec3 userDirection=normalize(-p);
	vec3 reflexion=reflect(-lightDirection, normal);

	float diffuse_term=0.8*clamp(abs(dot(normal, lightDirection)),0.0,1.0);
    float specular_term=0.5*pow(clamp(dot(reflexion, userDirection),0.0,1.0),64.0);
    float ambiant_term=0.4;

    color = (ambiant_term+diffuse_term)*fiberColor+specular_term*vec4(1.0, 1.0, 1.0, 0.0);//vec4(0.0, 0.5, 0.9, 1.0);  vec4(n, 0.0);//
}
