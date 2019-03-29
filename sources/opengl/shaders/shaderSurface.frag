#version 450

in vec3 n;
in vec3 p;

uniform vec3 singleColor;
layout (location = 0) out vec4 color;

uniform vec3 light=vec3(0.5,0.3,5.0);

void main(void)
{
	vec3 normal=normalize(n);
	vec3 lightDirection=normalize(light-p);
	vec3 userDirection=normalize(-p);
	vec3 reflexion=reflect(-lightDirection, normal);

	float diffuse_term=0.8*clamp(abs(dot(normal, lightDirection)),0.0,1.0);
    float specular_term=0.5*pow(clamp(dot(reflexion, userDirection),0.0,1.0),64.0);
    float ambiant_term=0.4;
    vec3 newColor = singleColor;
    color = (ambiant_term+diffuse_term)*vec4(newColor,1.0)+specular_term*vec4(1.0, 1.0, 1.0, 0.0);//vec4((normal+vec3(1.0,1.0,1.0))/2.0,0.0);//
}
