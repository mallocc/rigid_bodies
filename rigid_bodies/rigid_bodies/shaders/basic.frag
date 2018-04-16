#version 400 core



// ins
in vec3 o_color;



void main() 
{
// apply fragment color
	gl_FragColor = vec4(o_color,1);
}