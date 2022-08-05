#version 330

in block {
	vec2 uv; 
} In;

out vec4 out_color;

uniform vec3 color = vec3(1.0);
uniform bool invert = false;
uniform sampler2D fb_sampler;
uniform float gain = 1.0;

void main(){
	float tex_value = texture(fb_sampler, In.uv);
	if (invert){
		tex_value = 1 - tex_value;
	}
	out_color = vec4(tex_value * color * gain, 1.0);
}
