#version 330

in block {
	vec3 uv; 
} In;

out float out_color;

uniform sampler3D volume_sampler;
uniform float gain = 1.0;

void main(){
	out_color = texture(volume_sampler, In.uv) * gain;
}
