#version 330

in vec2 position;

out block {
	vec2 uv;
} Out;

void main(){
	Out.uv = (position + 1) / 2;
    gl_Position = vec4(position, 0.0, 1.0);
}
