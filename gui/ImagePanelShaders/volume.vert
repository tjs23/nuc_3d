#version 330

in vec2 position;

uniform mat4 model_world_xform = mat4(1.0);
uniform mat4 world_camera_xform = mat4(1.0);
uniform mat4 camera_clip_xform = mat4(1.0);

uniform int slices;
uniform int axis = 0;

out block {
	vec3 uv;
} Out;

void main(){
	vec3 pos = vec3(0);
	switch (axis){
		case 0:
			pos.yz = position;
			break;
		case 1:
			pos.xz = position;
			break;
		case 2:
			pos.xy = position;
			break;
	}
	pos[axis] = 2 * float(gl_InstanceID)/slices - 1;
	Out.uv = (pos / 2) + 0.5;

  gl_Position = ( camera_clip_xform
                * world_camera_xform
                * model_world_xform
                * vec4(pos, 1.0));
}
