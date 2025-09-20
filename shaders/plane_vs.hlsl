cbuffer _unused_name_uniform_buffer_global_layout_per_frame_binding : register(b0)
{
	// mesh
	column_major float4x4 model_transform;
	float3 dcolor;
	float _padding_dcolor;
	float3 scolor;
	float _padding_scolor;
	float roughness;
	float _padding_roughness_1;
	float _padding_roughness_2;
	float _padding_roughness_3;

	// camera
	column_major float4x4 view_transform;
	column_major float4x4 projection_transform;
	float3 eye_position;
	float _padding_eye_position;

	// light
	float4 rect_light_vetices[4];
	float intensity;
	float twoSided;
};

void main(
	in float3 in_position : POSITION0, 
	in float3 in_normal : NORMAL0, 
	out float4 d3d_Position : SV_POSITION, 
	out float3 out_position : TEXCOORD0, 
	out float3 out_normal : TEXCOORD1
	)
{
	float3 model_position = in_position;
	float3 world_position = mul(model_transform, float4(model_position, 1.0)).xyz;
	float4 clip_position = mul(projection_transform, mul(view_transform, float4(world_position, 1.0)));

	float3 model_normal = in_normal;
	// TODO: normal transform
	float3x3 tangent_transform = float3x3(model_transform[0].xyz, model_transform[1].xyz, model_transform[2].xyz);
	float3 world_normal = normalize(mul(tangent_transform, model_normal));

	out_position = world_position;
	out_normal = world_normal;
	d3d_Position = clip_position;
}