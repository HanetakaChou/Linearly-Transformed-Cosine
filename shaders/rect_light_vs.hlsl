cbuffer _unused_name_uniform_buffer_global_layout_per_frame_binding : register(b0)
{
	// camera
	column_major float4x4 view_transform;
	column_major float4x4 projection_transform;

	// light
	float4 rect_light_vetices[4];
	float intensity;
};

void main(
	in uint d3d_VertexID  : SV_VertexID, 
	out float4 d3d_Position : SV_POSITION
	)
{
	float3 world_position = rect_light_vetices[d3d_VertexID].xyz;
	float4 clip_position = mul(projection_transform, mul(view_transform, float4(world_position, 1.0)));

	d3d_Position = clip_position;
}