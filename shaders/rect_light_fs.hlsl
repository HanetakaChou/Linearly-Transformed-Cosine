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
	in float4 d3d_Position : SV_POSITION,
	out float4 out_color : SV_TARGET0
	)
{
	const float3 lcol = float3(intensity, intensity, intensity);

	float3 col = lcol;
	out_color = float4(col, 1.0);
}