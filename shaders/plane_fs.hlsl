// https://github.com/selfshadow/ltc_code/tree/master/webgl/shaders/ltc/ltc_quad.fs

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
	float _unused_culling_range;
};

SamplerState clamp_point_sampler : register(s0);

float3 ToLinear(float3 v)
{
	return pow(v, 2.2);
}

#include "LTC.hlsli"

#define INTERNAL_IRRADIANCE_THRESHOLD 1e-4

#define INTERNAL_ATTENUATION_THRESHOLD 1e-4

void main(
	in float4 d3d_Position
	: SV_POSITION,
	  in float3 in_position
	: TEXCOORD0,
	  in float3 in_normal
	: TEXCOORD1,
	  out float4 out_color
	: SV_TARGET0)
{
	float3 output_color = float3(0.0, 0.0, 0.0);

	const float3 quad_vertices_world_space[4] = {rect_light_vetices[0].xyz, rect_light_vetices[1].xyz, rect_light_vetices[2].xyz, rect_light_vetices[3].xyz};
	const float3 area_lighting_radiance = float3(intensity, intensity, intensity);

	float3 P = in_position;
	float3 N = in_normal;
	float3 V = normalize(eye_position - in_position);
	float3 diffuse_color = ToLinear(dcolor);
	float3 specular_color = ToLinear(scolor);

	brx_float quad_culling_range;
	{
		brx_float quad_area = brx_length(brx_cross(quad_vertices_world_space[1] - quad_vertices_world_space[0], quad_vertices_world_space[3] - quad_vertices_world_space[0]));   
		brx_float maximum_area_lighting_radiance = brx_max(brx_max(area_lighting_radiance.x, area_lighting_radiance.y), area_lighting_radiance.z);

		quad_culling_range = brx_sqrt((maximum_area_lighting_radiance * quad_area) / INTERNAL_IRRADIANCE_THRESHOLD);
	}

	float area_lighting_attenuation = brx_ltc_attenuation(quad_culling_range, P, quad_vertices_world_space);
	if (area_lighting_attenuation > INTERNAL_ATTENUATION_THRESHOLD)
	{
		brx_float3 diffuse_albedo = diffuse_color;
		
		brx_float3 specular_albedo = brx_hdr_specular_albedo(specular_color, roughness, N, V);
		
		brx_float3 radiance = brx_ltc_radiance(diffuse_albedo, specular_albedo, roughness, N, V, area_lighting_radiance, P, quad_vertices_world_space);

		output_color += area_lighting_attenuation * radiance;
	}

	out_color = float4(output_color, 1.0);
}