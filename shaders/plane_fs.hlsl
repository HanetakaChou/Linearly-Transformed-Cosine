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
	float culling_range;
};

SamplerState clamp_point_sampler : register(s0);

float3 ToLinear(float3 v)
{
	return pow(v, 2.2);
}

#include "LTC.hlsli"

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

	const float3 points[4] = {rect_light_vetices[0].xyz, rect_light_vetices[1].xyz, rect_light_vetices[2].xyz, rect_light_vetices[3].xyz};
	const float3 light_color = float3(intensity, intensity, intensity);

	float3 P = in_position;
	float3 N = in_normal;
	float3 V = normalize(eye_position - in_position);
	float3 diffuse_color = ToLinear(dcolor);
	float3 specular_color = ToLinear(scolor);

	float light_attenuation = EvaluateBRDFLTCLightAttenuation(culling_range, P, points);
	if (light_attenuation > 0.0)
	{
		output_color += light_color * light_attenuation * EvaluateBRDFLTC(diffuse_color, roughness, specular_color, P, N, V, points);
	}

	out_color = float4(output_color, 1.0);
}