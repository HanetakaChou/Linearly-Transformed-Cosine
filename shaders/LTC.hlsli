#ifndef _LTC_HLSLI_
#define _LTC_HLSLI_ 1

//
// Copyright (C) YuqiaoZhang(HanetakaChou)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//

#include "../thirdparty/Hemispherical-Directional-Reflectance/shaders/brx_hemispherical_directional_reflectance.bsli"
#include "../thirdparty/Linearly-Transformed-Cosine/shaders/brx_linearly_transformed_cosine_attenuation.bsli"
#include "../thirdparty/Linearly-Transformed-Cosine/shaders/brx_linearly_transformed_cosine_radiance.bsli"
#include "PreInt_HDR.hlsli"

#if 0

float3 DiffuseLambertLTC(float3 diffuse_color, float3 vertices_tangent_space[4])
{
	float form_factor_over_quad = EvaluateFormFactorOverQuad(vertices_tangent_space);

	float3 radiance_diffuse = Diffuse_Lambert(diffuse_color) * PI * form_factor_over_quad;
	return radiance_diffuse;
}

float3 DiffuseBurleyLTC(float3 diffuse_color, float roughness, float3 N, float3 V, float3 vertices_tangent_space[4])
{
	// TODO: In Unity3D, another LUT "LTC_DISNEY_DIFFUSE_MATRIX_INDEX" is provided.

	float3 vector_form_factor_over_quad = EvaluateVectorFormFactorOverQuad(vertices_tangent_space);
	float form_factor_over_quad = EvaluateFormFactorOverQuad(vertices_tangent_space);

	// UE4: RectIrradianceLambert
	float3 L = normalize(vector_form_factor_over_quad);
	float non_clamped_NdotL = form_factor_over_quad / length(vector_form_factor_over_quad);
	float non_clamped_NdotV = dot(N, V);

	// Real-Time Rendering Fourth Edition / 9.8 BRDF Models for Surface Reflection: "Hammon[657]. Hammon also shows a method to optimize the BRDF implementation by calculating n·h and l·h without calculating the vector h itself."
	// |L + V|^2 = L^2 + V^2 + 2L·V = 1 + 1 + 2L·V
	// N·H = N·((L + V)/(|L + V|)) = (N·L + N·V)/(|L + V|)
	// L·H = L·((L + V)/(|L + V|)) = (L^2 + L·V)/(|L + V|)
	// V·H = V·((L + V)/(|L + V|)) = (L·V + V^2)/(|L + V|)
	// 2L·H = 2V·H = L·H + V·H = (L^2 + L·V + L·V + V^2)/(|L + V|) = (|L + V|^2)/(|L + V|) = |L + V| = 1 + 1 + 2L·V
	// ⇒ L·H = 0.5 * |L + V| = (0.5 * |L + V|^2)/(|L + V|) =(0.5 * (1 + 1 + 2L·V))/(|L + V|) = 1/(|L + V|) + (L·V)/(|L + V|)
	// UE: [Init](https://github.com/EpicGames/UnrealEngine/blob/4.27/Engine/Shaders/Private/BRDF.ush#L31)
	// U3D: [GetBSDFAngle](https://github.com/Unity-Technologies/Graphics/blob/v10.8.0/com.unity.render-pipelines.core/ShaderLibrary/CommonLighting.hlsl#L361)
	float non_clamped_VdotL = dot(V, L);
	float invLenH = rsqrt(2.0 + 2.0 * non_clamped_VdotL);
	float NdotH = saturate((non_clamped_NdotL + non_clamped_NdotV) * invLenH);
	float LdotH = saturate(invLenH * non_clamped_VdotL + invLenH);

	// UE: [DefaultLitBxDF]https://github.com/EpicGames/UnrealEngine/blob/4.27/Engine/Shaders/Private/ShadingModels.ush#L218
	// U3D: [ClampNdotV](https://github.com/Unity-Technologies/Graphics/blob/v10.8.0/com.unity.render-pipelines.core/ShaderLibrary/CommonLighting.hlsl#L349)
	float NdotV = saturate(abs(non_clamped_NdotV) + 1e-5);
	float NdotL = saturate(non_clamped_NdotL);

	float3 radiance_diffuse = Diffuse_Burley(diffuse_color, roughness, NdotV, NdotL, LdotH) * PI * form_factor_over_quad;
	return radiance_diffuse;
}

float3 SpecularTRLTC(float3 specular_color, float roughness, float3 N, float3 V, float3 vertices_tangent_space[4])
{
	float3x3 linear_transform_inversed = LTC_MATRIX_TR(roughness, saturate(dot(N, V)));

	// [Hill 2016] [Stephen Hill. "LTC Fresnel Approximation." SIGGRAPH 2016.](https://blog.selfshadow.com/publications/s2016-advances/)
	float3 norm = PreIntegratedHDR_TR(specular_color, roughness, saturate(dot(N, V)));

	// LT "linear transform"
	float3 vertices_tangent_space_linear_transformed[4] = {
		mul(linear_transform_inversed, vertices_tangent_space[0]),
		mul(linear_transform_inversed, vertices_tangent_space[1]),
		mul(linear_transform_inversed, vertices_tangent_space[2]),
		mul(linear_transform_inversed, vertices_tangent_space[3])};

	float form_factor_over_quad = EvaluateFormFactorOverQuad(vertices_tangent_space_linear_transformed);

	float3 radiance_specular = norm * form_factor_over_quad;

	return radiance_specular;
}

float3 DualSpecularTRLTC(float material_roughness_0, float material_roughness_1, float material_lobe_mix, float3 specular_color, float roughness, float subsurface_mask, float3 N, float3 V, float3 vertices_tangent_space[4])
{
	float material_roughness_average = lerp(material_roughness_0, material_roughness_1, material_lobe_mix);
	float average_to_roughness_0 = material_roughness_0 / material_roughness_average;
	float average_to_roughness_1 = material_roughness_1 / material_roughness_average;

	float surface_roughness_average = roughness;
	float surface_roughness_0 = max(saturate(average_to_roughness_0 * surface_roughness_average), 0.02);
	float surface_roughness_1 = saturate(average_to_roughness_1 * surface_roughness_average);

	// UE4: SubsurfaceProfileBxDF
	surface_roughness_0 = lerp(1.0f, surface_roughness_0, saturate(10.0 * subsurface_mask));
	surface_roughness_1 = lerp(1.0f, surface_roughness_1, saturate(10.0 * subsurface_mask));

	float3 radiance_specular_0 = SpecularTRLTC(specular_color, surface_roughness_0, N, V, vertices_tangent_space);
	float3 radiance_specular_1 = SpecularTRLTC(specular_color, surface_roughness_1, N, V, vertices_tangent_space);
	float3 radiance_specular = lerp(radiance_specular_0, radiance_specular_1, material_lobe_mix);
	return radiance_specular;
}
#endif

Texture2DArray LTC_MATRIX_LUT : register(t0);

#define LTC_MATRIX_LUT_TR_INDEX 0

brx_int2 brx_ltc_application_bridge_get_specular_matrix_lut_dimension()
{
	uint out_width;
	uint out_height;
	uint out_elements;
	uint out_number_of_levels;
	LTC_MATRIX_LUT.GetDimensions(0, out_width, out_height, out_elements, out_number_of_levels);

	return brx_int2(out_width, out_height);
}

brx_float4 brx_ltc_application_bridge_get_specular_matrix_lut(in brx_float2 in_lut_uv)
{
	float4 ltc_matrix_lut_tr = LTC_MATRIX_LUT.SampleLevel(clamp_point_sampler, float3(in_lut_uv, float(LTC_MATRIX_LUT_TR_INDEX)), 0.0);
	return ltc_matrix_lut_tr;
}

#endif
