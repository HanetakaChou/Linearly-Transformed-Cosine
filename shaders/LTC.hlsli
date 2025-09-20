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

#include "BRDF.hlsli"
#include "PreInt_HDR.hlsli"

// This function is provided by the user
float3x3 LTC_MATRIX_TR(float roughness, float NoV);

// [in] culling_range: a finite range to support light culling 
// [in] P: The surface position in world space.
// [in] vertices_world_space: The vertices of the quad in world space. The facing of the quad is determined by the winding order of the vertices.
float EvaluateBRDFLTCLightAttenuation(float culling_range, float3 P, float3 vertices_world_space[4]);

// [in] P: The surface position in world space.
// [in] N: The surface normal in world space.
// [in] V: The outgoing direction in world space.
// [in] vertices_world_space: The vertices of the quad in world space. The facing of the quad is determined by the winding order of the vertices.
float3 EvaluateBRDFLTC(float3 diffuse_color, float roughness, float3 specular_color, float3 P, float3 N, float3 V, float3 vertices_world_space[4]);

// [in] vertices_tangent_space: The vertices of the quad in tangent space. The facing of the quad is determined by the winding order of the vertices.
float3 DiffuseLambertLTC(float3 diffuse_color, float3 vertices_tangent_space[4]);

// [in] vertices_tangent_space: The vertices of the quad in tangent space. The facing of the quad is determined by the winding order of the vertices.
float3 DiffuseBurleyLTC(float3 diffuse_color, float roughness, float3 N, float3 V, float3 vertices_tangent_space[4]);

// [in] vertices_tangent_space: The vertices of the quad in tangent space. The facing of the quad is determined by the winding order of the vertices.
float3 SpecularTRLTC(float3 specular_color, float roughness, float3 N, float3 V, float3 vertices_tangent_space[4]);

// [in] vertices_tangent_space: The vertices of the quad in tangent space. The facing of the quad is determined by the winding order of the vertices.
float3 DualSpecularTRLTC(float material_roughness_0, float material_roughness_1, float material_lobe_mix, float3 specular_color, float roughness, float subsurface_mask, float3 N, float3 V, float3 vertices_tangent_space[4]);

// [in] vertices_tangent_space: The vertices of the quad in tangent space. The facing of the quad is determined by the winding order of the vertices.
float EvaluateFormFactorOverQuad(float3 vertices_tangent_space[4]);

// [in] vertices_tangent_space: The vertices of the quad in tangent space. The facing of the quad is determined by the winding order of the vertices.
float3 EvaluateVectorFormFactorOverQuad(float3 vertices_tangent_space[4]);

// [in] v1: The first normalized vertex which is projected onto the sphere. The facing of the quad is determined by the winding order of the vertices.
// [in] v2: The second normalized vertex which is projected onto the sphere. The facing of the quad is determined by the winding order of the vertices.
float3 EvaluateVectorFormFactorOverQuadEdge(float3 v1, float3 v2);

// [in] cos_elevation_angle : implies the direction of the vector irrandiance of the sphere proxy
// [in] sin_angular_extent : implies the length of the vector irrandiance of the sphere proxy
float EvaluateFormFactorOverSphere(float cos_elevation_angle, float sin_angular_extent);

float EvaluateBRDFLTCLightAttenuation(float culling_range, float3 P, float3 vertices_world_space[4])
{
	float3 v0_to_v1 = vertices_world_space[1] - vertices_world_space[0];
	float3 v0_to_v3 = vertices_world_space[3] - vertices_world_space[0];
	float3 normal_back_face = normalize(cross(v0_to_v1, v0_to_v3));
	float3 L = vertices_world_space[0].xyz - P;
	float z = dot(L, normal_back_face);

	[branch]
	if (z > 0.0)
	{
		float x;
		float y;
		{
			// [FetchDiffuseFilteredTexture](http://blog.selfshadow.com/publications/ltc/ltc_demo.zip)
			float3 v0_to_p_projected = normal_back_face * z - L;

			// v0_to_p_projected = u * v0_to_v1 + v * v0_to_v3
			// ⇒
			// v0_to_p_projected_dot_v0_to_v1 =  u * v0_to_v1_2 + v * v0_to_v3_dot_v0_to_v1
			// v0_to_p_projected_dot_v0_to_v3 =  u * v0_to_v1_dot_v0_to_v3 + v * v0_to_v3_2
			// ⇒
			// u = (v0_to_p_projected_dot_v0_to_v1 * v0_to_v3_2 - v0_to_p_projected_dot_v0_to_v3 * v0_to_v3_dot_v0_to_v1) / (v0_to_v1_2 * v0_to_v3_2 - v0_to_v1_dot_v0_to_v3 * v0_to_v3_dot_v0_to_v1)
			// v = (v0_to_p_projected_dot_v0_to_v3 * v0_to_v1_2 - v0_to_p_projected_dot_v0_to_v1 * v0_to_v1_dot_v0_to_v3) / (v0_to_v1_2 * v0_to_v3_2 - v0_to_v1_dot_v0_to_v3 * v0_to_v3_dot_v0_to_v1)
			float v0_to_p_projected_dot_v0_to_v1 = dot(v0_to_p_projected, v0_to_v1);
			float v0_to_p_projected_dot_v0_to_v3 = dot(v0_to_p_projected, v0_to_v3);
			float v0_to_v1_2 = dot(v0_to_v1, v0_to_v1);
			float v0_to_v3_2 = dot(v0_to_v3, v0_to_v3);
			float v0_to_v1_dot_v0_to_v3 = dot(v0_to_v1, v0_to_v3);
			float v0_to_v3_dot_v0_to_v1 = v0_to_v1_dot_v0_to_v3;
			float determinant = v0_to_v1_2 * v0_to_v3_2 - v0_to_v1_dot_v0_to_v3 * v0_to_v3_dot_v0_to_v1;
			float u = (v0_to_p_projected_dot_v0_to_v1 * v0_to_v3_2 - v0_to_p_projected_dot_v0_to_v3 * v0_to_v3_dot_v0_to_v1) / determinant;
			float v = (v0_to_p_projected_dot_v0_to_v3 * v0_to_v1_2 - v0_to_p_projected_dot_v0_to_v1 * v0_to_v1_dot_v0_to_v3) / determinant;
			x = u * length(v0_to_v1);
			y = v * length(v0_to_v3);
		}

		float factor_2;
		{
			// [EllipsoidalDistanceAttenuation](https://github.com/Unity-Technologies/Graphics/blob/v10.8.0/com.unity.render-pipelines.high-definition/Runtime/Material/Lit/Lit.hlsl#L1591)
			// [BoxDistanceAttenuation](https://github.com/Unity-Technologies/Graphics/blob/v10.8.0/com.unity.render-pipelines.high-definition/Runtime/Material/Lit/Lit.hlsl#L1597)
			float normalized_x = x / (culling_range + 0.5 * length(v0_to_v1));
			float normalized_y = y / (culling_range + 0.5 * length(v0_to_v3));
			float normalized_z = z / culling_range;
			factor_2 = normalized_x * normalized_x + normalized_y * normalized_y + normalized_z * normalized_z;
		}

		float windowing_function;
		{
			// [Sebastian Lagarde, Charles Rousiers. "Moving Frostbite to PBR." SIGGRAPH 2014.](https://www.ea.com/frostbite/news/moving-frostbite-to-pb)
			// factor = "x (distance) / light_radius"
			// float factor_2 = factor * factor;
			float factor_4 = factor_2 * factor_2;
			float distance_factor = saturate(1.0 - factor_4);
			float smooth_distance_factor = distance_factor * distance_factor;
			windowing_function = smooth_distance_factor;
		}

		return windowing_function;
	}
	else
	{
		// TODO: some parts of the quad may still in the upper hemisphere
		return 0.0;
	}
}

float3 EvaluateBRDFLTC(float3 diffuse_color, float roughness, float3 specular_color, float3 P, float3 N, float3 V, float3 vertices_world_space[4])
{
	float3 radiance = float3(0.0, 0.0, 0.0);

	// Transform the vertices to the tangent space of the current shading position.
	float3 vertices_tangent_space[4];
	{
		// The LUTs are precomputed by assuming that the outgoing direction V is in the XOZ plane, since the TR BRDF is isotropic.
		float3 T1 = normalize(V - N * dot(V, N));

		float3 T2 = cross(N, T1);

		float4x4 world_to_tangent_transform = float4x4(
			float4(T1, dot(T1, -P)),   // row 0
			float4(T2, dot(T2, -P)),   // row 1
			float4(N, dot(N, -P)),	   // row 2
			float4(0.0, 0.0, 0.0, 1.0) // row 3
			);

		vertices_tangent_space[0] = mul(world_to_tangent_transform, float4(vertices_world_space[0], 1.0)).xyz;
		vertices_tangent_space[1] = mul(world_to_tangent_transform, float4(vertices_world_space[1], 1.0)).xyz;
		vertices_tangent_space[2] = mul(world_to_tangent_transform, float4(vertices_world_space[2], 1.0)).xyz;
		vertices_tangent_space[3] = mul(world_to_tangent_transform, float4(vertices_world_space[3], 1.0)).xyz;
	}

	// radiance += DiffuseLambertLTC(diffuse_color, vertices_tangent_space);
	radiance += DiffuseBurleyLTC(diffuse_color, roughness, N, V, vertices_tangent_space);

	// radiance += SpecularTRLTC(specular_color, roughness, N, V, vertices_tangent_space);
	radiance += DualSpecularTRLTC(0.75, 1.30, 0.85, specular_color, roughness, 1.0, N, V, vertices_tangent_space);

	return radiance;
}

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
		mul(linear_transform_inversed, vertices_tangent_space[3]) };

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

float EvaluateFormFactorOverQuad(float3 vertices_tangent_space[4])
{
	// [Hill 2016] [Stephen Hill, Eric Heitz. "Real-Time Area Lighting: a Journey from Research to Production." SIGGRAPH 2016.](https://blog.selfshadow.com/publications/s2016-advances/)
	// Theory & Implementation / 3. Clip Polygon to upper hemisphere

	// The vector form factor can be calculated even if the quad id NOT horizon-clipped
	float3 vector_form_factor_over_quad = EvaluateVectorFormFactorOverQuad(vertices_tangent_space);

	// Introduce the sphere proxy with the same vector form factor
	float cos_elevation_angle = normalize(vector_form_factor_over_quad).z;
	float sin_angular_extent = sqrt(length(vector_form_factor_over_quad));

	float form_factor_over_sphere = EvaluateFormFactorOverSphere(cos_elevation_angle, sin_angular_extent);

	return form_factor_over_sphere;
}

float3 EvaluateVectorFormFactorOverQuad(float3 vertices_tangent_space[4])
{
	// [Heitz 2017] [Eric Heitz. "Geometric Derivation of the Irradiance of Polygonal Lights." Technical report 2017.](https://hal.archives-ouvertes.fr/hal-01458129)

	float3 vertices_normalized[4] = {
		normalize(vertices_tangent_space[0]),
		normalize(vertices_tangent_space[1]),
		normalize(vertices_tangent_space[2]),
		normalize(vertices_tangent_space[3]) };

	float3 vector_form_factor_over_quad = float3(0.0, 0.0, 0.0);
	vector_form_factor_over_quad += EvaluateVectorFormFactorOverQuadEdge(vertices_normalized[0], vertices_normalized[1]);
	vector_form_factor_over_quad += EvaluateVectorFormFactorOverQuadEdge(vertices_normalized[1], vertices_normalized[2]);
	vector_form_factor_over_quad += EvaluateVectorFormFactorOverQuadEdge(vertices_normalized[2], vertices_normalized[3]);
	vector_form_factor_over_quad += EvaluateVectorFormFactorOverQuadEdge(vertices_normalized[3], vertices_normalized[0]);

	return vector_form_factor_over_quad;
}

float3 EvaluateVectorFormFactorOverQuadEdge(float3 v1, float3 v2)
{
	// [Hill 2016] [Stephen Hill, Eric Heitz. "Real-Time Area Lighting: a Journey from Research to Production." SIGGRAPH 2016.](https://blog.selfshadow.com/publications/s2016-advances/)
	// Theory & Implementation / 4. Compute edge intergrals

	// acos(dot(v1, v2)) * normalized(cross(v1, v2)) * (1 / 2PI)
	// = acos(dot(v1, v2)) * cross(v1, v2)
	// = (acos(dot(v1, v2)) * (1 / sin(acos(dot(v1, v2)))) * (1 / 2PI)) * cross(v1, v2)

	// cubic rational fit
	// theta_sintheta ≈ (acos(dot(v1, v2)) * (1 / sin(acos(dot(v1, v2)))) * (1 / 2PI))

	float x = dot(v1, v2);
	float y = abs(x);

	// (1 / 2PI) has been multiplied here
	float a = 0.8543985 + (0.4965155 + 0.0145206 * y) * y;
	float b = 3.4175940 + (4.1616724 + y) * y;
	float v = a / b;

	float theta_sintheta = (x > 0.0) ? v : 0.5 * rsqrt(max(1.0 - x * x, 1e-7)) - v;

	return cross(v1, v2) * theta_sintheta;
}

float EvaluateFormFactorOverSphere(float cos_omega, float sin_sigma)
{
	// [Snyder 1996]. [John Snyder. "Area Light Sources for Real-Time Graphics." Technical Report 1996.](https://www.microsoft.com/en-us/research/publication/area-light-sources-for-real-time-graphics/)

	float form_factor_over_sphere;
#if 0
	// UE4: SphereHorizonCosWrap  
	// ω ∈ [0, π/2 - σ]
	if(cos_omega > sin_sigma)
	{
		form_factor_over_sphere = cos_omega * sin_sigma * sin_sigma;
	}
	// ω ∈ [π/2 - σ, π]
	else
	{
		float tmp = (sin_sigma + max(cos_omega, -sin_sigma));

		form_factor_over_sphere = tmp * tmp / (4 * sin_sigma) * sin_sigma * sin_sigma;
	}
#else
	form_factor_over_sphere = sin_sigma * sin_sigma * (sin_sigma * sin_sigma + cos_omega) / (sin_sigma * sin_sigma + 1.0);
#endif

	return form_factor_over_sphere;
}

Texture2DArray LTC_MATRIX_LUT : register(t0);

#define LTC_MATRIX_LUT_TR_INDEX 0

float3x3 LTC_MATRIX_TR(float roughness, float NoV)
{
	// float out_width;
	// float out_height;
	// float out_elements;
	// float out_number_of_levels;
	// LTC_MATRIX_LUT.GetDimensions(0, out_width, out_height, out_elements, out_number_of_levels);

	// float2 lut_texture_size = float2(out_width, out_height);
	// float2 lut_bias = float2(0.5, 0.5) / float2(lut_texture_size.x, lut_texture_size.y);
	// float2 lut_scale = float2(1.0, 1.0) - float2(1.0, 1.0) / float2(lut_texture_size.x, lut_texture_size.y);
	// float2 lut_uv = lut_bias + lut_scale * float2(roughness, sqrt(1.0 - NoV));

	const float alpha_min = 0.02;
	float alpha = max(alpha_min, roughness * roughness);

	float2 lut_uv = float2(sqrt(max(0.0, 1.0 - NoV)), log(alpha) / log(alpha_min));

	float4 ltc_matrix_lut_tr = LTC_MATRIX_LUT.SampleLevel(clamp_point_sampler, float3(lut_uv, float(LTC_MATRIX_LUT_TR_INDEX)), 0.0).rgba;

	float3x3 linear_transform_inversed = float3x3(
		float3(ltc_matrix_lut_tr.x, 0.0, ltc_matrix_lut_tr.z), // row 0
		float3(0.0, 1.0, 0.0),								   // row 1
		float3(ltc_matrix_lut_tr.y, 0.0, ltc_matrix_lut_tr.w)  // row 2
	);

	return linear_transform_inversed;
}


#endif