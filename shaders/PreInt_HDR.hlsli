#ifndef _PREINT_HDR_HLSLI_
#define _PREINT_HDR_HLSLI_ 1

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

// TODO: term "HDR" is ambiguous
// HDR (Hemispherical Directional Reflectance)
// HDR (High Dynamic Range)

Texture2DArray PreIntegratedHDR_LUT : register(t1);

#define PreIntegratedHDR_LUT_TR_INDEX 0

float3 PreIntegratedHDR_TR(float3 f0, float f90, float roughness, float NdotV)
{
    // UE: [EnvBRDF](https://github.com/EpicGames/UnrealEngine/blob/4.27/Engine/Shaders/Private/BRDF.ush#L471)
    // U3D: [GetPreIntegratedFGDGGXAndDisneyDiffuse](https://github.com/Unity-Technologies/Graphics/blob/v10.8.0/com.unity.render-pipelines.high-definition/Runtime/Material/PreIntegratedFGD/PreIntegratedFGD.hlsl#L8)

#if 0
    // Remap: [0, 1] -> [0.5/size, 1.0 - 0.5/size]
    // U3D: [Remap01ToHalfTexelCoord](https://github.com/Unity-Technologies/Graphics/blob/v10.8.0/com.unity.render-pipelines.core/ShaderLibrary/Common.hlsl#L661)
    // UE: [N/A](https://github.com/EpicGames/UnrealEngine/blob/4.27/Engine/Shaders/Private/RectLight.ush#L450)
    float out_width;
	float out_height;
    float out_elements;
	float out_number_of_levels;
	PreIntegratedHDR_LUT.GetDimensions(0, out_width, out_height, out_elements, out_number_of_levels);
    float2 lut_texture_size = float2(out_width, out_height);
    float2 lut_bias = float2(0.5, 0.5) / float2(lut_texture_size.x, lut_texture_size.y);
    float2 lut_scale = float2(1.0, 1.0) - float2(1.0, 1.0) / float2(lut_texture_size.x, lut_texture_size.y);
    float2 lut_uv = lut_bias + lut_scale * float2(NdotV, roughness);
#endif

    const float alpha_min = 0.02;
    float alpha = max(alpha_min, roughness * roughness);

    float2 lut_uv = float2(sqrt(max(0.0, 1.0 - NdotV)), log(alpha) / log(alpha_min));

    float2 fresnel_factor = PreIntegratedHDR_LUT.SampleLevel(clamp_point_sampler, float3(lut_uv, float(PreIntegratedHDR_LUT_TR_INDEX)), 0.0).rg;
    float f0_factor = fresnel_factor.x;
    float f90_factor = fresnel_factor.y;

    float3 hdr = f0 * f0_factor + float3(f90, f90, f90) * f90_factor;
    return hdr;
}

float3 PreIntegratedHDR_TR(float3 f0, float roughness, float NdotV)
{
    // UE: [EnvBRDF](https://github.com/EpicGames/UnrealEngine/blob/4.27/Engine/Shaders/Private/BRDF.ush#L476)
    // Anything less than 2% is physically impossible and is instead considered to be shadowing
    float f90 = saturate(50.0 * f0.g);

    return PreIntegratedHDR_TR(f0, f90, roughness, NdotV);
}

#endif
