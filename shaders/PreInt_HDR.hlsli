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

brx_int2 brx_hdr_application_bridge_get_specular_fresnel_factor_lut_dimension()
{
    uint out_width;
    uint out_height;
    uint out_elements;
    uint out_number_of_levels;
    PreIntegratedHDR_LUT.GetDimensions(0, out_width, out_height, out_elements, out_number_of_levels);
    return brx_int2(out_width, out_height);
}

brx_float2 brx_hdr_application_bridge_get_specular_fresnel_factor_lut(in brx_float2 in_lut_uv)
{
    return PreIntegratedHDR_LUT.SampleLevel(clamp_point_sampler, float3(in_lut_uv, float(PreIntegratedHDR_LUT_TR_INDEX)), 0.0).rg;
}

#endif
