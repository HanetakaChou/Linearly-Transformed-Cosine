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

#ifndef _BRX_LINEARLY_TRANSFORMED_COSINE_OPTIMIZATION_H_
#define _BRX_LINEARLY_TRANSFORMED_COSINE_OPTIMIZATION_H_ 1

#include <DirectXMath.h>
#include <cmath>
#include <cfloat>
#include <algorithm>
#include "../../Brioche-Shader-Language/include/brx_low_discrepancy_sequence.h"
#include "../../Brioche-Shader-Language/include/brx_brdf.h"
#include "../../Nelder-Mead/include/brx_fminsearch.h"
#include "../../McRT-Malloc/include/mcrt_vector.h"
#include "../../McRT-Malloc/include/mcrt_parallel_map.h"
#include "../../McRT-Malloc/include/mcrt_parallel_reduce.h"

// UE4：[128](https://github.com/EpicGames/UnrealEngine/blob/4.27/Engine/Source/Runtime/Renderer/Private/SystemTextures.cpp#L322)
// U3D: [4096](https://github.com/Unity-Technologies/Graphics/blob/v10.8.1/com.unity.render-pipelines.core/ShaderLibrary/ImageBasedLighting.hlsl#L340)
static constexpr uint32_t const INTERNAL_BRX_LTC_OPTIMIZATION_MONTE_CARLO_SAMPLE_COUNT = 4096U;

static constexpr uint32_t const INTERNAL_BRX_LTC_OPTIMIZATION_MONTE_CARLO_GRAIN_SIZE = 64U;

static constexpr float const INTERNAL_BRX_LTC_OPTIMIZATION_LENGTH_SQUARE_MINIMUM = 1E-5F;

static constexpr float const INTERNAL_BRX_LTC_OPTIMIZATION_SCALE_MINIMUM = 1E-5F;

static inline void internal_brx_linearly_transformed_cosine_optimization_compute_matrix(uint32_t const lut_width_index, uint32_t const lut_height_index, uint32_t const lut_width, uint32_t const lut_height, float &inout_M2_m11, float &inout_M2_m22, float &inout_M2_m13, DirectX::XMFLOAT4 &out_packed_M);

static inline void brx_linearly_transformed_cosine_optimization_compute_matrices(DirectX::XMFLOAT4 *const out_lut, uint32_t const lut_width, uint32_t const lut_height)
{
    // UE4: [128x32](https://github.com/EpicGames/UnrealEngine/blob/4.27/Engine/Source/Runtime/Renderer/Private/SystemTextures.cpp#L289)
    // U3D: [64x64] (https://github.com/Unity-Technologies/Graphics/blob/v10.8.1/com.unity.render-pipelines.high-definition/Runtime/Material/PreIntegratedFGD/PreIntegratedFGD.cs.hlsl#L10)

    mcrt_vector<float> ndotv_one_M2_m11(static_cast<size_t>(lut_height));
    mcrt_vector<float> ndotv_one_M2_m22(static_cast<size_t>(lut_height));

    for (uint32_t lut_height_index = 0U; lut_height_index < lut_height; ++lut_height_index)
    {
        float M2_m11;
        float M2_m22;
        // roughness -> 1
        // init the hemisphere in which the distribution is fitted
        if (0U == lut_height_index)
        {
            M2_m11 = 1.0F;
            M2_m22 = 1.0F;
        }
        // init with roughness of previous fit (of the same ndotv)
        else
        {
            assert(lut_height_index >= 1U);
            M2_m11 = ndotv_one_M2_m11[lut_height_index - 1U];
            M2_m22 = ndotv_one_M2_m22[lut_height_index - 1U];
        }

        float M2_m13 = 0.0F;

        DirectX::XMFLOAT4 packed_M;

        internal_brx_linearly_transformed_cosine_optimization_compute_matrix(0U, lut_height_index, lut_width, lut_height, M2_m11, M2_m22, M2_m13, packed_M);

        ndotv_one_M2_m11[lut_height_index] = M2_m11;
        ndotv_one_M2_m22[lut_height_index] = M2_m22;
        // isotropic
        assert(std::abs(M2_m13) < 2E-3F);

        out_lut[lut_width * lut_height_index] = packed_M;
    }

    struct parallel_map_user_data
    {
        DirectX::XMFLOAT4 *const out_lut;
        uint32_t const lut_width;
        uint32_t const lut_height;
        float const *ndotv_one_M2_m11;
        float const *ndotv_one_M2_m22;
    };

    parallel_map_user_data user_data = {out_lut, lut_width, lut_height, ndotv_one_M2_m11.data(), ndotv_one_M2_m22.data()};

    mcrt_parallel_map(
        0U,
        lut_height * (lut_width - 1U),
        (lut_width - 1U),
        [](uint32_t begin, uint32_t end, void *wrapped_user_data) -> void
        {
            parallel_map_user_data *unwrapped_user_data = static_cast<parallel_map_user_data *>(wrapped_user_data);
            DirectX::XMFLOAT4 *const out_lut = unwrapped_user_data->out_lut;
            uint32_t const lut_width = unwrapped_user_data->lut_width;
            uint32_t const lut_height = unwrapped_user_data->lut_height;
            float const *const ndotv_one_M2_m11 = unwrapped_user_data->ndotv_one_M2_m11;
            float const *const ndotv_one_M2_m22 = unwrapped_user_data->ndotv_one_M2_m22;

            uint32_t const lut_height_index = (begin / (lut_width - 1U));

            assert((begin + (lut_width - 1U)) == end);

            // use previous configuration (of the same roughness) as first guess
            float M2_m11 = ndotv_one_M2_m11[lut_height_index];
            float M2_m22 = ndotv_one_M2_m22[lut_height_index];
            float M2_m13 = 0.0F;

            for (uint32_t lut_width_index = 1U; lut_width_index < lut_width; ++lut_width_index)
            {
                DirectX::XMFLOAT4 packed_M;

                internal_brx_linearly_transformed_cosine_optimization_compute_matrix(lut_width_index, lut_height_index, lut_width, lut_height, M2_m11, M2_m22, M2_m13, packed_M);

                out_lut[lut_width * lut_height_index + lut_width_index] = packed_M;
            }
        },
        &user_data);
}

static inline DirectX::XMFLOAT3 internal_brx_linearly_transformed_cosine_sample_omega(DirectX::XMFLOAT2 const &xi, DirectX::XMFLOAT3X3 const &M)
{
    // PBR Book V3: [13.6.3 Cosine-Weighted Hemisphere Sampling](https://www.pbr-book.org/3ed-2018/Monte_Carlo_Integration/2D_Sampling_with_Multidimensional_Transformations#Cosine-WeightedHemisphereSampling)
    // PBRT-V3: [CosineSampleHemisphere](https://github.com/mmp/pbrt-v3/blob/book/src/core/sampling.h#L155)
    // PBR Book V4: [A.5.3 Cosine-Weighted Hemisphere Sampling](https://www.pbr-book.org/4ed/Sampling_Algorithms/Sampling_Multidimensional_Functions#Cosine-WeightedHemisphereSampling)
    // PBRT-V4: [SampleCosineHemisphere](https://github.com/mmp/pbrt-v4/blob/master/src/pbrt/util/sampling.h#L409)
    // UE4: [CosineSampleHemisphere](https://github.com/EpicGames/UnrealEngine/blob/4.27/Engine/Shaders/Private/MonteCarlo.ush#L241)
    // U3D: [SampleHemisphereCosine](https://github.com/Unity-Technologies/Graphics/blob/v10.8.1/com.unity.render-pipelines.core/ShaderLibrary/Sampling/Sampling.hlsl#L157)

    // TODO: "ConcentricSampleDisk" does NOT work, but "UniformSampleDisk" works

    DirectX::XMFLOAT2 d;
    {
#if 0
        // Map uniform random numbers to $[-1, 1]^2$
        DirectX::XMFLOAT2 u_offset(xi.x * 2.0F - 1.0F, xi.y * 2.0F - 1.0F);

        if (0.0F == u_offset.x && 0.0F == u_offset.y)
        {
            // Handle degeneracy at the origin
            d = DirectX::XMFLOAT2(0.0F, 0.0F);
        }
        else
        {
            // Apply concentric mapping to point
            float r;
            float theta;

            if (std::abs(u_offset.x) > std::abs(u_offset.y))
            {
                r = u_offset.x;
                theta = (DirectX::XM_PI / 4.0F) * (u_offset.y / u_offset.x);
            }
            else
            {
                r = u_offset.y;
                theta = (DirectX::XM_PI / 2.0F) - (DirectX::XM_PI / 4.0F) * (u_offset.x / u_offset.y);
            }

#if 0
            float cos_theta;
            float sin_theta;
            DirectX::XMScalarSinCos(&sin_theta, &cos_theta, theta);
#else
            float cos_theta = std::cos(theta);
            float sin_theta = std::sin(theta);
#endif

            d = DirectX::XMFLOAT2(r * cos_theta, r * sin_theta);
        }
#else
        float r = std::sqrt(xi.x);
        float theta = 2.0F * DirectX::XM_PI * xi.y;
#if 0
        float cos_theta;
        float sin_theta;
        DirectX::XMScalarSinCos(&sin_theta, &cos_theta, theta);
#else
        float cos_theta = std::cos(theta);
        float sin_theta = std::sin(theta);
#endif

        d = DirectX::XMFLOAT2(r * cos_theta, r * sin_theta);
#endif
    }

    DirectX::XMFLOAT3 omega;
    {
        float z = std::sqrt(std::max(0.0F, 1.0F - DirectX::XMVectorGetX(DirectX::XMVector2Dot(DirectX::XMLoadFloat2(&d), DirectX::XMLoadFloat2(&d)))));

        omega = DirectX::XMFLOAT3(d.x, d.y, z);
    }

    DirectX::XMFLOAT3 linearly_transformed_omega;
    {
        DirectX::XMStoreFloat3(&linearly_transformed_omega, DirectX::XMVector3Normalize(DirectX::XMVector3TransformNormal(DirectX::XMLoadFloat3(&omega), DirectX::XMLoadFloat3x3(&M))));
    }

    return linearly_transformed_omega;
}

static inline float internal_brx_linearly_transformed_cosine_pdf(DirectX::XMFLOAT3X3 const &M, DirectX::XMFLOAT3X3 const &inverse_M, float const absolute_determinant_M, DirectX::XMFLOAT3 const &linearly_transformed_omega)
{
    DirectX::XMFLOAT3 omega;
    {
        DirectX::XMStoreFloat3(&omega, DirectX::XMVector3Normalize(DirectX::XMVector3TransformNormal(DirectX::XMLoadFloat3(&linearly_transformed_omega), DirectX::XMLoadFloat3x3(&inverse_M))));
    }

    float D = (1.0F / DirectX::XM_PI) * std::max(0.0F, omega.z);

    float Jacobian;
    {
        // reciprocal of "Equation 18" of \[Heitz 2016\] [Eric Heitz, Jonathan Dupuy, Stephen Hill, David Neubelt. "Real-Time Polygonal-Light Shading with Linearly Transformed Cosines." SIGGRAPH 2016.](https://eheitzresearch.wordpress.com/415-2/)

        float linearly_transformed_omega_length = DirectX::XMVectorGetX(DirectX::XMVector3Length(DirectX::XMVector3TransformNormal(DirectX::XMLoadFloat3(&omega), DirectX::XMLoadFloat3x3(&M))));

        Jacobian = (linearly_transformed_omega_length * linearly_transformed_omega_length * linearly_transformed_omega_length) / absolute_determinant_M;
    }

    return D * Jacobian;
}

static inline void internal_brx_linearly_transformed_cosine_optimization_compute_monochromatic_brdf_cosine_spherical_moments(float const raw_alpha, DirectX::XMFLOAT3 const &raw_omega_o, float &out_norm, DirectX::XMFLOAT3 &out_median_vector)
{
    assert(raw_alpha >= BRX_TROWBRIDGE_REITZ_ALPHA_MINIMUM);
    float const alpha = std::max(BRX_TROWBRIDGE_REITZ_ALPHA_MINIMUM, raw_alpha);

    assert(raw_omega_o.z >= BRX_TROWBRIDGE_REITZ_NDOTV_MINIMUM);
    DirectX::XMFLOAT3 const omega_o(raw_omega_o.x, raw_omega_o.y, std::max(BRX_TROWBRIDGE_REITZ_NDOTV_MINIMUM, raw_omega_o.z));

    struct parallel_reduce_user_data
    {
        DirectX::XMFLOAT3 const omega_o;
        float const alpha;
    };

    parallel_reduce_user_data user_data = {omega_o, alpha};

    mcrt_double4 norm_and_median_vector = mcrt_parallel_reduce_double4(
        0U,
        INTERNAL_BRX_LTC_OPTIMIZATION_MONTE_CARLO_SAMPLE_COUNT,
        INTERNAL_BRX_LTC_OPTIMIZATION_MONTE_CARLO_GRAIN_SIZE,
        [](uint32_t begin, uint32_t end, void *wrapped_user_data) -> mcrt_double4
        {
            parallel_reduce_user_data *unwrapped_user_data = static_cast<parallel_reduce_user_data *>(wrapped_user_data);
            DirectX::XMFLOAT3 const omega_o = unwrapped_user_data->omega_o;
            float const alpha = unwrapped_user_data->alpha;

            double norm = 0.0;
            double median_vector_x = 0.0;
            double median_vector_y = 0.0;
            double median_vector_z = 0.0;

            for (uint32_t sample_index = begin; sample_index < end; ++sample_index)
            {
                DirectX::XMFLOAT2 xi = brx_hammersley_2d(sample_index, INTERNAL_BRX_LTC_OPTIMIZATION_MONTE_CARLO_SAMPLE_COUNT);

                DirectX::XMFLOAT3 omega_h = brx_trowbridge_reitz_sample_omega_h(xi, alpha, omega_o);

                DirectX::XMFLOAT3 omega_i;
                {
                    DirectX::XMStoreFloat3(&omega_i, DirectX::XMVector3Normalize(DirectX::XMVector3Reflect(DirectX::XMVectorNegate(DirectX::XMLoadFloat3(&omega_o)), DirectX::XMLoadFloat3(&omega_h))));
                }

                float NdotL = std::max(0.0F, omega_i.z);

                float monochromatic_throughput = brx_trowbridge_reitz_throughput_without_fresnel(alpha, NdotL);

                assert(monochromatic_throughput >= 0.0F);
                norm += ((1.0 / static_cast<double>(INTERNAL_BRX_LTC_OPTIMIZATION_MONTE_CARLO_SAMPLE_COUNT)) * static_cast<double>(monochromatic_throughput));
                median_vector_x += ((1.0 / static_cast<double>(INTERNAL_BRX_LTC_OPTIMIZATION_MONTE_CARLO_SAMPLE_COUNT)) * static_cast<double>(monochromatic_throughput) * static_cast<double>(omega_i.x));
                median_vector_y += ((1.0 / static_cast<double>(INTERNAL_BRX_LTC_OPTIMIZATION_MONTE_CARLO_SAMPLE_COUNT)) * static_cast<double>(monochromatic_throughput) * static_cast<double>(omega_i.y));
                median_vector_z += ((1.0 / static_cast<double>(INTERNAL_BRX_LTC_OPTIMIZATION_MONTE_CARLO_SAMPLE_COUNT)) * static_cast<double>(monochromatic_throughput) * static_cast<double>(omega_i.z));
            }

            return mcrt_double4{norm, median_vector_x, median_vector_y, median_vector_z};
        },
        &user_data);

    // For isotropic BRDF, y component must be zero
    assert(std::abs(norm_and_median_vector.z / std::sqrt(norm_and_median_vector.y * norm_and_median_vector.y + norm_and_median_vector.z * norm_and_median_vector.z + norm_and_median_vector.w * norm_and_median_vector.w)) < 1E-3F);

    out_norm = static_cast<float>(norm_and_median_vector.x);

    {
        DirectX::XMFLOAT3 median_vector_non_unit(static_cast<float>(norm_and_median_vector.y), 0.0F, static_cast<float>(norm_and_median_vector.w));
        DirectX::XMStoreFloat3(&out_median_vector, DirectX::XMVector3Normalize(DirectX::XMLoadFloat3(&median_vector_non_unit)));
    }
}

static inline float internal_brx_linearly_transformed_cosine_optimization_compute_monochromatic_brdf_cosine_mean_cubed_error(float const raw_alpha, DirectX::XMFLOAT3 const &raw_omega_o, float monochromatic_brdf_cosine_norm, DirectX::XMFLOAT3X3 const &ltc_M)
{
    assert(raw_alpha >= BRX_TROWBRIDGE_REITZ_ALPHA_MINIMUM);
    float const alpha = std::max(BRX_TROWBRIDGE_REITZ_ALPHA_MINIMUM, raw_alpha);

    assert(raw_omega_o.z >= BRX_TROWBRIDGE_REITZ_NDOTV_MINIMUM);
    DirectX::XMFLOAT3 const omega_o(raw_omega_o.x, raw_omega_o.y, std::max(BRX_TROWBRIDGE_REITZ_NDOTV_MINIMUM, raw_omega_o.z));

    DirectX::XMFLOAT3X3 inverse_ltc_M;
    float absolute_determinant_ltc_M;
    {
        DirectX::XMVECTOR determinant_ltc_M;
        DirectX::XMStoreFloat3x3(&inverse_ltc_M, DirectX::XMMatrixInverse(&determinant_ltc_M, DirectX::XMLoadFloat3x3(&ltc_M)));

        absolute_determinant_ltc_M = std::abs(DirectX::XMVectorGetX(determinant_ltc_M));
    }

    struct parallel_reduce_user_data
    {
        float const alpha;
        DirectX::XMFLOAT3 const omega_o;
        float const monochromatic_brdf_cosine_norm;
        DirectX::XMFLOAT3X3 const ltc_M;
        DirectX::XMFLOAT3X3 const inverse_ltc_M;
        float absolute_determinant_ltc_M;
    };

    parallel_reduce_user_data user_data = {alpha, omega_o, monochromatic_brdf_cosine_norm, ltc_M, inverse_ltc_M, absolute_determinant_ltc_M};

    double mean_cubed_error = mcrt_parallel_reduce_double(
        0U,
        INTERNAL_BRX_LTC_OPTIMIZATION_MONTE_CARLO_SAMPLE_COUNT,
        INTERNAL_BRX_LTC_OPTIMIZATION_MONTE_CARLO_GRAIN_SIZE,
        [](uint32_t begin, uint32_t end, void *wrapped_user_data) -> double
        {
            parallel_reduce_user_data *unwrapped_user_data = static_cast<parallel_reduce_user_data *>(wrapped_user_data);
            float const alpha = unwrapped_user_data->alpha;
            DirectX::XMFLOAT3 const omega_o = unwrapped_user_data->omega_o;
            float const monochromatic_brdf_cosine_norm = unwrapped_user_data->monochromatic_brdf_cosine_norm;
            DirectX::XMFLOAT3X3 const &ltc_M = unwrapped_user_data->ltc_M;
            DirectX::XMFLOAT3X3 const &inverse_ltc_M = unwrapped_user_data->inverse_ltc_M;
            float const absolute_determinant_ltc_M = unwrapped_user_data->absolute_determinant_ltc_M;

            double mean_cubed_error = 0.0;

            for (uint32_t sample_index = begin; sample_index < end; ++sample_index)
            {
                DirectX::XMFLOAT2 xi = brx_hammersley_2d(sample_index, INTERNAL_BRX_LTC_OPTIMIZATION_MONTE_CARLO_SAMPLE_COUNT);

                {
                    DirectX::XMFLOAT3 omega_i = internal_brx_linearly_transformed_cosine_sample_omega(xi, ltc_M);

                    float ltc_pdf = internal_brx_linearly_transformed_cosine_pdf(ltc_M, inverse_ltc_M, absolute_determinant_ltc_M, omega_i);

                    float ltc_value = ltc_pdf;

                    float brdf_pdf;

                    float monochromatic_brdf_cosine_value;
                    {
                        float NdotL = omega_i.z;

                        float NdotV = omega_o.z;

                        float NdotH;
                        {
#if 0
                            float VdotL = DirectX::XMVectorGetX(DirectX::XMVector3Dot(DirectX::XMLoadFloat3(&omega_o), DirectX::XMLoadFloat3(&omega_i)));

                            // Real-Time Rendering Fourth Edition / 9.8 BRDF Models for Surface Reflection / [Hammon 2017]
                            // UE4: [Init](https://github.com/EpicGames/UnrealEngine/blob/4.27/Engine/Shaders/Private/BRDF.ush#L31)
                            // U3D: [GetBSDFAngle](https://github.com/Unity-Technologies/Graphics/blob/v10.8.0/com.unity.render-pipelines.core/ShaderLibrary/CommonLighting.hlsl#L361)
                            float invLenH = 1.0F / std::sqrt(std::max(INTERNAL_BRX_LTC_OPTIMIZATION_LENGTH_SQUARE_MINIMUM, 2.0F + 2.0F * VdotL));
                            NdotH = std::min(std::max(0.0F, (NdotL + NdotV) * invLenH), 1.0F);
#else

                            // [BrdfGGX::eval](https://github.com/selfshadow/ltc_code/blob/master/fit/brdf_ggx.h#L34)
                            // The precision is higher when calculated based on slope, but why???

                            DirectX::XMFLOAT3 omega_h;
                            DirectX::XMStoreFloat3(&omega_h, DirectX::XMVector3Normalize(DirectX::XMVectorAdd(DirectX::XMLoadFloat3(&omega_o), DirectX::XMLoadFloat3(&omega_i))));

                            if (omega_h.z > 0.0F)
                            {
                                NdotH = std::sqrt(1.0F / (1.0F + (omega_h.x / omega_h.z) * (omega_h.x / omega_h.z) + (omega_h.y / omega_h.z) * (omega_h.y / omega_h.z)));
                            }
                            else
                            {
                                NdotH = 0.0F;
                            }
#endif
                        }

                        NdotL = std::max(0.0F, NdotL);

                        NdotV = std::max(0.0F, NdotV);

                        brdf_pdf = brx_trowbridge_reitz_pdf_omega_i(alpha, NdotV, NdotH);

                        monochromatic_brdf_cosine_value = brx_trowbridge_reitz_brdf_without_fresnel(alpha, NdotH, NdotV, NdotL) * NdotL;
                    }

                    float absolute_distance = std::abs(monochromatic_brdf_cosine_value - monochromatic_brdf_cosine_norm * ltc_value);

                    assert(ltc_pdf >= 0.0F);
                    assert(ltc_value >= 0.0F);
                    assert(brdf_pdf >= 0.0F);
                    assert(monochromatic_brdf_cosine_value >= 0.0F);
                    assert((brdf_pdf + ltc_pdf) > 0.0F);
                    mean_cubed_error += ((1.0 / static_cast<double>(INTERNAL_BRX_LTC_OPTIMIZATION_MONTE_CARLO_SAMPLE_COUNT)) * (static_cast<double>(absolute_distance) * static_cast<double>(absolute_distance) * static_cast<double>(absolute_distance)) / (static_cast<double>(ltc_pdf) + static_cast<double>(brdf_pdf)));
                }

                {
                    DirectX::XMFLOAT3 omega_h = brx_trowbridge_reitz_sample_omega_h(xi, alpha, omega_o);

                    DirectX::XMFLOAT3 omega_i;
                    {
                        DirectX::XMStoreFloat3(&omega_i, DirectX::XMVector3Normalize(DirectX::XMVector3Reflect(DirectX::XMVectorNegate(DirectX::XMLoadFloat3(&omega_o)), DirectX::XMLoadFloat3(&omega_h))));
                    }

                    float ltc_pdf = internal_brx_linearly_transformed_cosine_pdf(ltc_M, inverse_ltc_M, absolute_determinant_ltc_M, omega_i);

                    float ltc_value = ltc_pdf;

                    float brdf_pdf;

                    float monochromatic_brdf_cosine_value;
                    {
                        float NdotL = std::max(0.0F, omega_i.z);

                        float NdotV = std::max(0.0F, omega_o.z);

                        float NdotH = std::max(0.0F, omega_h.z);

                        brdf_pdf = brx_trowbridge_reitz_pdf_omega_i(alpha, NdotV, NdotH);

                        monochromatic_brdf_cosine_value = brx_trowbridge_reitz_brdf_without_fresnel(alpha, NdotH, NdotV, NdotL) * NdotL;
                    }

                    float absolute_distance = std::abs(monochromatic_brdf_cosine_value - monochromatic_brdf_cosine_norm * ltc_value);

                    assert(ltc_pdf >= 0.0F);
                    assert(ltc_value >= 0.0F);
                    assert(brdf_pdf >= 0.0F);
                    assert(monochromatic_brdf_cosine_value >= 0.0F);
                    assert((brdf_pdf + ltc_pdf) > 0.0F);
                    mean_cubed_error += ((1.0 / static_cast<double>(INTERNAL_BRX_LTC_OPTIMIZATION_MONTE_CARLO_SAMPLE_COUNT)) * (static_cast<double>(absolute_distance) * static_cast<double>(absolute_distance) * static_cast<double>(absolute_distance)) / (static_cast<double>(ltc_pdf) + static_cast<double>(brdf_pdf)));
                }
            }
            return mean_cubed_error;
        },
        &user_data);

    return static_cast<float>(mean_cubed_error);
}

static inline float internal_brx_linearly_transformed_cosine_optimization_expit(float const x)
{
    return static_cast<float>(1.0 / (1.0 + std::exp(static_cast<double>(-x))));
}

static inline float internal_brx_linearly_transformed_cosine_optimization_logit(float const p)
{
    return static_cast<float>(std::log(static_cast<double>(p)) - std::log1p(-static_cast<double>(p)));
}

static inline void internal_brx_linearly_transformed_cosine_optimization_compute_matrix(uint32_t const lut_width_index, uint32_t const lut_height_index, uint32_t const lut_width, uint32_t const lut_height, float &inout_M2_m11, float &inout_M2_m22, float &inout_M2_m13, DirectX::XMFLOAT4 &out_packed_M)
{
    // Remap: [0, 1] -> [0.5/size, 1.0 - 0.5/size]
    // U3D: [Remap01ToHalfTexelCoord](https://github.com/Unity-Technologies/Graphics/blob/v10.8.0/com.unity.render-pipelines.core/ShaderLibrary/Common.hlsl#L661)
    // UE4: [N/A](https://github.com/EpicGames/UnrealEngine/blob/4.27/Engine/Shaders/Private/RectLight.ush#L450)

    assert(lut_width_index < lut_width);
    float texcoord_u = static_cast<float>(lut_width_index) / static_cast<float>(lut_width - 1U);

    assert(lut_height_index < lut_height);
    float texcoord_v = static_cast<float>(lut_height_index) / static_cast<float>(lut_height - 1U);

    // u = 1 - cos_theta
    // cos_theta = 1 - u
    DirectX::XMFLOAT3 omega_o;
    {
        float const cos_theta_o = std::max(BRX_TROWBRIDGE_REITZ_NDOTV_MINIMUM, 1.0F - texcoord_u);
        omega_o = DirectX::XMFLOAT3(std::sqrt(std::max(0.0F, 1.0F - cos_theta_o * cos_theta_o)), 0.0F, cos_theta_o);
    }

    // v = 1 - alpha
    // alpha = 1 - v
    float const alpha = std::max(BRX_TROWBRIDGE_REITZ_ALPHA_MINIMUM, 1.0F - texcoord_v);

    float norm;
    DirectX::XMFLOAT3 median_vector;
    internal_brx_linearly_transformed_cosine_optimization_compute_monochromatic_brdf_cosine_spherical_moments(alpha, omega_o, norm, median_vector);

    assert(0.0F == median_vector.y);
    assert(0.0F == omega_o.y);

    DirectX::XMFLOAT3 T1;
    DirectX::XMFLOAT3 T2;
    {
        DirectX::XMVECTOR simd_omega_o = DirectX::XMLoadFloat3(&omega_o);
        DirectX::XMVECTOR simd_median_vector = DirectX::XMLoadFloat3(&median_vector);
        DirectX::XMVECTOR simd_T1_non_unit = DirectX::XMVectorSubtract(simd_omega_o, DirectX::XMVectorScale(simd_median_vector, DirectX::XMVectorGetX(DirectX::XMVector3Dot(simd_omega_o, simd_median_vector))));
        float T1_length_square = DirectX::XMVectorGetX(DirectX::XMVector3Dot(simd_T1_non_unit, simd_T1_non_unit));
        if (T1_length_square > BRX_TROWBRIDGE_REITZ_TANGENT_SPACE_T1_LENGTH_SQUARE_MINIMUM)
        {
            DirectX::XMStoreFloat3(&T1, DirectX::XMVectorScale(simd_T1_non_unit, 1.0F / std::sqrt(T1_length_square)));
        }
        else
        {
            T1 = DirectX::XMFLOAT3(1.0F, 0.0F, 0.0F);
        }

        DirectX::XMStoreFloat3(&T2, DirectX::XMVector3Normalize(DirectX::XMVector3Cross(simd_median_vector, DirectX::XMLoadFloat3(&T1))));
    }

    // if theta == 0 the lobe is rotationally symmetric and aligned with Z = (0 0 1)
    // if (lut_height_index == 0)
    // assert(std::abs(X.x - 1.0F) < 1E-3F);
    // assert(std::abs(X.y - 0.0F) < 1E-3F);
    // assert(std::abs(X.z - 0.0F) < 1E-3F);

    // assert(std::abs(Y.x - 0.0F) < 1E-3F);
    // assert(std::abs(Y.y - 1.0F) < 1E-3F);
    // assert(std::abs(Y.z - 0.0F) < 1E-3F);

    // assert(std::abs(Z.x - 0.0F) < 1E-3F);
    // assert(std::abs(Z.y - 0.0F) < 1E-3F);
    // assert(std::abs(Z.z - 1.0F) < 1E-3F);

    DirectX::XMFLOAT3X3 const M1 = {
        T1.x, T1.y, T1.z,
        T2.x, T2.y, T2.z,
        median_vector.x, median_vector.y, median_vector.z};

    struct fminsearch_3d_user_data
    {
        float const alpha;
        DirectX::XMFLOAT3 const omega_o;

        float const norm;
        DirectX::XMFLOAT3X3 const M1;
    };

    fminsearch_3d_user_data user_data = {alpha, omega_o, norm, M1};

    DirectX::XMFLOAT3 inout_x = {2.0F * internal_brx_linearly_transformed_cosine_optimization_logit(inout_M2_m11 * 0.5F), 2.0F * internal_brx_linearly_transformed_cosine_optimization_logit(inout_M2_m22 * 0.5F), 2.0F * internal_brx_linearly_transformed_cosine_optimization_logit((1.0F + inout_M2_m13) * 0.5F)};

    brx_fminsearch_3d(
        [](DirectX::XMFLOAT3 const &x, void *wrapped_user_data) -> float
        {
            fminsearch_3d_user_data *unwrapped_user_data = static_cast<fminsearch_3d_user_data *>(wrapped_user_data);
            float const alpha = unwrapped_user_data->alpha;
            DirectX::XMFLOAT3 const omega_o = unwrapped_user_data->omega_o;
            float const norm = unwrapped_user_data->norm;
            DirectX::XMFLOAT3X3 const M1 = unwrapped_user_data->M1;

            float M2_m11 = 2.0F * internal_brx_linearly_transformed_cosine_optimization_expit(x.x * 0.5F);
            float M2_m22 = 2.0F * internal_brx_linearly_transformed_cosine_optimization_expit(x.y * 0.5F);
            float M2_m13 = 2.0F * internal_brx_linearly_transformed_cosine_optimization_expit(x.z * 0.5F) - 1.0F;

            DirectX::XMFLOAT3X3 const M2 = {
                M2_m11, 0.0F, 0.0F,
                0.0F, M2_m22, 0.0F,
                M2_m13, 0.0F, 1.0F};

            DirectX::XMFLOAT3X3 M;
            DirectX::XMStoreFloat3x3(&M, DirectX::XMMatrixMultiply(DirectX::XMLoadFloat3x3(&M2), DirectX::XMLoadFloat3x3(&M1)));

            return internal_brx_linearly_transformed_cosine_optimization_compute_monochromatic_brdf_cosine_mean_cubed_error(alpha, omega_o, norm, M);
        },
        &user_data,
        inout_x);

    // Update LTC with best fitting values
    inout_M2_m11 = 2.0F * internal_brx_linearly_transformed_cosine_optimization_expit(inout_x.x * 0.5F);
    inout_M2_m22 = 2.0F * internal_brx_linearly_transformed_cosine_optimization_expit(inout_x.y * 0.5F);
    inout_M2_m13 = 2.0F * internal_brx_linearly_transformed_cosine_optimization_expit(inout_x.z * 0.5F) - 1.0F;

    DirectX::XMFLOAT3X3 const M2 = {
        inout_M2_m11, 0.0F, 0.0F,
        0.0F, inout_M2_m22, 0.0F,
        inout_M2_m13, 0.0F, 1.0F};

    DirectX::XMFLOAT3X3 inverse_M;
    {
        DirectX::XMMATRIX M = DirectX::XMMatrixMultiply(DirectX::XMLoadFloat3x3(&M2), DirectX::XMLoadFloat3x3(&M1));

        DirectX::XMVECTOR unused_determinant;
        DirectX::XMStoreFloat3x3(&inverse_M, DirectX::XMMatrixInverse(&unused_determinant, M));
    }

    // normalize by the middle element
    // ltc scale invariant
    assert(std::abs(inverse_M.m[1][1]) >= INTERNAL_BRX_LTC_OPTIMIZATION_SCALE_MINIMUM);
    float reciprocal_inverse_M_m22 = 1.0F / inverse_M.m[1][1];

    out_packed_M.x = inverse_M.m[0][0] * reciprocal_inverse_M_m22;
    out_packed_M.y = inverse_M.m[0][2] * reciprocal_inverse_M_m22;
    out_packed_M.z = inverse_M.m[2][0] * reciprocal_inverse_M_m22;
    out_packed_M.w = inverse_M.m[2][2] * reciprocal_inverse_M_m22;
}

#endif
