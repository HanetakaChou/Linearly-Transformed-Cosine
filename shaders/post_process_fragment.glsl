#version 310 es

#extension GL_GOOGLE_include_directive : enable

#define USE_FRAMEBUFFER_SET 1
#include "unified_pipeline_layout.h"

layout(location = 0) out highp vec4 fragment_output_color;

highp vec3 aces_fitted(highp vec3 color);

highp vec3 ToSRGB(highp vec3 v);

void main()
{
    highp vec3 col = subpassLoad(framebuffer_set_input_attachment_backup).rgb;

    col = aces_fitted(col);

    col = ToSRGB(col);

    fragment_output_color = vec4(col, 1.0);
}

highp float saturate(highp float v)
{
    return clamp(v, 0.0, 1.0);
}

highp vec3 saturate(highp vec3 v)
{
    return vec3(saturate(v.x), saturate(v.y), saturate(v.z));
}

highp vec3 rrt_odt_fit(highp vec3 v)
{
    highp vec3 a = v * (v + 0.0245786) - 0.000090537;
    highp vec3 b = v * (0.983729 * v + 0.4329510) + 0.238081;
    return a / b;
}

highp vec3 aces_fitted(highp vec3 color)
{
    highp mat3x3 ACES_INPUT_MAT = mat3x3(
        vec3(0.59719, 0.07600, 0.02840),
        vec3(0.35458, 0.90834, 0.13383),
        vec3(0.04823, 0.01566, 0.83777));

    highp mat3x3 ACES_OUTPUT_MAT = mat3x3(
        vec3(1.60475, -0.10208, -0.00327),
        vec3(-0.53108, 1.10813, -0.07276),
        vec3(-0.07367, -0.00605, 1.07602));

    color = ACES_INPUT_MAT * color;

    // Apply RRT and ODT
    color = rrt_odt_fit(color);

    color = ACES_OUTPUT_MAT * color;

    // Clamp to [0, 1]
    color = saturate(color);

    return color;
}

highp vec3 PowVec3(highp vec3 v, highp float p)
{
    return vec3(pow(v.x, p), pow(v.y, p), pow(v.z, p));
}

highp vec3 ToSRGB(highp vec3 v)
{
    const highp float gamma = 2.2;
    return PowVec3(v, 1.0 / gamma);
}