#ifndef _FORWARD_SHADING_PIPELINE_LAYOUT_H_
#define _FORWARD_SHADING_PIPELINE_LAYOUT_H_ 1

#if defined(__STDC__) || defined(__cplusplus)

struct forward_shading_layout_global_set_frame_uniform_buffer_binding
{
    // camera
    DirectX::XMFLOAT4X4 view_transform;
    DirectX::XMFLOAT4X4 projection_transform;
    DirectX::XMFLOAT3 eye_position;
    float _padding_eye_position;

    // light
    DirectX::XMFLOAT4 rect_light_vetices[4];
    float intensity;
    float culling_range;
};

struct forward_shading_layout_global_set_object_uniform_buffer_binding
{
    DirectX::XMFLOAT4X4 model_transform;
};

struct forward_shading_layout_material_set_uniform_buffer_binding
{
    // mesh
    DirectX::XMFLOAT3 dcolor;
    float _padding_dcolor;
    DirectX::XMFLOAT3 scolor;
    float _padding_scolor;
    float roughness;
    float _padding_roughness_1;
    float _padding_roughness_2;
    float _padding_roughness_3;
};

#elif defined(GL_SPIRV) || defined(VULKAN)

layout(set = 0, binding = 0, column_major) uniform _global_set_frame_uniform_buffer_binding
{
    // camera
    highp mat4x4 view_transform;
    highp mat4x4 projection_transform;
    highp vec3 eye_position;
    highp float _padding_eye_position;

    // light
    highp vec4 rect_light_vetices[4];
    highp float intensity;
    highp float culling_range;
};

layout(set = 0, binding = 1, column_major) uniform _global_set_object_uniform_buffer_binding
{
    highp mat4x4 model_transform;
};

layout(set = 0, binding = 2) uniform highp sampler2D ltc_matrix_lut;

layout(set = 0, binding = 3) uniform highp sampler2D preintegrated_hdr_lut;

layout(set = 1, binding = 0, column_major) uniform _material_set_uniform_buffer_binding
{
    // mesh
    highp vec3 dcolor;
    highp float _padding_dcolor;
    highp vec3 scolor;
    highp float _padding_scolor;
    highp float roughness;
    highp float _padding_roughness_1;
    highp float _padding_roughness_2;
    highp float _padding_roughness_3;
};

#if defined(GL_FRAGMENT_SHADER) || defined(USE_FRAMEBUFFER_SET)
layout(set = 2, binding = 0, input_attachment_index = 0) uniform highp subpassInput framebuffer_set_input_attachment_backup;
#endif

#else
#error Unknown Compiler
#endif

#endif