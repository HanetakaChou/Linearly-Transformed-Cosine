#version 310 es

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

#extension GL_GOOGLE_include_directive : enable

#include "unified_pipeline_layout.h"

layout(location = 0) in highp vec3 vertex_input_position;
layout(location = 1) in highp vec3 vertex_input_normal;

layout(location = 0) out highp vec3 vertex_output_position;
layout(location = 1) out highp vec3 vertex_output_normal;

void main()
{
    highp vec3 model_position = vertex_input_position;
    highp vec3 world_position = (model_transform * vec4(model_position, 1.0)).xyz;
    highp vec4 clip_position = projection_transform * view_transform * vec4(world_position, 1.0);

    highp vec3 model_normal = vertex_input_normal;
    // TODO: normal transform
    highp mat3x3 tangent_transform = mat3x3(model_transform[0].xyz, model_transform[1].xyz, model_transform[2].xyz);
    highp vec3 world_normal = normalize(tangent_transform * model_normal);

    vertex_output_position = world_position;
    vertex_output_normal = world_normal;
    gl_Position = clip_position;
}
