
#include <sdkddkver.h>
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX 1
#include <windows.h>
#include <stdint.h>
#include <assert.h>
#include <cmath>
#include <algorithm>
#include <DirectXMath.h>
#include <dxgi.h>
#include <d3d11.h>
#include "support_d3d11/resolution.h"
#include "support_d3d11/camera_controller.h"
#include "demo_d3d11.h"
#include "../thirdparty/Linearly-Transformed-Cosine/include/brx_linearly_transformed_cosine_look_up_table_matrices.h"
#include "../thirdparty/Hemispherical-Directional-Reflectance/include/brx_hemispherical_directional_reflectance_look_up_table_norms.h"
#include "../shaders/plane_vs.hlsl.inl"
#include "../shaders/plane_fs.hlsl.inl"
#include "../shaders/rect_light_vs.hlsl.inl"
#include "../shaders/rect_light_fs.hlsl.inl"
#include "../shaders/post_process_vs.hlsl.inl"
#include "../shaders/post_process_fs.hlsl.inl"
#include "../assets/plane.h"

struct plane_uniform_buffer_per_frame_binding_t
{
	// mesh
	DirectX::XMFLOAT4X4 model_transform;
	DirectX::XMFLOAT3 dcolor;
	float _padding_dcolor;
	DirectX::XMFLOAT3 scolor;
	float _padding_scolor;
	float roughness;
	float _padding_roughness_1;
	float _padding_roughness_2;
	float _padding_roughness_3;

	// camera
	DirectX::XMFLOAT4X4 view_transform;
	DirectX::XMFLOAT4X4 projection_transform;
	DirectX::XMFLOAT3 eye_position;
	float _padding_eye_position;

	// light
	DirectX::XMFLOAT4 rect_light_vetices[4];
	float intensity;
	float culling_range;
	float __padding_align16_1;
	float __padding_align16_2;
};

struct rect_light_uniform_buffer_per_frame_binding_t
{
	// camera
	DirectX::XMFLOAT4X4 view_transform;
	DirectX::XMFLOAT4X4 projection_transform;

	// light
	DirectX::XMFLOAT4 rect_light_vetices[4];
	float intensity;
	float __padding_align16_1;
	float __padding_align16_2;
	float __padding_align16_3;
};

void Demo::Init(ID3D11Device *d3d_device, ID3D11DeviceContext *d3d_device_context, IDXGISwapChain *dxgi_swap_chain)
{
	m_attachment_backbuffer_rtv = NULL;
	{
		ID3D11Texture2D *attachment_backbuffer = NULL;
		HRESULT res_dxgi_swap_chain_get_buffer = dxgi_swap_chain->GetBuffer(0U, IID_PPV_ARGS(&attachment_backbuffer));
		assert(SUCCEEDED(res_dxgi_swap_chain_get_buffer));

		D3D11_RENDER_TARGET_VIEW_DESC d3d_render_target_view_desc;
		d3d_render_target_view_desc.Format = DXGI_FORMAT_B8G8R8A8_UNORM;
		d3d_render_target_view_desc.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE2D;
		d3d_render_target_view_desc.Texture2D.MipSlice = 0U;

		HRESULT res_d3d_device_create_render_target_view = d3d_device->CreateRenderTargetView(attachment_backbuffer, &d3d_render_target_view_desc, &m_attachment_backbuffer_rtv);
		assert(SUCCEEDED(res_d3d_device_create_render_target_view));
	}

	m_plane_vb_position = NULL;
	{

		D3D11_BUFFER_DESC d3d_buffer_desc;
		d3d_buffer_desc.ByteWidth = sizeof(plane_vertex_postion);
		d3d_buffer_desc.Usage = D3D11_USAGE_DEFAULT;
		d3d_buffer_desc.BindFlags = D3D11_BIND_VERTEX_BUFFER;
		d3d_buffer_desc.CPUAccessFlags = 0U;
		d3d_buffer_desc.MiscFlags = 0U;
		d3d_buffer_desc.StructureByteStride = 0U;

		D3D11_SUBRESOURCE_DATA d3d_subresource_data;
		d3d_subresource_data.pSysMem = plane_vertex_postion;
		d3d_subresource_data.SysMemPitch = sizeof(plane_vertex_postion);
		d3d_subresource_data.SysMemSlicePitch = sizeof(plane_vertex_postion);

		HRESULT res_d3d_device_create_buffer = d3d_device->CreateBuffer(&d3d_buffer_desc, &d3d_subresource_data, &m_plane_vb_position);
		assert(SUCCEEDED(res_d3d_device_create_buffer));
	}

	m_plane_vb_varying = NULL;
	{
		D3D11_BUFFER_DESC d3d_buffer_desc;
		d3d_buffer_desc.ByteWidth = sizeof(plane_vertex_normal);
		d3d_buffer_desc.Usage = D3D11_USAGE_DEFAULT;
		d3d_buffer_desc.BindFlags = D3D11_BIND_VERTEX_BUFFER;
		d3d_buffer_desc.CPUAccessFlags = 0U;
		d3d_buffer_desc.MiscFlags = 0U;
		d3d_buffer_desc.StructureByteStride = 0U;

		D3D11_SUBRESOURCE_DATA d3d_subresource_data;
		d3d_subresource_data.pSysMem = plane_vertex_normal;
		d3d_subresource_data.SysMemPitch = sizeof(plane_vertex_normal);
		d3d_subresource_data.SysMemSlicePitch = sizeof(plane_vertex_normal);

		HRESULT res_d3d_device_create_buffer = d3d_device->CreateBuffer(&d3d_buffer_desc, &d3d_subresource_data, &m_plane_vb_varying);
		assert(SUCCEEDED(res_d3d_device_create_buffer));
	}

	m_plane_vao = NULL;
	{
		D3D11_INPUT_ELEMENT_DESC d3d_input_elements_desc[] =
			{
				{"POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0, 0, D3D11_INPUT_PER_VERTEX_DATA, 0},
				{"NORMAL", 0, DXGI_FORMAT_R32G32B32_FLOAT, 1, 0, D3D11_INPUT_PER_VERTEX_DATA, 0},
			};

		HRESULT res_d3d_create_input_layout = d3d_device->CreateInputLayout(
			d3d_input_elements_desc,
			sizeof(d3d_input_elements_desc) / sizeof(d3d_input_elements_desc[0]),
			plane_vs_bytecode,
			sizeof(plane_vs_bytecode),
			&m_plane_vao);
		assert(SUCCEEDED(res_d3d_create_input_layout));
	}

	m_plane_vs = NULL;
	{
		HRESULT res_d3d_device_create_vertex_shader = d3d_device->CreateVertexShader(plane_vs_bytecode, sizeof(plane_vs_bytecode), NULL, &m_plane_vs);
		assert(SUCCEEDED(res_d3d_device_create_vertex_shader));
	}

	m_plane_fs = NULL;
	{
		HRESULT res_d3d_device_create_pixel_shader = d3d_device->CreatePixelShader(plane_fs_bytecode, sizeof(plane_fs_bytecode), NULL, &m_plane_fs);
		assert(SUCCEEDED(res_d3d_device_create_pixel_shader));
	}

	m_plane_rs = NULL;
	{
		D3D11_RASTERIZER_DESC d3d_rasterizer_desc;
		d3d_rasterizer_desc.FillMode = D3D11_FILL_SOLID;
		d3d_rasterizer_desc.CullMode = D3D11_CULL_BACK;
		d3d_rasterizer_desc.FrontCounterClockwise = TRUE;
		d3d_rasterizer_desc.DepthBias = 0;
		d3d_rasterizer_desc.SlopeScaledDepthBias = 0.0f;
		d3d_rasterizer_desc.DepthBiasClamp = 0.0f;
		d3d_rasterizer_desc.DepthClipEnable = TRUE;
		d3d_rasterizer_desc.ScissorEnable = FALSE;
		d3d_rasterizer_desc.MultisampleEnable = FALSE;
		d3d_rasterizer_desc.AntialiasedLineEnable = FALSE;
		d3d_rasterizer_desc.MultisampleEnable = FALSE;
		HRESULT res_d3d_device_create_rasterizer_state = d3d_device->CreateRasterizerState(&d3d_rasterizer_desc, &m_plane_rs);
		assert(SUCCEEDED(res_d3d_device_create_rasterizer_state));
	}

	m_plane_uniform_buffer_per_frame_binding = NULL;
	{
		D3D11_BUFFER_DESC d3d_buffer_desc;
		d3d_buffer_desc.ByteWidth = sizeof(plane_uniform_buffer_per_frame_binding_t);
		d3d_buffer_desc.Usage = D3D11_USAGE_DEFAULT;
		d3d_buffer_desc.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
		d3d_buffer_desc.CPUAccessFlags = 0U;
		d3d_buffer_desc.MiscFlags = 0U;
		d3d_buffer_desc.StructureByteStride = 0U;

		HRESULT res_d3d_device_create_buffer = d3d_device->CreateBuffer(&d3d_buffer_desc, NULL, &m_plane_uniform_buffer_per_frame_binding);
		assert(SUCCEEDED(res_d3d_device_create_buffer));
	}

	m_rect_light_vao = NULL;

	m_rect_light_vs = NULL;
	{
		HRESULT res_d3d_device_create_vertex_shader = d3d_device->CreateVertexShader(rect_light_vs_bytecode, sizeof(rect_light_vs_bytecode), NULL, &m_rect_light_vs);
		assert(SUCCEEDED(res_d3d_device_create_vertex_shader));
	}

	m_rect_light_fs = NULL;
	{
		HRESULT res_d3d_device_create_pixel_shader = d3d_device->CreatePixelShader(rect_light_fs_bytecode, sizeof(rect_light_fs_bytecode), NULL, &m_rect_light_fs);
		assert(SUCCEEDED(res_d3d_device_create_pixel_shader));
	}

	m_rect_light_rs = NULL;
	{
		D3D11_RASTERIZER_DESC d3d_rasterizer_desc;
		d3d_rasterizer_desc.FillMode = D3D11_FILL_SOLID;
		d3d_rasterizer_desc.CullMode = D3D11_CULL_BACK;
		d3d_rasterizer_desc.FrontCounterClockwise = FALSE;
		d3d_rasterizer_desc.DepthBias = 0;
		d3d_rasterizer_desc.SlopeScaledDepthBias = 0.0f;
		d3d_rasterizer_desc.DepthBiasClamp = 0.0f;
		d3d_rasterizer_desc.DepthClipEnable = TRUE;
		d3d_rasterizer_desc.ScissorEnable = FALSE;
		d3d_rasterizer_desc.MultisampleEnable = FALSE;
		d3d_rasterizer_desc.AntialiasedLineEnable = FALSE;
		d3d_rasterizer_desc.MultisampleEnable = FALSE;
		HRESULT res_d3d_device_create_rasterizer_state = d3d_device->CreateRasterizerState(&d3d_rasterizer_desc, &m_rect_light_rs);
		assert(SUCCEEDED(res_d3d_device_create_rasterizer_state));
	}

	m_rect_light_uniform_buffer_per_frame_binding = NULL;
	{
		D3D11_BUFFER_DESC d3d_buffer_desc;
		d3d_buffer_desc.ByteWidth = sizeof(rect_light_uniform_buffer_per_frame_binding_t);
		d3d_buffer_desc.Usage = D3D11_USAGE_DEFAULT;
		d3d_buffer_desc.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
		d3d_buffer_desc.CPUAccessFlags = 0U;
		d3d_buffer_desc.MiscFlags = 0U;
		d3d_buffer_desc.StructureByteStride = 0U;

		HRESULT res_d3d_device_create_buffer = d3d_device->CreateBuffer(&d3d_buffer_desc, NULL, &m_rect_light_uniform_buffer_per_frame_binding);
		assert(SUCCEEDED(res_d3d_device_create_buffer));
	}

	m_post_process_vao = NULL;

	m_post_process_vs = NULL;
	{
		HRESULT res_d3d_device_create_vertex_shader = d3d_device->CreateVertexShader(post_process_vs_bytecode, sizeof(post_process_vs_bytecode), NULL, &m_post_process_vs);
		assert(SUCCEEDED(res_d3d_device_create_vertex_shader));
	}

	m_post_process_fs = NULL;
	{
		HRESULT res_d3d_device_create_pixel_shader = d3d_device->CreatePixelShader(post_process_fs_bytecode, sizeof(post_process_fs_bytecode), NULL, &m_post_process_fs);
		assert(SUCCEEDED(res_d3d_device_create_pixel_shader));
	}

	m_post_process_rs = NULL;
	{
		D3D11_RASTERIZER_DESC d3d_rasterizer_desc;
		d3d_rasterizer_desc.FillMode = D3D11_FILL_SOLID;
		d3d_rasterizer_desc.CullMode = D3D11_CULL_BACK;
		d3d_rasterizer_desc.FrontCounterClockwise = FALSE;
		d3d_rasterizer_desc.DepthBias = 0;
		d3d_rasterizer_desc.SlopeScaledDepthBias = 0.0f;
		d3d_rasterizer_desc.DepthBiasClamp = 0.0f;
		d3d_rasterizer_desc.DepthClipEnable = TRUE;
		d3d_rasterizer_desc.ScissorEnable = FALSE;
		d3d_rasterizer_desc.MultisampleEnable = FALSE;
		d3d_rasterizer_desc.AntialiasedLineEnable = FALSE;
		d3d_rasterizer_desc.MultisampleEnable = FALSE;
		HRESULT res_d3d_device_create_rasterizer_state = d3d_device->CreateRasterizerState(&d3d_rasterizer_desc, &m_post_process_rs);
		assert(SUCCEEDED(res_d3d_device_create_rasterizer_state));
	}

	m_ltc_lut_sampler = NULL;
	{
		D3D11_SAMPLER_DESC d3d_sampler_desc;
		d3d_sampler_desc.Filter = D3D11_FILTER_MIN_POINT_MAG_LINEAR_MIP_POINT;
		d3d_sampler_desc.AddressU = D3D11_TEXTURE_ADDRESS_CLAMP;
		d3d_sampler_desc.AddressV = D3D11_TEXTURE_ADDRESS_CLAMP;
		d3d_sampler_desc.AddressW = D3D11_TEXTURE_ADDRESS_CLAMP;
		d3d_sampler_desc.MipLODBias = 0.0;
		d3d_sampler_desc.MaxAnisotropy = 0U;
		d3d_sampler_desc.ComparisonFunc = D3D11_COMPARISON_ALWAYS;
		d3d_sampler_desc.BorderColor[0] = 0.0;
		d3d_sampler_desc.BorderColor[1] = 0.0;
		d3d_sampler_desc.BorderColor[2] = 0.0;
		d3d_sampler_desc.BorderColor[3] = 1.0;
		d3d_sampler_desc.MinLOD = 0.0;
		d3d_sampler_desc.MaxLOD = 4096.0;

		HRESULT res_d3d_device_create_sampler = d3d_device->CreateSamplerState(&d3d_sampler_desc, &m_ltc_lut_sampler);
		assert(SUCCEEDED(res_d3d_device_create_sampler));
	}

	m_ltc_matrix_lut = NULL;
	{
		static_assert((4U * g_brx_ltc_lut_width * g_brx_ltc_lut_height) == (sizeof(g_brx_ltc_lut_matrices) / sizeof(g_brx_ltc_lut_matrices[0])), "");

		D3D11_TEXTURE2D_DESC d3d_texture2d_desc;
		d3d_texture2d_desc.Width = g_brx_ltc_lut_width;
		d3d_texture2d_desc.Height = g_brx_ltc_lut_height;
		d3d_texture2d_desc.MipLevels = 1U;
		d3d_texture2d_desc.ArraySize = 1U;
		d3d_texture2d_desc.Format = DXGI_FORMAT_R16G16B16A16_FLOAT;
		d3d_texture2d_desc.SampleDesc.Count = 1U;
		d3d_texture2d_desc.SampleDesc.Quality = 0U;
		d3d_texture2d_desc.Usage = D3D11_USAGE_DEFAULT;
		d3d_texture2d_desc.BindFlags = D3D11_BIND_SHADER_RESOURCE;
		d3d_texture2d_desc.CPUAccessFlags = 0U;
		d3d_texture2d_desc.MiscFlags = 0U;

		D3D11_SUBRESOURCE_DATA d3d_subresource_data;
		d3d_subresource_data.pSysMem = g_brx_ltc_lut_matrices;
		d3d_subresource_data.SysMemPitch = sizeof(uint16_t) * 4U * g_brx_ltc_lut_width;
		d3d_subresource_data.SysMemSlicePitch = sizeof(uint16_t) * 4U * g_brx_ltc_lut_width * g_brx_ltc_lut_height;

		HRESULT res_d3d_device_create_texture = d3d_device->CreateTexture2D(&d3d_texture2d_desc, &d3d_subresource_data, &m_ltc_matrix_lut);
		assert(SUCCEEDED(res_d3d_device_create_texture));
	}

	m_ltc_matrix_lut_srv = NULL;
	{
		D3D11_SHADER_RESOURCE_VIEW_DESC d3d_shader_resource_view_desc;
		d3d_shader_resource_view_desc.Format = DXGI_FORMAT_R16G16B16A16_FLOAT;
		d3d_shader_resource_view_desc.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE2DARRAY;
		d3d_shader_resource_view_desc.Texture2DArray.MostDetailedMip = 0U;
		d3d_shader_resource_view_desc.Texture2DArray.MipLevels = 1U;
		d3d_shader_resource_view_desc.Texture2DArray.FirstArraySlice = 0U;
		d3d_shader_resource_view_desc.Texture2DArray.ArraySize = 1U;

		HRESULT res_d3d_device_create_shader_resource_view = d3d_device->CreateShaderResourceView(m_ltc_matrix_lut, &d3d_shader_resource_view_desc, &m_ltc_matrix_lut_srv);
		assert(SUCCEEDED(res_d3d_device_create_shader_resource_view));
	}

	m_ltc_norm_lut = NULL;
	{
		static_assert((2u * g_brx_hdr_lut_width * g_brx_hdr_lut_height) == (sizeof(g_brx_hdr_lut_norms) / sizeof(g_brx_hdr_lut_norms[0])), "");

		D3D11_TEXTURE2D_DESC d3d_texture2d_desc;
		d3d_texture2d_desc.Width = g_brx_hdr_lut_width;
		d3d_texture2d_desc.Height = g_brx_hdr_lut_height;
		d3d_texture2d_desc.MipLevels = 1U;
		d3d_texture2d_desc.ArraySize = 1U;
		d3d_texture2d_desc.Format = DXGI_FORMAT_R16G16_FLOAT;
		d3d_texture2d_desc.SampleDesc.Count = 1U;
		d3d_texture2d_desc.SampleDesc.Quality = 0U;
		d3d_texture2d_desc.Usage = D3D11_USAGE_DEFAULT;
		d3d_texture2d_desc.BindFlags = D3D11_BIND_SHADER_RESOURCE;
		d3d_texture2d_desc.CPUAccessFlags = 0U;
		d3d_texture2d_desc.MiscFlags = 0U;

		D3D11_SUBRESOURCE_DATA d3d_subresource_data;
		d3d_subresource_data.pSysMem = g_brx_hdr_lut_norms;
		d3d_subresource_data.SysMemPitch = sizeof(uint16_t) * 2U * g_brx_hdr_lut_width;
		d3d_subresource_data.SysMemSlicePitch = sizeof(uint16_t) * 2U * g_brx_hdr_lut_width * g_brx_hdr_lut_height;

		HRESULT res_d3d_device_create_texture = d3d_device->CreateTexture2D(&d3d_texture2d_desc, &d3d_subresource_data, &m_ltc_norm_lut);
		assert(SUCCEEDED(res_d3d_device_create_texture));
	}

	m_ltc_norm_lut_srv = NULL;
	{
		D3D11_SHADER_RESOURCE_VIEW_DESC d3d_shader_resource_view_desc;
		d3d_shader_resource_view_desc.Format = DXGI_FORMAT_R16G16_FLOAT;
		d3d_shader_resource_view_desc.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE2D;
		d3d_shader_resource_view_desc.Texture2D.MostDetailedMip = 0U;
		d3d_shader_resource_view_desc.Texture2D.MipLevels = 1U;

		HRESULT res_d3d_device_create_shader_resource_view = d3d_device->CreateShaderResourceView(m_ltc_norm_lut, &d3d_shader_resource_view_desc, &m_ltc_norm_lut_srv);
		assert(SUCCEEDED(res_d3d_device_create_shader_resource_view));
	}

	m_attachment_backup_odd = NULL;
	{
		D3D11_TEXTURE2D_DESC d3d_texture2d_desc;
		d3d_texture2d_desc.Width = g_resolution_width;
		d3d_texture2d_desc.Height = g_resolution_height;
		d3d_texture2d_desc.MipLevels = 1U;
		d3d_texture2d_desc.ArraySize = 1U;
		d3d_texture2d_desc.Format = DXGI_FORMAT_R16G16B16A16_FLOAT;
		d3d_texture2d_desc.SampleDesc.Count = 1U;
		d3d_texture2d_desc.SampleDesc.Quality = 0U;
		d3d_texture2d_desc.Usage = D3D11_USAGE_DEFAULT;
		d3d_texture2d_desc.BindFlags = D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_RENDER_TARGET;
		d3d_texture2d_desc.CPUAccessFlags = 0U;
		d3d_texture2d_desc.MiscFlags = 0U;

		HRESULT res_d3d_device_create_buffer = d3d_device->CreateTexture2D(&d3d_texture2d_desc, NULL, &m_attachment_backup_odd);
		assert(SUCCEEDED(res_d3d_device_create_buffer));
	}

	m_attachment_backup_odd_rtv = NULL;
	{
		D3D11_RENDER_TARGET_VIEW_DESC d3d_render_target_view_desc;
		d3d_render_target_view_desc.Format = DXGI_FORMAT_R16G16B16A16_FLOAT;
		d3d_render_target_view_desc.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE2D;
		d3d_render_target_view_desc.Texture2D.MipSlice = 0U;

		HRESULT res_d3d_device_create_render_target_view = d3d_device->CreateRenderTargetView(m_attachment_backup_odd, &d3d_render_target_view_desc, &m_attachment_backup_odd_rtv);
		assert(SUCCEEDED(res_d3d_device_create_render_target_view));
	}

	m_attachment_backup_odd_srv = NULL;
	{
		D3D11_SHADER_RESOURCE_VIEW_DESC d3d_shader_resource_view_desc;
		d3d_shader_resource_view_desc.Format = DXGI_FORMAT_R16G16B16A16_FLOAT;
		d3d_shader_resource_view_desc.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE2D;
		d3d_shader_resource_view_desc.Texture2D.MostDetailedMip = 0U;
		d3d_shader_resource_view_desc.Texture2D.MipLevels = 1U;

		HRESULT res_d3d_device_create_shader_resource_view = d3d_device->CreateShaderResourceView(m_attachment_backup_odd, &d3d_shader_resource_view_desc, &m_attachment_backup_odd_srv);
		assert(SUCCEEDED(res_d3d_device_create_shader_resource_view));
	}

	m_attachment_depth = NULL;
	{
		D3D11_TEXTURE2D_DESC d3d_texture2d_desc;
		d3d_texture2d_desc.Width = g_resolution_width;
		d3d_texture2d_desc.Height = g_resolution_width;
		d3d_texture2d_desc.MipLevels = 1U;
		d3d_texture2d_desc.ArraySize = 1U;
		d3d_texture2d_desc.Format = DXGI_FORMAT_R32_TYPELESS;
		d3d_texture2d_desc.SampleDesc.Count = 1U;
		d3d_texture2d_desc.SampleDesc.Quality = 0U;
		d3d_texture2d_desc.Usage = D3D11_USAGE_DEFAULT;
		d3d_texture2d_desc.BindFlags = D3D11_BIND_DEPTH_STENCIL;
		d3d_texture2d_desc.CPUAccessFlags = 0U;
		d3d_texture2d_desc.MiscFlags = 0U;

		HRESULT res_d3d_device_create_buffer = d3d_device->CreateTexture2D(&d3d_texture2d_desc, NULL, &m_attachment_depth);
		assert(SUCCEEDED(res_d3d_device_create_buffer));
	}

	m_attachment_depth_dsv = NULL;
	{
		D3D11_DEPTH_STENCIL_VIEW_DESC d3d_depth_stencil_view_desc;
		d3d_depth_stencil_view_desc.Format = DXGI_FORMAT_D32_FLOAT;
		d3d_depth_stencil_view_desc.ViewDimension = D3D11_DSV_DIMENSION_TEXTURE2D;
		d3d_depth_stencil_view_desc.Flags = 0U;
		d3d_depth_stencil_view_desc.Texture2D.MipSlice = 0U;

		HRESULT res_d3d_device_create_shader_resource_view = d3d_device->CreateDepthStencilView(m_attachment_depth, &d3d_depth_stencil_view_desc, &m_attachment_depth_dsv);
		assert(SUCCEEDED(res_d3d_device_create_shader_resource_view));
	}

	g_camera_controller.m_eye_position = DirectX::XMFLOAT3(0.00000000, 6.00000000, -0.500000000);
	g_camera_controller.m_eye_direction = DirectX::XMFLOAT3(0.00000000, 0.174311504, 1.99238944);
	g_camera_controller.m_up_direction = DirectX::XMFLOAT3(0.0, 1.0, 0.0);
}

void Demo::Tick(ID3D11Device *d3d_device, ID3D11DeviceContext *d3d_device_context, IDXGISwapChain *dxgi_swap_chain)
{
	// Upload
	plane_uniform_buffer_per_frame_binding_t plane_uniform_buffer_data_per_frame_binding;
	rect_light_uniform_buffer_per_frame_binding_t rect_light_uniform_buffer_data_per_frame_binding;
	{
		// camera
		{
			DirectX::XMFLOAT3 eye_position = g_camera_controller.m_eye_position;
			DirectX::XMFLOAT3 eye_direction = g_camera_controller.m_eye_direction;
			DirectX::XMFLOAT3 up_direction = g_camera_controller.m_up_direction;

			DirectX::XMMATRIX tmp_view_transform = DirectX::XMMatrixLookToRH(DirectX::XMLoadFloat3(&eye_position), DirectX::XMLoadFloat3(&eye_direction), DirectX::XMLoadFloat3(&up_direction));
			DirectX::XMFLOAT4X4 view_transform;
			DirectX::XMStoreFloat4x4(&view_transform, tmp_view_transform);

			float fov_angle_y = 2.0 * atan((1.0 / 2.0));
			DirectX::XMMATRIX tmp_projection_transform = DirectX::XMMatrixPerspectiveFovRH(fov_angle_y, 1.0, 7.0, 7777.0);
			DirectX::XMFLOAT4X4 projection_transform;
			DirectX::XMStoreFloat4x4(&projection_transform, tmp_projection_transform);

			plane_uniform_buffer_data_per_frame_binding.view_transform = view_transform;
			plane_uniform_buffer_data_per_frame_binding.projection_transform = projection_transform;
			plane_uniform_buffer_data_per_frame_binding.eye_position = eye_position;

			rect_light_uniform_buffer_data_per_frame_binding.view_transform = view_transform;
			rect_light_uniform_buffer_data_per_frame_binding.projection_transform = projection_transform;
		}

		// light
		{
			plane_uniform_buffer_data_per_frame_binding.rect_light_vetices[0] = DirectX::XMFLOAT4(-4.0, 2.0, 32.0, 1.0);
			plane_uniform_buffer_data_per_frame_binding.rect_light_vetices[1] = DirectX::XMFLOAT4(4.0, 2.0, 32.0, 1.0);
			plane_uniform_buffer_data_per_frame_binding.rect_light_vetices[2] = DirectX::XMFLOAT4(4.0, 10.0, 32.0, 1.0);
			plane_uniform_buffer_data_per_frame_binding.rect_light_vetices[3] = DirectX::XMFLOAT4(-4.0, 10.0, 32.0, 1.0);
			plane_uniform_buffer_data_per_frame_binding.culling_range = 20.0;
			plane_uniform_buffer_data_per_frame_binding.intensity = 4.0;

			rect_light_uniform_buffer_data_per_frame_binding.rect_light_vetices[0] = DirectX::XMFLOAT4(-4.0, 2.0, 32.0, 1.0);
			rect_light_uniform_buffer_data_per_frame_binding.rect_light_vetices[1] = DirectX::XMFLOAT4(4.0, 2.0, 32.0, 1.0);
			rect_light_uniform_buffer_data_per_frame_binding.rect_light_vetices[2] = DirectX::XMFLOAT4(-4.0, 10.0, 32.0, 1.0);
			rect_light_uniform_buffer_data_per_frame_binding.rect_light_vetices[3] = DirectX::XMFLOAT4(4.0, 10.0, 32.0, 1.0);
			rect_light_uniform_buffer_data_per_frame_binding.intensity = 4.0;
		}

		// mesh
		{
			DirectX::XMFLOAT4X4 model_transform;
			DirectX::XMStoreFloat4x4(&model_transform, DirectX::XMMatrixIdentity());
			plane_uniform_buffer_data_per_frame_binding.model_transform = model_transform;
			plane_uniform_buffer_data_per_frame_binding.dcolor = DirectX::XMFLOAT3(1.0, 1.0, 1.0);
			plane_uniform_buffer_data_per_frame_binding.scolor = DirectX::XMFLOAT3(0.23, 0.23, 0.23);
			plane_uniform_buffer_data_per_frame_binding.roughness = 0.25;
		}
	}
	d3d_device_context->UpdateSubresource(m_plane_uniform_buffer_per_frame_binding, 0U, NULL, &plane_uniform_buffer_data_per_frame_binding, sizeof(plane_uniform_buffer_per_frame_binding_t), sizeof(plane_uniform_buffer_per_frame_binding_t));
	d3d_device_context->UpdateSubresource(m_rect_light_uniform_buffer_per_frame_binding, 0U, NULL, &rect_light_uniform_buffer_data_per_frame_binding, sizeof(rect_light_uniform_buffer_per_frame_binding_t), sizeof(rect_light_uniform_buffer_per_frame_binding_t));

	// Light Pass
	{
		d3d_device_context->OMSetRenderTargets(1U, &m_attachment_backup_odd_rtv, m_attachment_depth_dsv);

		FLOAT color[4] = {0.0f, 0.0f, 0.0f, 0.0f};
		d3d_device_context->ClearRenderTargetView(m_attachment_backup_odd_rtv, color);

		FLOAT depth = 1.0f;
		d3d_device_context->ClearDepthStencilView(m_attachment_depth_dsv, D3D11_CLEAR_DEPTH, depth, 0U);

		// glEnable(GL_DEPTH_TEST);
		// glDepthFunc(GL_LESS);
		// glColorMaski(0, GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

		D3D11_VIEWPORT viewport = {0.0f, 0.0f, g_resolution_width, g_resolution_height, 0.0f, 1.0f};
		d3d_device_context->RSSetViewports(1U, &viewport);

		// Draw Plane
		{
			d3d_device_context->RSSetState(m_plane_rs);

			d3d_device_context->VSSetShader(m_plane_vs, NULL, 0U);
			d3d_device_context->PSSetShader(m_plane_fs, NULL, 0U);

			d3d_device_context->VSSetConstantBuffers(0U, 1U, &m_plane_uniform_buffer_per_frame_binding);
			d3d_device_context->PSSetConstantBuffers(0U, 1U, &m_plane_uniform_buffer_per_frame_binding);

			d3d_device_context->PSSetSamplers(0U, 1U, &m_ltc_lut_sampler);
			d3d_device_context->PSSetShaderResources(0U, 1U, &m_ltc_matrix_lut_srv);
			d3d_device_context->PSSetShaderResources(1U, 1U, &m_ltc_norm_lut_srv);

			d3d_device_context->IASetInputLayout(m_plane_vao);

			ID3D11Buffer *vertex_buffers[2] = {m_plane_vb_position, m_plane_vb_varying};
			UINT strides[2] = {sizeof(float) * 3, sizeof(float) * 3};
			UINT offsets[2] = {0U, 0U};
			d3d_device_context->IASetVertexBuffers(0U, 2U, vertex_buffers, strides, offsets);
			d3d_device_context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLESTRIP);
			d3d_device_context->DrawInstanced(4U, 1U, 0U, 0U);
		}

		// Draw Rect Light
		{
			d3d_device_context->RSSetState(m_rect_light_rs);

			d3d_device_context->VSSetShader(m_rect_light_vs, NULL, 0U);
			d3d_device_context->PSSetShader(m_rect_light_fs, NULL, 0U);

			d3d_device_context->VSSetConstantBuffers(0U, 1U, &m_rect_light_uniform_buffer_per_frame_binding);
			d3d_device_context->PSSetConstantBuffers(0U, 1U, &m_rect_light_uniform_buffer_per_frame_binding);

			d3d_device_context->IASetInputLayout(m_rect_light_vao);
			d3d_device_context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLESTRIP);
			d3d_device_context->DrawInstanced(4U, 1U, 0U, 0U);
		}
	}

	// Post Process Pass
	{
		d3d_device_context->OMSetRenderTargets(1U, &m_attachment_backbuffer_rtv, NULL);

		d3d_device_context->RSSetState(m_post_process_rs);

		// glDisable(GL_DEPTH_TEST);

		d3d_device_context->VSSetShader(m_post_process_vs, NULL, 0U);
		d3d_device_context->PSSetShader(m_post_process_fs, NULL, 0U);

		d3d_device_context->PSSetSamplers(0U, 1U, &m_ltc_lut_sampler);
		d3d_device_context->PSSetShaderResources(0U, 1U, &m_attachment_backup_odd_srv);

		d3d_device_context->IASetInputLayout(m_post_process_vao);
		d3d_device_context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLESTRIP);
		d3d_device_context->DrawInstanced(3U, 1U, 0U, 0U);

		// unbind
		ID3D11ShaderResourceView *shader_resource_views = {NULL};
		d3d_device_context->PSSetShaderResources(0U, 1U, &shader_resource_views);
	}

	HRESULT res_dxgi_swap_chain_present = dxgi_swap_chain->Present(1U, 0U);
	assert(SUCCEEDED(res_dxgi_swap_chain_present));
}
