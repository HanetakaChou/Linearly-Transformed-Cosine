#ifndef _DEMO_D3D11_H_
#define _DEMO_D3D11_H_ 1

#include <sdkddkver.h>
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN 1
#endif
#ifndef NOMINMAX
#define NOMINMAX 1
#endif
#include <windows.h>

#include <dxgi.h>
#include <d3d11.h>

class Demo
{
	ID3D11RenderTargetView *m_attachment_backbuffer_rtv;

	ID3D11Buffer *m_plane_vb_position;
	ID3D11Buffer *m_plane_vb_varying;
	ID3D11InputLayout *m_plane_vao;
	ID3D11VertexShader *m_plane_vs;
	ID3D11PixelShader *m_plane_fs;
	ID3D11RasterizerState *m_plane_rs;
	ID3D11Buffer *m_plane_uniform_buffer_per_frame_binding;

	ID3D11InputLayout *m_rect_light_vao;
	ID3D11VertexShader *m_rect_light_vs;
	ID3D11PixelShader *m_rect_light_fs;
	ID3D11RasterizerState *m_rect_light_rs;
	ID3D11Buffer *m_rect_light_uniform_buffer_per_frame_binding;

	ID3D11InputLayout *m_post_process_vao;
	ID3D11VertexShader *m_post_process_vs;
	ID3D11PixelShader *m_post_process_fs;
	ID3D11RasterizerState *m_post_process_rs;

	ID3D11SamplerState *m_ltc_lut_sampler;
	ID3D11Texture2D *m_ltc_matrix_lut;
	ID3D11ShaderResourceView *m_ltc_matrix_lut_srv;
	ID3D11Texture2D *m_ltc_norm_lut;
	ID3D11ShaderResourceView *m_ltc_norm_lut_srv;

	ID3D11Texture2D *m_attachment_backup_odd;
	ID3D11RenderTargetView *m_attachment_backup_odd_rtv;
	ID3D11ShaderResourceView *m_attachment_backup_odd_srv;
	ID3D11Texture2D *m_attachment_depth;
	ID3D11DepthStencilView *m_attachment_depth_dsv;

public:
	void Init(ID3D11Device *d3d11_device, ID3D11DeviceContext *d3d11_device_context, IDXGISwapChain *dxgi_swap_chain);
	void Tick(ID3D11Device *d3d11_device, ID3D11DeviceContext *d3d11_device_context, IDXGISwapChain *dxgi_swap_chain);
};

#endif