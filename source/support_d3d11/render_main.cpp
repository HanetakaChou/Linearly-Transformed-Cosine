#define WIN32_LEAN_AND_MEAN
#define NOMINMAX 1
#include <sdkddkver.h>
#include <windows.h>

#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include <dxgi.h>
#include <d3d11.h>

#include "resolution.h"
#include "window_main.h"
#include "render_main.h"

#include "../demo_d3d11.h"

unsigned __stdcall render_main(void *pVoid)
{
	HWND hWnd = static_cast<HWND>(pVoid);

	ID3D11Device *d3d11_device = NULL;
	ID3D11DeviceContext *d3d11_device_context = NULL;
	IDXGISwapChain *dxgi_swap_chain = NULL;
	{
		IDXGIFactory1 *dxgi_factory = NULL;
		HRESULT res_create_dxgi_factory = CreateDXGIFactory1(IID_PPV_ARGS(&dxgi_factory));
		assert(SUCCEEDED(res_create_dxgi_factory));

		IDXGIAdapter1 *dxgi_adapter = NULL;
		// "0U" means "default" adapter // perhaps user-select
		HRESULT res_dxgi_factory_enum_adapters = dxgi_factory->EnumAdapters1(0U, &dxgi_adapter);
		assert(SUCCEEDED(res_dxgi_factory_enum_adapters));

		D3D_FEATURE_LEVEL pFeatureLevels[1] = {D3D_FEATURE_LEVEL_11_0};
		HRESULT res_d3d11_create_device = D3D11CreateDevice(
			dxgi_adapter,
			D3D_DRIVER_TYPE_UNKNOWN,
			NULL,
#ifndef NDEBUG
			D3D11_CREATE_DEVICE_SINGLETHREADED | D3D11_CREATE_DEVICE_DEBUG,
#else
			D3D11_CREATE_DEVICE_SINGLETHREADED,
#endif
			pFeatureLevels,
			sizeof(pFeatureLevels) / sizeof(pFeatureLevels[0]),
			D3D11_SDK_VERSION,
			&d3d11_device,
			NULL,
			&d3d11_device_context);
		assert(SUCCEEDED(res_d3d11_create_device));

		DXGI_SWAP_CHAIN_DESC dxgi_swap_chain_desc;
		dxgi_swap_chain_desc.BufferDesc.Width = g_resolution_width;
		dxgi_swap_chain_desc.BufferDesc.Height = g_resolution_height;
		dxgi_swap_chain_desc.BufferDesc.RefreshRate.Numerator = 60U;
		dxgi_swap_chain_desc.BufferDesc.RefreshRate.Denominator = 1U;
		dxgi_swap_chain_desc.BufferDesc.Format = DXGI_FORMAT_B8G8R8A8_UNORM;
		dxgi_swap_chain_desc.BufferDesc.ScanlineOrdering = DXGI_MODE_SCANLINE_ORDER_UNSPECIFIED;
		dxgi_swap_chain_desc.BufferDesc.Scaling = DXGI_MODE_SCALING_UNSPECIFIED;
		dxgi_swap_chain_desc.SampleDesc.Count = 1U;
		dxgi_swap_chain_desc.SampleDesc.Quality = 0U;
		dxgi_swap_chain_desc.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
		dxgi_swap_chain_desc.BufferCount = 1U;
		dxgi_swap_chain_desc.OutputWindow = hWnd;
		dxgi_swap_chain_desc.Windowed = TRUE;
		dxgi_swap_chain_desc.SwapEffect = DXGI_SWAP_EFFECT_DISCARD;
		dxgi_swap_chain_desc.Flags = 0U;
		HRESULT res_dxgi_factory_create_swap_chain = dxgi_factory->CreateSwapChain(d3d11_device, &dxgi_swap_chain_desc, &dxgi_swap_chain);
		assert(SUCCEEDED(res_dxgi_factory_create_swap_chain));

		dxgi_adapter->Release();
		dxgi_factory->Release();
	}

	class Demo demo;

	demo.Init(d3d11_device, d3d11_device_context, dxgi_swap_chain);

	while (!g_window_quit)
	{
		demo.Tick(d3d11_device, d3d11_device_context, dxgi_swap_chain);
	}

	d3d11_device_context->Release();
	d3d11_device->Release();

	return 0U;
}
