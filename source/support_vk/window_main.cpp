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

#define NOMINMAX 1
#define WIN32_LEAN_AND_MEAN 1
#include <sdkddkver.h>
#include <Windows.h>
#include <windowsx.h>
#include <process.h>
#include "camera_controller.h"
#include "resolution.h"
#include "render_main.h"
#include "window_main.h"

bool volatile g_window_quit = false;

static HANDLE render_main_thread = NULL;

static LRESULT CALLBACK wnd_proc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);

int APIENTRY wWinMain(_In_ HINSTANCE hInstance, _In_opt_ HINSTANCE, _In_ LPWSTR lpCmdLine, _In_ int nCmdShow)
{
	ATOM hWndCls;
	{
		WNDCLASSEXW Desc = {
			sizeof(WNDCLASSEX),
			CS_OWNDC,
			wnd_proc,
			0,
			0,
			hInstance,
			LoadIconW(NULL, IDI_APPLICATION),
			LoadCursorW(NULL, IDC_ARROW),
			(HBRUSH)(COLOR_WINDOW + 1),
			NULL,
			L"Demo:0XFFFFFFFF",
			LoadIconW(NULL, IDI_APPLICATION),
		};
		hWndCls = RegisterClassExW(&Desc);
	}

	HWND hWnd;
	{
		HWND hDesktop = GetDesktopWindow();
		HMONITOR hMonitor = MonitorFromWindow(hDesktop, MONITOR_DEFAULTTONEAREST);
		MONITORINFOEXW MonitorInfo;
		MonitorInfo.cbSize = sizeof(MONITORINFOEXW);
		GetMonitorInfoW(hMonitor, &MonitorInfo);

		RECT rect = {(MonitorInfo.rcWork.left + MonitorInfo.rcWork.right) / 2 - g_resolution_width / 2,
					 (MonitorInfo.rcWork.bottom + MonitorInfo.rcWork.top) / 2 - g_resolution_height / 2,
					 (MonitorInfo.rcWork.left + MonitorInfo.rcWork.right) / 2 + g_resolution_width / 2,
					 (MonitorInfo.rcWork.bottom + MonitorInfo.rcWork.top) / 2 + g_resolution_height / 2};
		AdjustWindowRectEx(&rect, WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU | WS_MINIMIZEBOX | WS_MAXIMIZEBOX | WS_THICKFRAME | WS_VISIBLE | WS_CLIPCHILDREN | WS_CLIPSIBLINGS, FALSE, WS_EX_APPWINDOW);

		hWnd = CreateWindowExW(WS_EX_APPWINDOW,
							   MAKEINTATOM(hWndCls),
							   L"Demo",
							   WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU | WS_MINIMIZEBOX | WS_MAXIMIZEBOX | WS_THICKFRAME | WS_VISIBLE | WS_CLIPCHILDREN | WS_CLIPSIBLINGS,
							   rect.left,
							   rect.top,
							   rect.right - rect.left,
							   rect.bottom - rect.top,
							   hDesktop,
							   NULL,
							   hInstance,
							   NULL);
	}

	render_main_thread = reinterpret_cast<HANDLE>(_beginthreadex(NULL, 0U, render_main, hWnd, 0U, NULL));

	MSG msg;
	while (GetMessageW(&msg, NULL, 0, 0))
	{
		TranslateMessage(&msg);
		DispatchMessageW(&msg);
	}

	return (int)msg.wParam;
}

static LRESULT CALLBACK wnd_proc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	static float g_SpeedTranslate = 0.5f;

	switch (message)
	{
	case WM_DESTROY:
	{
		g_window_quit = true;

		DWORD res_wait_for_single_object = WaitForSingleObjectEx(render_main_thread, INFINITE, FALSE);
		assert(WAIT_OBJECT_0 == res_wait_for_single_object);
		CloseHandle(render_main_thread);

		PostQuitMessage(0);
	}
		return 0;
	case WM_ERASEBKGND:
		return 1;
	case WM_KEYDOWN:
	{
		switch (wParam)
		{
		case 'W':
		{
			g_camera_controller.MoveForward();
		}
		break;
		case 'S':
		{
			g_camera_controller.MoveBack();
		}
		break;
		case 'A':
		{
			g_camera_controller.MoveLeft();
		}
		break;
		case 'D':
		{
			g_camera_controller.MoveRight();
		}
		break;
		case 'Q':
		{
			g_camera_controller.MoveDown();
		}
		break;
		case 'E':
		{
			g_camera_controller.MoveUp();
		}
		break;
		}
	}
		return 0;
	case WM_MOUSEMOVE:
	{
		int Current_X = GET_X_LPARAM(lParam);
		int Current_Y = GET_Y_LPARAM(lParam);

		float CurrentNormalized_X = static_cast<float>(Current_X) / g_resolution_width;
		float CurrentNormalized_Y = static_cast<float>(Current_Y) / g_resolution_height;

		g_camera_controller.OnMouseMove(CurrentNormalized_X, CurrentNormalized_Y, (0 != (wParam & MK_RBUTTON)));
	}
		return 0;
	default:
		return DefWindowProc(hWnd, message, wParam, lParam);
	}
}