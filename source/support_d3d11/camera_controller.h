#ifndef _CAMERA_CONTROLLER_H_
#define _CAMERA_CONTROLLER_H_ 1

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN 1
#endif
#ifndef NOMINMAX
#define NOMINMAX 1
#endif
#include <sdkddkver.h>
#include <windows.h>
#include <DirectXMath.h>

class CameraController
{
	float Previous_X;
	float Previous_Y;

public:
	DirectX::XMFLOAT3 m_eye_position;
	DirectX::XMFLOAT3 m_eye_direction;
	DirectX::XMFLOAT3 m_up_direction;

	void MoveForward();
	void MoveBack();
	void MoveLeft();
	void MoveRight();
	void MoveUp();
	void MoveDown();

	void OnMouseMove(float x, float y, bool hold);
};

extern class CameraController g_camera_controller;

#endif