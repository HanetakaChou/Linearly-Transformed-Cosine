#define NOMINMAX 1
#define WIN32_LEAN_AND_MEAN 1
#include <sdkddkver.h>
#include <Windows.h>
#include <cmath>
#include <algorithm>
#include "camera_controller.h"

class CameraController g_camera_controller;

static float g_SpeedTranslate = 0.5f;

void CameraController::MoveForward()
{
	m_eye_position.x += m_eye_direction.x * g_SpeedTranslate;
	m_eye_position.y += m_eye_direction.y * g_SpeedTranslate;
	m_eye_position.z += m_eye_direction.z * g_SpeedTranslate;
}

void CameraController::MoveBack()
{
	m_eye_position.x -= m_eye_direction.x * g_SpeedTranslate;
	m_eye_position.y -= m_eye_direction.y * g_SpeedTranslate;
	m_eye_position.z -= m_eye_direction.z * g_SpeedTranslate;
}

void CameraController::MoveLeft()
{
	DirectX::XMFLOAT3 EyeDirection = {m_eye_direction.x, m_eye_direction.y, m_eye_direction.z};
	DirectX::XMFLOAT3 UpDirection = {m_up_direction.x, m_up_direction.y, m_up_direction.z};

	DirectX::XMVECTOR NegEyeDirection = DirectX::XMVectorNegate(DirectX::XMLoadFloat3(&EyeDirection));
	DirectX::XMVECTOR AxisForwardDirection = DirectX::XMVector3Normalize(NegEyeDirection);
	DirectX::XMVECTOR AxisRightDirection = DirectX::XMVector3Normalize(DirectX::XMVector3Cross(DirectX::XMLoadFloat3(&UpDirection), AxisForwardDirection));
	// DirectX::XMVECTOR AxisViewUpDirection = DirectX::XMVector3Cross(AxisForwardDirection, AxisRightDirection);

	DirectX::XMFLOAT3 AxisRightDirectionF;
	DirectX::XMStoreFloat3(&AxisRightDirectionF, AxisRightDirection);

	m_eye_position.x -= AxisRightDirectionF.x * g_SpeedTranslate;
	m_eye_position.y -= AxisRightDirectionF.y * g_SpeedTranslate;
	m_eye_position.z -= AxisRightDirectionF.z * g_SpeedTranslate;
}

void CameraController::MoveRight()
{
	DirectX::XMFLOAT3 EyeDirection = {m_eye_direction.x, m_eye_direction.y, m_eye_direction.z};
	DirectX::XMFLOAT3 UpDirection = {m_up_direction.x, m_up_direction.y, m_up_direction.z};

	DirectX::XMVECTOR NegEyeDirection = DirectX::XMVectorNegate(DirectX::XMLoadFloat3(&EyeDirection));
	DirectX::XMVECTOR AxisForwardDirection = DirectX::XMVector3Normalize(NegEyeDirection);
	DirectX::XMVECTOR AxisRightDirection = DirectX::XMVector3Normalize(DirectX::XMVector3Cross(DirectX::XMLoadFloat3(&UpDirection), AxisForwardDirection));
	// DirectX::XMVECTOR AxisViewUpDirection = DirectX::XMVector3Cross(AxisForwardDirection, AxisRightDirection);

	DirectX::XMFLOAT3 AxisRightDirectionF;
	DirectX::XMStoreFloat3(&AxisRightDirectionF, AxisRightDirection);

	m_eye_position.x += AxisRightDirectionF.x * g_SpeedTranslate;
	m_eye_position.y += AxisRightDirectionF.y * g_SpeedTranslate;
	m_eye_position.z += AxisRightDirectionF.z * g_SpeedTranslate;
}

void CameraController::MoveUp()
{
	DirectX::XMFLOAT3 EyeDirection = {m_eye_direction.x, m_eye_direction.y, m_eye_direction.z};
	DirectX::XMFLOAT3 UpDirection = {m_up_direction.x, m_up_direction.y, m_up_direction.z};

	DirectX::XMVECTOR NegEyeDirection = DirectX::XMVectorNegate(DirectX::XMLoadFloat3(&EyeDirection));
	DirectX::XMVECTOR AxisForwardDirection = DirectX::XMVector3Normalize(NegEyeDirection);
	DirectX::XMVECTOR AxisRightDirection = DirectX::XMVector3Normalize(DirectX::XMVector3Cross(DirectX::XMLoadFloat3(&UpDirection), AxisForwardDirection));
	DirectX::XMVECTOR AxisViewUpDirection = DirectX::XMVector3Cross(AxisForwardDirection, AxisRightDirection);

	DirectX::XMFLOAT3 AxisViewUpDirectionF;
	DirectX::XMStoreFloat3(&AxisViewUpDirectionF, AxisViewUpDirection);

	m_eye_position.x += AxisViewUpDirectionF.x * g_SpeedTranslate;
	m_eye_position.y += AxisViewUpDirectionF.y * g_SpeedTranslate;
	m_eye_position.z += AxisViewUpDirectionF.z * g_SpeedTranslate;
}

void CameraController::MoveDown()
{
	DirectX::XMFLOAT3 EyeDirection = {m_eye_direction.x, m_eye_direction.y, m_eye_direction.z};
	DirectX::XMFLOAT3 UpDirection = {m_up_direction.x, m_up_direction.y, m_up_direction.z};

	DirectX::XMVECTOR NegEyeDirection = DirectX::XMVectorNegate(DirectX::XMLoadFloat3(&EyeDirection));
	DirectX::XMVECTOR AxisForwardDirection = DirectX::XMVector3Normalize(NegEyeDirection);
	DirectX::XMVECTOR AxisRightDirection = DirectX::XMVector3Normalize(DirectX::XMVector3Cross(DirectX::XMLoadFloat3(&UpDirection), AxisForwardDirection));
	DirectX::XMVECTOR AxisViewUpDirection = DirectX::XMVector3Cross(AxisForwardDirection, AxisRightDirection);

	DirectX::XMFLOAT3 AxisViewUpDirectionF;
	DirectX::XMStoreFloat3(&AxisViewUpDirectionF, AxisViewUpDirection);

	m_eye_position.x -= AxisViewUpDirectionF.x * g_SpeedTranslate;
	m_eye_position.y -= AxisViewUpDirectionF.y * g_SpeedTranslate;
	m_eye_position.z -= AxisViewUpDirectionF.z * g_SpeedTranslate;
}

void CameraController::OnMouseMove(float x, float y, bool hold)
{
	float Current_X = std::min(std::max(0.0f, x), 1.0f);
	float Current_Y = std::min(std::max(0.0f, y), 1.0f);

	float Offset_X = Current_X - Previous_X;
	float Offset_Y = Current_Y - Previous_Y;

	if (hold)
	{
		float LengthNormalized = ::sqrtf(Offset_X * Offset_X + Offset_Y * Offset_Y);

		if ((Offset_Y < Offset_X) && (Offset_Y > -Offset_X))
		{
			// right
			DirectX::XMFLOAT3 EyeDirection = {m_eye_direction.x, m_eye_direction.y, m_eye_direction.z};
			DirectX::XMFLOAT3 UpDirection = {m_up_direction.x, m_up_direction.y, m_up_direction.z};

			DirectX::XMVECTOR NegEyeDirection = DirectX::XMVectorNegate(DirectX::XMLoadFloat3(&EyeDirection));
			DirectX::XMVECTOR AxisForwardDirection = DirectX::XMVector3Normalize(NegEyeDirection);
			DirectX::XMVECTOR AxisRightDirection = DirectX::XMVector3Normalize(DirectX::XMVector3Cross(DirectX::XMLoadFloat3(&UpDirection), AxisForwardDirection));
			DirectX::XMVECTOR AxisViewUpDirection = DirectX::XMVector3Cross(AxisForwardDirection, AxisRightDirection);

			DirectX::XMFLOAT3 EyeDirectionNew;
			DirectX::XMStoreFloat3(&EyeDirectionNew, DirectX::XMVectorNegate(DirectX::XMVector3Transform(AxisForwardDirection, DirectX::XMMatrixRotationAxis(AxisViewUpDirection, -2.0f * LengthNormalized))));
			m_eye_direction.x = EyeDirectionNew.x;
			m_eye_direction.y = EyeDirectionNew.y;
			m_eye_direction.z = EyeDirectionNew.z;
		}
		else if ((Offset_Y < -Offset_X) && (Offset_Y > Offset_X))
		{
			// left
			DirectX::XMFLOAT3 EyeDirection = {m_eye_direction.x, m_eye_direction.y, m_eye_direction.z};
			DirectX::XMFLOAT3 UpDirection = {m_up_direction.x, m_up_direction.y, m_up_direction.z};

			DirectX::XMVECTOR NegEyeDirection = DirectX::XMVectorNegate(DirectX::XMLoadFloat3(&EyeDirection));
			DirectX::XMVECTOR AxisForwardDirection = DirectX::XMVector3Normalize(NegEyeDirection);
			DirectX::XMVECTOR AxisRightDirection = DirectX::XMVector3Normalize(DirectX::XMVector3Cross(DirectX::XMLoadFloat3(&UpDirection), AxisForwardDirection));
			DirectX::XMVECTOR AxisViewUpDirection = DirectX::XMVector3Cross(AxisForwardDirection, AxisRightDirection);

			DirectX::XMFLOAT3 EyeDirectionNew;
			DirectX::XMStoreFloat3(&EyeDirectionNew, DirectX::XMVectorNegate(DirectX::XMVector3Transform(AxisForwardDirection, DirectX::XMMatrixRotationAxis(AxisViewUpDirection, 2.0f * LengthNormalized))));
			m_eye_direction.x = EyeDirectionNew.x;
			m_eye_direction.y = EyeDirectionNew.y;
			m_eye_direction.z = EyeDirectionNew.z;
		}
		else if ((Offset_Y < Offset_X) && (Offset_Y < -Offset_X))
		{
			// up
			DirectX::XMFLOAT3 EyeDirection = {m_eye_direction.x, m_eye_direction.y, m_eye_direction.z};
			DirectX::XMFLOAT3 UpDirection = {m_up_direction.x, m_up_direction.y, m_up_direction.z};

			DirectX::XMVECTOR NegEyeDirection = DirectX::XMVectorNegate(DirectX::XMLoadFloat3(&EyeDirection));
			DirectX::XMVECTOR AxisForwardDirection = DirectX::XMVector3Normalize(NegEyeDirection);
			DirectX::XMVECTOR AxisRightDirection = DirectX::XMVector3Normalize(DirectX::XMVector3Cross(DirectX::XMLoadFloat3(&UpDirection), AxisForwardDirection));
			DirectX::XMVECTOR AxisViewUpDirection = DirectX::XMVector3Cross(AxisForwardDirection, AxisRightDirection);

			DirectX::XMFLOAT3 EyeDirectionNew;
			DirectX::XMStoreFloat3(&EyeDirectionNew, DirectX::XMVectorNegate(DirectX::XMVector3Transform(AxisForwardDirection, DirectX::XMMatrixRotationAxis(AxisRightDirection, 2.0f * LengthNormalized))));
			m_eye_direction.x = EyeDirectionNew.x;
			m_eye_direction.y = EyeDirectionNew.y;
			m_eye_direction.z = EyeDirectionNew.z;
		}
		else if ((Offset_Y > Offset_X) && (Offset_Y > -Offset_X))
		{
			// down
			DirectX::XMFLOAT3 EyeDirection = {m_eye_direction.x, m_eye_direction.y, m_eye_direction.z};
			DirectX::XMFLOAT3 UpDirection = {m_up_direction.x, m_up_direction.y, m_up_direction.z};

			DirectX::XMVECTOR NegEyeDirection = DirectX::XMVectorNegate(DirectX::XMLoadFloat3(&EyeDirection));
			DirectX::XMVECTOR AxisForwardDirection = DirectX::XMVector3Normalize(NegEyeDirection);
			DirectX::XMVECTOR AxisRightDirection = DirectX::XMVector3Normalize(DirectX::XMVector3Cross(DirectX::XMLoadFloat3(&UpDirection), AxisForwardDirection));
			DirectX::XMVECTOR AxisViewUpDirection = DirectX::XMVector3Cross(AxisForwardDirection, AxisRightDirection);

			DirectX::XMFLOAT3 EyeDirectionNew;
			DirectX::XMStoreFloat3(&EyeDirectionNew, DirectX::XMVectorNegate(DirectX::XMVector3Transform(AxisForwardDirection, DirectX::XMMatrixRotationAxis(AxisRightDirection, -2.0f * LengthNormalized))));
			m_eye_direction.x = EyeDirectionNew.x;
			m_eye_direction.y = EyeDirectionNew.y;
			m_eye_direction.z = EyeDirectionNew.z;
		}
	}

	Previous_X = Current_X;
	Previous_Y = Current_Y;
}
