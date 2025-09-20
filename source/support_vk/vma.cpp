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

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX 1
#define VK_USE_PLATFORM_WIN32_KHR 1
#include <sdkddkver.h>
#include <windows.h>
#include "../../thirdparty/Vulkan-Headers/include/vulkan/vulkan.h"
#include "../../thirdparty/VulkanMemoryAllocator/include/vk_mem_alloc.h"

#define VMA_IMPLEMENTATION 1
#define VMA_STATIC_VULKAN_FUNCTIONS 0
#define VMA_DYNAMIC_VULKAN_FUNCTIONS 1
#include "../../thirdparty/VulkanMemoryAllocator/include/vk_mem_alloc.h"
