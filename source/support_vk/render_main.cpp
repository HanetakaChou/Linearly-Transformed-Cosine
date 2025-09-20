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
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include "../../thirdparty/Vulkan-Headers/include/vulkan/vulkan.h"
#include "../../thirdparty/VulkanMemoryAllocator/include/vk_mem_alloc.h"
#include "window_main.h"
#include "render_main.h"
#include "streaming.h"
#include "frame_throttling.h"
#include "resolution.h"
#include "../demo_vk.h"

static VkBool32 VKAPI_CALL vulkan_debug_utils_messenger_callback(VkDebugUtilsMessageSeverityFlagBitsEXT messageSeverity, VkDebugUtilsMessageTypeFlagsEXT messageTypes, const VkDebugUtilsMessengerCallbackDataEXT *pCallbackData, void *pUserData);

static inline uint32_t __internal_find_lowest_memory_type_index(struct VkPhysicalDeviceMemoryProperties const *physical_device_memory_properties, VkDeviceSize memory_requirements_size, uint32_t memory_requirements_memory_type_bits, VkMemoryPropertyFlags required_property_flags);

static inline uint32_t __internal_find_lowest_memory_type_index(struct VkPhysicalDeviceMemoryProperties const *physical_device_memory_properties, VkDeviceSize memory_requirements_size, uint32_t memory_requirements_memory_type_bits, VkMemoryPropertyFlags required_property_flags, VkMemoryPropertyFlags preferred_property_flags);

static inline void vulkan_query_frame_buffer_extent(
	VkInstance vulkan_instance, PFN_vkGetInstanceProcAddr pfn_get_instance_proc_addr, VkPhysicalDevice vulkan_physical_device,
	VkSurfaceKHR vulkan_surface,
	uint32_t &vulkan_framebuffer_width, uint32_t &vulkan_framebuffer_height);

static inline void vulkan_create_swapchain(
	VkDevice vulkan_device, PFN_vkGetDeviceProcAddr pfn_get_device_proc_addr, VkAllocationCallbacks *vulkan_allocation_callbacks,
	VkSurfaceKHR vulkan_surface, VkFormat vulkan_swapchain_image_format,
	uint32_t vulkan_framebuffer_width, uint32_t vulkan_framebuffer_height,
	VkSwapchainKHR &vulkan_swapchain,
	uint32_t &vulkan_swapchain_image_count,
	std::vector<VkImageView> &vulkan_swapchain_image_views);

static inline void vulkan_destroy_swapchain(
	VkDevice vulkan_device, PFN_vkGetDeviceProcAddr pfn_get_device_proc_addr, VkAllocationCallbacks *vulkan_allocation_callbacks,
	uint32_t &vulkan_framebuffer_width, uint32_t &vulkan_framebuffer_height,
	VkSwapchainKHR &vulkan_swapchain,
	uint32_t &vulkan_swapchain_image_count,
	std::vector<VkImageView> &vulkan_swapchain_image_views);

unsigned __stdcall render_main(void *pVoid)
{
	HWND hWnd = static_cast<HWND>(pVoid);

	VkAllocationCallbacks *vulkan_allocation_callbacks = NULL;

#ifndef NDEBUG
	{
		WCHAR file_name[4096];
		DWORD res_get_file_name = GetModuleFileNameW(NULL, file_name, sizeof(file_name) / sizeof(file_name[0]));
		assert(0U != res_get_file_name);

		for (int i = (res_get_file_name - 1); i > 0; --i)
		{
			if (L'\\' == file_name[i])
			{
				file_name[i] = L'\0';
				break;
			}
		}

		BOOL res_set_environment_variable = SetEnvironmentVariableW(L"VK_LAYER_PATH", file_name);
		assert(FALSE != res_set_environment_variable);
	}
#endif

	PFN_vkGetInstanceProcAddr pfn_get_instance_proc_addr = vkGetInstanceProcAddr;

	uint32_t vulkan_api_version = VK_API_VERSION_1_0;

	VkInstance vulkan_instance = VK_NULL_HANDLE;
	{
		PFN_vkCreateInstance pfn_vk_create_instance = reinterpret_cast<PFN_vkCreateInstance>(pfn_get_instance_proc_addr(VK_NULL_HANDLE, "vkCreateInstance"));
		assert(NULL != pfn_vk_create_instance);

		VkApplicationInfo application_info = {
			VK_STRUCTURE_TYPE_APPLICATION_INFO,
			NULL,
			"Windows-Vulkan-Demo",
			0,
			"Windows-Vulkan-Demo",
			0,
			VK_API_VERSION_1_0};

#ifndef NDEBUG
		char const *enabled_layer_names[] = {
			"VK_LAYER_KHRONOS_validation"};
		char const *enabled_extension_names[] = {
			VK_KHR_SURFACE_EXTENSION_NAME,
			VK_KHR_WIN32_SURFACE_EXTENSION_NAME,
			VK_EXT_DEBUG_UTILS_EXTENSION_NAME};
#else
		char const *enabled_extension_names[] = {
			VK_KHR_SURFACE_EXTENSION_NAME,
			VK_KHR_WIN32_SURFACE_EXTENSION_NAME};
#endif

		VkInstanceCreateInfo instance_create_info = {
			VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO,
			NULL,
			0U,
			&application_info,
#ifndef NDEBUG
			sizeof(enabled_layer_names) / sizeof(enabled_layer_names[0]),
			enabled_layer_names,
#else
			0U,
			NULL,
#endif
			sizeof(enabled_extension_names) / sizeof(enabled_extension_names[0]),
			enabled_extension_names};

		VkResult res_create_instance = pfn_vk_create_instance(&instance_create_info, vulkan_allocation_callbacks, &vulkan_instance);
		assert(VK_SUCCESS == res_create_instance);
	}

	pfn_get_instance_proc_addr = reinterpret_cast<PFN_vkGetInstanceProcAddr>(pfn_get_instance_proc_addr(vulkan_instance, "vkGetInstanceProcAddr"));
	assert(NULL != pfn_get_instance_proc_addr);

#ifndef NDEBUG
	VkDebugUtilsMessengerEXT vulkan_messenge;
	{
		PFN_vkCreateDebugUtilsMessengerEXT pfn_vk_create_debug_utils_messenger = reinterpret_cast<PFN_vkCreateDebugUtilsMessengerEXT>(pfn_get_instance_proc_addr(vulkan_instance, "vkCreateDebugUtilsMessengerEXT"));
		assert(NULL != pfn_vk_create_debug_utils_messenger);

		VkDebugUtilsMessengerCreateInfoEXT debug_utils_messenger_create_info =
			{
				VK_STRUCTURE_TYPE_DEBUG_UTILS_MESSENGER_CREATE_INFO_EXT,
				0U,
				0U,
				VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT | VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT,
				VK_DEBUG_UTILS_MESSAGE_TYPE_VALIDATION_BIT_EXT | VK_DEBUG_UTILS_MESSAGE_TYPE_PERFORMANCE_BIT_EXT,
				vulkan_debug_utils_messenger_callback,
				NULL};

		VkResult res_create_debug_utils_messenger = pfn_vk_create_debug_utils_messenger(vulkan_instance, &debug_utils_messenger_create_info, vulkan_allocation_callbacks, &vulkan_messenge);
		assert(VK_SUCCESS == res_create_debug_utils_messenger);
	}
#endif

	VkPhysicalDevice vulkan_physical_device = VK_NULL_HANDLE;
	uint32_t vulkan_min_uniform_buffer_offset_alignment;
	uint32_t vulkan_min_storage_buffer_offset_alignment;
	uint32_t vulkan_optimal_buffer_copy_offset_alignment;
	uint32_t vulkan_optimal_buffer_copy_row_pitch_alignment;
	{
		PFN_vkEnumeratePhysicalDevices pfn_enumerate_physical_devices = reinterpret_cast<PFN_vkEnumeratePhysicalDevices>(pfn_get_instance_proc_addr(vulkan_instance, "vkEnumeratePhysicalDevices"));
		assert(NULL != pfn_enumerate_physical_devices);
		PFN_vkGetPhysicalDeviceProperties pfn_get_physical_device_properties = reinterpret_cast<PFN_vkGetPhysicalDeviceProperties>(pfn_get_instance_proc_addr(vulkan_instance, "vkGetPhysicalDeviceProperties"));
		assert(NULL != pfn_get_physical_device_properties);

		uint32_t physical_device_count_1 = uint32_t(-1);
		VkResult res_enumerate_physical_devices_1 = pfn_enumerate_physical_devices(vulkan_instance, &physical_device_count_1, NULL);
		assert(VK_SUCCESS == res_enumerate_physical_devices_1 && 0U < physical_device_count_1);

		VkPhysicalDevice *physical_devices = static_cast<VkPhysicalDevice *>(malloc(sizeof(VkPhysicalDevice) * physical_device_count_1));
		assert(NULL != physical_devices);

		uint32_t physical_device_count_2 = physical_device_count_1;
		VkResult res_enumerate_physical_devices_2 = pfn_enumerate_physical_devices(vulkan_instance, &physical_device_count_2, physical_devices);
		assert(VK_SUCCESS == res_enumerate_physical_devices_2 && physical_device_count_1 == physical_device_count_2);

		// The lower index may imply the user preference (e.g. VK_LAYER_MESA_device_select)
		uint32_t physical_device_index_first_discrete_gpu = uint32_t(-1);
		VkDeviceSize min_uniform_buffer_offset_alignment_first_discrete_gpu = VkDeviceSize(-1);
		VkDeviceSize min_storage_buffer_offset_alignment_first_discrete_gpu = VkDeviceSize(-1);
		VkDeviceSize optimal_buffer_copy_offset_alignment_first_discrete_gpu = VkDeviceSize(-1);
		VkDeviceSize optimal_buffer_copy_row_pitch_alignment_first_discrete_gpu = VkDeviceSize(-1);
		uint32_t physical_device_index_first_integrated_gpu = uint32_t(-1);
		VkDeviceSize min_uniform_buffer_offset_alignment_first_integrated_gpu = VkDeviceSize(-1);
		VkDeviceSize min_storage_buffer_offset_alignment_first_integrated_gpu = VkDeviceSize(-1);
		VkDeviceSize optimal_buffer_copy_offset_alignment_first_integrated_gpu = VkDeviceSize(-1);
		VkDeviceSize optimal_buffer_copy_row_pitch_alignment_first_integrated_gpu = VkDeviceSize(-1);
		for (uint32_t physical_device_index = 0; (uint32_t(-1) == physical_device_index_first_discrete_gpu) && (physical_device_index < physical_device_count_2); ++physical_device_index)
		{
			VkPhysicalDeviceProperties physical_device_properties;
			pfn_get_physical_device_properties(physical_devices[physical_device_index], &physical_device_properties);

			if (VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU == physical_device_properties.deviceType)
			{
				physical_device_index_first_discrete_gpu = physical_device_index;
				min_uniform_buffer_offset_alignment_first_discrete_gpu = physical_device_properties.limits.minUniformBufferOffsetAlignment;
				min_storage_buffer_offset_alignment_first_discrete_gpu = physical_device_properties.limits.minStorageBufferOffsetAlignment;
				optimal_buffer_copy_offset_alignment_first_discrete_gpu = physical_device_properties.limits.optimalBufferCopyOffsetAlignment;
				optimal_buffer_copy_row_pitch_alignment_first_discrete_gpu = physical_device_properties.limits.optimalBufferCopyRowPitchAlignment;
			}
			else if (VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU == physical_device_properties.deviceType)
			{
				physical_device_index_first_integrated_gpu = physical_device_index;
				min_uniform_buffer_offset_alignment_first_integrated_gpu = physical_device_properties.limits.minUniformBufferOffsetAlignment;
				min_storage_buffer_offset_alignment_first_integrated_gpu = physical_device_properties.limits.minStorageBufferOffsetAlignment;
				optimal_buffer_copy_offset_alignment_first_integrated_gpu = physical_device_properties.limits.optimalBufferCopyOffsetAlignment;
				optimal_buffer_copy_row_pitch_alignment_first_integrated_gpu = physical_device_properties.limits.optimalBufferCopyRowPitchAlignment;
			}
		}

		// The discrete gpu is preferred
		if (uint32_t(-1) != physical_device_index_first_discrete_gpu)
		{
			vulkan_physical_device = physical_devices[physical_device_index_first_discrete_gpu];
			assert(min_uniform_buffer_offset_alignment_first_discrete_gpu < static_cast<VkDeviceSize>(UINT32_MAX));
			assert(min_storage_buffer_offset_alignment_first_discrete_gpu < static_cast<VkDeviceSize>(UINT32_MAX));
			assert(optimal_buffer_copy_offset_alignment_first_discrete_gpu < static_cast<VkDeviceSize>(UINT32_MAX));
			assert(optimal_buffer_copy_row_pitch_alignment_first_discrete_gpu < static_cast<VkDeviceSize>(UINT32_MAX));
			vulkan_min_uniform_buffer_offset_alignment = static_cast<uint32_t>(min_uniform_buffer_offset_alignment_first_discrete_gpu);
			vulkan_min_storage_buffer_offset_alignment = static_cast<uint32_t>(min_storage_buffer_offset_alignment_first_discrete_gpu);
			vulkan_optimal_buffer_copy_offset_alignment = static_cast<uint32_t>(optimal_buffer_copy_offset_alignment_first_discrete_gpu);
			vulkan_optimal_buffer_copy_row_pitch_alignment = static_cast<uint32_t>(optimal_buffer_copy_row_pitch_alignment_first_discrete_gpu);
		}
		else if (uint32_t(-1) != physical_device_index_first_integrated_gpu)
		{
			vulkan_physical_device = physical_devices[physical_device_index_first_integrated_gpu];
			assert(min_uniform_buffer_offset_alignment_first_integrated_gpu < static_cast<VkDeviceSize>(UINT32_MAX));
			assert(min_storage_buffer_offset_alignment_first_integrated_gpu < static_cast<VkDeviceSize>(UINT32_MAX));
			assert(optimal_buffer_copy_offset_alignment_first_integrated_gpu < static_cast<VkDeviceSize>(UINT32_MAX));
			assert(optimal_buffer_copy_row_pitch_alignment_first_integrated_gpu < static_cast<VkDeviceSize>(UINT32_MAX));
			vulkan_min_uniform_buffer_offset_alignment = static_cast<uint32_t>(min_uniform_buffer_offset_alignment_first_integrated_gpu);
			vulkan_min_storage_buffer_offset_alignment = static_cast<uint32_t>(min_storage_buffer_offset_alignment_first_integrated_gpu);
			vulkan_optimal_buffer_copy_offset_alignment = static_cast<uint32_t>(optimal_buffer_copy_offset_alignment_first_integrated_gpu);
			vulkan_optimal_buffer_copy_row_pitch_alignment = static_cast<uint32_t>(optimal_buffer_copy_row_pitch_alignment_first_integrated_gpu);
		}
		else
		{
			assert(false);
		}

		free(physical_devices);
	}

	// https://github.com/ValveSoftware/dxvk
	// src/dxvk/dxvk_device.h
	// DxvkDevice::hasDedicatedTransferQueue
	bool vulkan_has_dedicated_transfer_queue = false;
	// nvpro-samples/shared_sources/nvvk/context_vk.cpp
	// Context::initDevice
	// m_queue_GCT
	// m_queue_CT
	// m_queue_transfer
	uint32_t vulkan_queue_graphics_family_index = VK_QUEUE_FAMILY_IGNORED;
	uint32_t vulkan_queue_transfer_family_index = VK_QUEUE_FAMILY_IGNORED;
	uint32_t vulkan_queue_graphics_queue_index = uint32_t(-1);
	uint32_t vulkan_queue_transfer_queue_index = uint32_t(-1);
	{
		PFN_vkGetPhysicalDeviceQueueFamilyProperties pfn_vk_get_physical_device_queue_family_properties = reinterpret_cast<PFN_vkGetPhysicalDeviceQueueFamilyProperties>(pfn_get_instance_proc_addr(vulkan_instance, "vkGetPhysicalDeviceQueueFamilyProperties"));
		assert(NULL != pfn_vk_get_physical_device_queue_family_properties);
		PFN_vkGetPhysicalDeviceWin32PresentationSupportKHR pfn_vk_get_physical_device_win32_presentation_support = reinterpret_cast<PFN_vkGetPhysicalDeviceWin32PresentationSupportKHR>(pfn_get_instance_proc_addr(vulkan_instance, "vkGetPhysicalDeviceWin32PresentationSupportKHR"));
		assert(NULL != pfn_vk_get_physical_device_win32_presentation_support);

		uint32_t queue_family_property_count_1 = uint32_t(-1);
		pfn_vk_get_physical_device_queue_family_properties(vulkan_physical_device, &queue_family_property_count_1, NULL);

		VkQueueFamilyProperties *queue_family_properties = static_cast<VkQueueFamilyProperties *>(malloc(sizeof(VkQueueFamilyProperties) * queue_family_property_count_1));
		assert(NULL != queue_family_properties);

		uint32_t queue_family_property_count_2 = queue_family_property_count_1;
		pfn_vk_get_physical_device_queue_family_properties(vulkan_physical_device, &queue_family_property_count_2, queue_family_properties);
		assert(queue_family_property_count_1 == queue_family_property_count_2);

		// TODO: support seperated present queue
		// src/dxvk/dxvk_adapter.cpp
		// DxvkAdapter::findQueueFamilies
		// src/d3d11/d3d11_swapchain.cpp
		// D3D11SwapChain::CreatePresenter

		for (uint32_t queue_family_index = 0U; queue_family_index < queue_family_property_count_2; ++queue_family_index)
		{
			VkQueueFamilyProperties queue_family_property = queue_family_properties[queue_family_index];

			if ((queue_family_property.queueFlags & VK_QUEUE_GRAPHICS_BIT) && pfn_vk_get_physical_device_win32_presentation_support(vulkan_physical_device, queue_family_index))
			{
				vulkan_queue_graphics_family_index = queue_family_index;
				vulkan_queue_graphics_queue_index = 0U;
				break;
			}
		}

		if (VK_QUEUE_FAMILY_IGNORED != vulkan_queue_graphics_family_index)
		{
			// Find transfer queues
			for (uint32_t queue_family_index = 0U; queue_family_index < queue_family_property_count_2; ++queue_family_index)
			{
				if ((vulkan_queue_graphics_family_index != queue_family_index) && (VK_QUEUE_TRANSFER_BIT == (queue_family_properties[queue_family_index].queueFlags & (VK_QUEUE_GRAPHICS_BIT | VK_QUEUE_COMPUTE_BIT | VK_QUEUE_TRANSFER_BIT))))
				{
					vulkan_queue_transfer_family_index = queue_family_index;
					vulkan_queue_transfer_queue_index = 0U;
					vulkan_has_dedicated_transfer_queue = true;
					break;
				}
			}

			// Fallback to other graphics/compute queues
			// By vkspec, "either GRAPHICS or COMPUTE implies TRANSFER". This means TRANSFER is optional.
			if (VK_QUEUE_FAMILY_IGNORED == vulkan_queue_transfer_family_index)
			{
				for (uint32_t queue_family_index = 0U; queue_family_index < queue_family_property_count_2; ++queue_family_index)
				{
					if ((vulkan_queue_graphics_family_index != queue_family_index) && (0 != (queue_family_properties[queue_family_index].queueFlags & (VK_QUEUE_GRAPHICS_BIT | VK_QUEUE_COMPUTE_BIT))))
					{
						vulkan_queue_transfer_family_index = queue_family_index;
						vulkan_queue_transfer_queue_index = 0U;
						vulkan_has_dedicated_transfer_queue = true;
						break;
					}
				}
			}

			// Try the same queue family
			if (VK_QUEUE_FAMILY_IGNORED == vulkan_queue_transfer_family_index)
			{
				if (2U <= queue_family_properties[vulkan_queue_graphics_family_index].queueCount)
				{
					vulkan_queue_transfer_family_index = vulkan_queue_graphics_family_index;
					vulkan_queue_transfer_queue_index = 1U;
					vulkan_has_dedicated_transfer_queue = true;
				}
				else
				{
					vulkan_queue_transfer_family_index = VK_QUEUE_FAMILY_IGNORED;
					vulkan_queue_transfer_queue_index = uint32_t(-1);
					vulkan_has_dedicated_transfer_queue = false;
				}
			}
		}

		assert(VK_QUEUE_FAMILY_IGNORED != vulkan_queue_graphics_family_index && (!vulkan_has_dedicated_transfer_queue || VK_QUEUE_FAMILY_IGNORED != vulkan_queue_transfer_family_index));

		free(queue_family_properties);
	}

	VkSurfaceKHR vulkan_surface = VK_NULL_HANDLE;
	{
		PFN_vkCreateWin32SurfaceKHR pfn_create_win32_surface = reinterpret_cast<PFN_vkCreateWin32SurfaceKHR>(pfn_get_instance_proc_addr(vulkan_instance, "vkCreateWin32SurfaceKHR"));

		VkWin32SurfaceCreateInfoKHR win32_surface_create_info;
		win32_surface_create_info.sType = VK_STRUCTURE_TYPE_WIN32_SURFACE_CREATE_INFO_KHR;
		win32_surface_create_info.pNext = NULL;
		win32_surface_create_info.flags = 0U;
		win32_surface_create_info.hinstance = reinterpret_cast<HINSTANCE>(GetClassLongPtrW(hWnd, GCLP_HMODULE));
		win32_surface_create_info.hwnd = hWnd;
		VkResult res_create_win32_surface = pfn_create_win32_surface(vulkan_instance, &win32_surface_create_info, vulkan_allocation_callbacks, &vulkan_surface);
		assert(VK_SUCCESS == res_create_win32_surface);

		PFN_vkGetPhysicalDeviceSurfaceSupportKHR pfn_get_physical_device_surface_support = reinterpret_cast<PFN_vkGetPhysicalDeviceSurfaceSupportKHR>(pfn_get_instance_proc_addr(vulkan_instance, "vkGetPhysicalDeviceSurfaceSupportKHR"));
		VkBool32 supported;
		VkResult res_get_physical_device_surface_support = pfn_get_physical_device_surface_support(vulkan_physical_device, vulkan_queue_graphics_queue_index, vulkan_surface, &supported);
		assert(VK_SUCCESS == res_get_physical_device_surface_support);
		assert(VK_FALSE != supported);
	}

	VkFormat vulkan_swapchain_image_format = VK_FORMAT_UNDEFINED;
	{
		PFN_vkGetPhysicalDeviceSurfaceFormatsKHR pfn_get_physical_device_surface_formats = reinterpret_cast<PFN_vkGetPhysicalDeviceSurfaceFormatsKHR>(pfn_get_instance_proc_addr(vulkan_instance, "vkGetPhysicalDeviceSurfaceFormatsKHR"));
		assert(NULL != pfn_get_physical_device_surface_formats);

		uint32_t surface_format_count_1 = uint32_t(-1);
		pfn_get_physical_device_surface_formats(vulkan_physical_device, vulkan_surface, &surface_format_count_1, NULL);

		VkSurfaceFormatKHR *surface_formats = static_cast<VkSurfaceFormatKHR *>(malloc(sizeof(VkQueueFamilyProperties) * surface_format_count_1));
		assert(NULL != surface_formats);

		uint32_t surface_format_count_2 = surface_format_count_1;
		pfn_get_physical_device_surface_formats(vulkan_physical_device, vulkan_surface, &surface_format_count_2, surface_formats);
		assert(surface_format_count_1 == surface_format_count_2);

		assert(surface_format_count_2 >= 1U);
		if (VK_FORMAT_UNDEFINED != surface_formats[0].format)
		{
			vulkan_swapchain_image_format = surface_formats[0].format;
		}
		else
		{
			vulkan_swapchain_image_format = VK_FORMAT_B8G8R8A8_UNORM;
		}

		assert(VK_COLOR_SPACE_SRGB_NONLINEAR_KHR == surface_formats[0].colorSpace);

		free(surface_formats);
	}

	bool physical_device_feature_texture_compression_ASTC_LDR;
	bool physical_device_feature_texture_compression_BC;
	VkDevice vulkan_device = VK_NULL_HANDLE;
	{
		PFN_vkGetPhysicalDeviceFeatures pfn_vk_get_physical_device_features = reinterpret_cast<PFN_vkGetPhysicalDeviceFeatures>(pfn_get_instance_proc_addr(vulkan_instance, "vkGetPhysicalDeviceFeatures"));
		assert(NULL != pfn_vk_get_physical_device_features);
		PFN_vkCreateDevice pfn_vk_create_device = reinterpret_cast<PFN_vkCreateDevice>(pfn_get_instance_proc_addr(vulkan_instance, "vkCreateDevice"));
		assert(NULL != pfn_vk_create_device);

		VkPhysicalDeviceFeatures physical_device_features;
		pfn_vk_get_physical_device_features(vulkan_physical_device, &physical_device_features);

		physical_device_feature_texture_compression_ASTC_LDR = (physical_device_features.textureCompressionASTC_LDR != VK_FALSE) ? true : false;
		physical_device_feature_texture_compression_BC = (physical_device_features.textureCompressionBC != VK_FALSE) ? true : false;

		struct VkDeviceCreateInfo device_create_info;
		device_create_info.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO;
		device_create_info.pNext = NULL;
		device_create_info.flags = 0U;
		float queue_graphics_priority = 1.0F;
		float queue_transfer_priority = 1.0F;
		float queue_graphics_transfer_priorities[2] = {queue_graphics_priority, queue_transfer_priority};
		struct VkDeviceQueueCreateInfo device_queue_create_infos[2];
		if (vulkan_has_dedicated_transfer_queue)
		{
			if (vulkan_queue_graphics_family_index != vulkan_queue_transfer_family_index)
			{
				device_queue_create_infos[0].sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
				device_queue_create_infos[0].pNext = NULL;
				device_queue_create_infos[0].flags = 0U;
				device_queue_create_infos[0].queueFamilyIndex = vulkan_queue_graphics_family_index;
				assert(0U == vulkan_queue_graphics_queue_index);
				device_queue_create_infos[0].queueCount = 1U;
				device_queue_create_infos[0].pQueuePriorities = &queue_graphics_priority;

				device_queue_create_infos[1].sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
				device_queue_create_infos[1].pNext = NULL;
				device_queue_create_infos[1].flags = 0U;
				device_queue_create_infos[1].queueFamilyIndex = vulkan_queue_transfer_family_index;
				assert(0U == vulkan_queue_transfer_queue_index);
				device_queue_create_infos[1].queueCount = 1U;
				device_queue_create_infos[1].pQueuePriorities = &queue_transfer_priority;

				device_create_info.pQueueCreateInfos = device_queue_create_infos;
				device_create_info.queueCreateInfoCount = 2U;
			}
			else
			{
				device_queue_create_infos[0].sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
				device_queue_create_infos[0].pNext = NULL;
				device_queue_create_infos[0].flags = 0U;
				device_queue_create_infos[0].queueFamilyIndex = vulkan_queue_graphics_family_index;
				assert(0U == vulkan_queue_graphics_queue_index);
				assert(1U == vulkan_queue_transfer_queue_index);
				device_queue_create_infos[0].queueCount = 2U;
				device_queue_create_infos[0].pQueuePriorities = queue_graphics_transfer_priorities;
				device_create_info.pQueueCreateInfos = device_queue_create_infos;
				device_create_info.queueCreateInfoCount = 1U;
			}
		}
		else
		{
			device_queue_create_infos[0].sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
			device_queue_create_infos[0].pNext = NULL;
			device_queue_create_infos[0].flags = 0U;
			device_queue_create_infos[0].queueFamilyIndex = vulkan_queue_graphics_family_index;
			assert(0U == vulkan_queue_graphics_queue_index);
			device_queue_create_infos[0].queueCount = 1U;
			device_queue_create_infos[0].pQueuePriorities = &queue_graphics_priority;
			device_create_info.pQueueCreateInfos = device_queue_create_infos;
			device_create_info.queueCreateInfoCount = 1U;
		}
		device_create_info.enabledLayerCount = 0U;
		device_create_info.ppEnabledLayerNames = NULL;
		char const *enabled_extension_names[1] = {VK_KHR_SWAPCHAIN_EXTENSION_NAME};
		device_create_info.enabledExtensionCount = 1U;
		device_create_info.ppEnabledExtensionNames = enabled_extension_names;
		device_create_info.pEnabledFeatures = &physical_device_features;

		VkResult res_create_device = pfn_vk_create_device(vulkan_physical_device, &device_create_info, vulkan_allocation_callbacks, &vulkan_device);
		assert(VK_SUCCESS == res_create_device);
	}

	PFN_vkGetDeviceProcAddr pfn_get_device_proc_addr = reinterpret_cast<PFN_vkGetDeviceProcAddr>(pfn_get_instance_proc_addr(vulkan_instance, "vkGetDeviceProcAddr"));
	assert(NULL != pfn_get_device_proc_addr);

	pfn_get_device_proc_addr = reinterpret_cast<PFN_vkGetDeviceProcAddr>(pfn_get_device_proc_addr(vulkan_device, "vkGetDeviceProcAddr"));
	assert(NULL != pfn_get_device_proc_addr);

	VkQueue vulkan_queue_graphics = VK_NULL_HANDLE;
	VkQueue vulkan_queue_transfer = VK_NULL_HANDLE;
	{
		PFN_vkGetDeviceQueue vk_get_device_queue = reinterpret_cast<PFN_vkGetDeviceQueue>(pfn_get_device_proc_addr(vulkan_device, "vkGetDeviceQueue"));
		assert(NULL != vk_get_device_queue);

		vk_get_device_queue(vulkan_device, vulkan_queue_graphics_family_index, vulkan_queue_graphics_queue_index, &vulkan_queue_graphics);

		if (vulkan_has_dedicated_transfer_queue)
		{
			assert(VK_QUEUE_FAMILY_IGNORED != vulkan_queue_transfer_family_index);
			assert(uint32_t(-1) != vulkan_queue_transfer_queue_index);
			vk_get_device_queue(vulkan_device, vulkan_queue_transfer_family_index, vulkan_queue_transfer_queue_index, &vulkan_queue_transfer);
		}
	}
	assert(VK_NULL_HANDLE != vulkan_queue_graphics && (!vulkan_has_dedicated_transfer_queue || VK_NULL_HANDLE != vulkan_queue_transfer));

	VkCommandPool vulkan_command_pools[FRAME_THROTTLING_COUNT];
	VkCommandBuffer vulkan_command_buffers[FRAME_THROTTLING_COUNT];
	VkSemaphore vulkan_semaphores_acquire_next_image[FRAME_THROTTLING_COUNT];
	VkSemaphore vulkan_semaphores_queue_submit[FRAME_THROTTLING_COUNT];
	VkFence vulkan_fences[FRAME_THROTTLING_COUNT];
	{
		PFN_vkCreateCommandPool pfn_create_command_pool = reinterpret_cast<PFN_vkCreateCommandPool>(pfn_get_device_proc_addr(vulkan_device, "vkCreateCommandPool"));
		assert(NULL != pfn_create_command_pool);
		PFN_vkAllocateCommandBuffers pfn_allocate_command_buffers = reinterpret_cast<PFN_vkAllocateCommandBuffers>(pfn_get_device_proc_addr(vulkan_device, "vkAllocateCommandBuffers"));
		assert(NULL != pfn_allocate_command_buffers);
		PFN_vkCreateFence pfn_create_fence = reinterpret_cast<PFN_vkCreateFence>(pfn_get_device_proc_addr(vulkan_device, "vkCreateFence"));
		assert(NULL != pfn_create_fence);
		PFN_vkCreateSemaphore pfn_create_semaphore = reinterpret_cast<PFN_vkCreateSemaphore>(pfn_get_device_proc_addr(vulkan_device, "vkCreateSemaphore"));
		assert(NULL != pfn_create_semaphore);

		for (uint32_t frame_throtting_index = 0U; frame_throtting_index < FRAME_THROTTLING_COUNT; ++frame_throtting_index)
		{
			VkCommandPoolCreateInfo command_pool_create_info;
			command_pool_create_info.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
			command_pool_create_info.pNext = NULL;
			command_pool_create_info.flags = 0U;
			command_pool_create_info.queueFamilyIndex = vulkan_queue_graphics_family_index;
			VkResult res_create_command_pool = pfn_create_command_pool(vulkan_device, &command_pool_create_info, vulkan_allocation_callbacks, &vulkan_command_pools[frame_throtting_index]);
			assert(VK_SUCCESS == res_create_command_pool);

			VkCommandBufferAllocateInfo command_buffer_allocate_info;
			command_buffer_allocate_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
			command_buffer_allocate_info.pNext = NULL;
			command_buffer_allocate_info.commandPool = vulkan_command_pools[frame_throtting_index];
			command_buffer_allocate_info.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
			command_buffer_allocate_info.commandBufferCount = 1U;
			VkResult res_allocate_command_buffers = pfn_allocate_command_buffers(vulkan_device, &command_buffer_allocate_info, &vulkan_command_buffers[frame_throtting_index]);
			assert(VK_SUCCESS == res_allocate_command_buffers);

			VkSemaphoreCreateInfo semaphore_create_info_acquire_next_image;
			semaphore_create_info_acquire_next_image.sType = VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO;
			semaphore_create_info_acquire_next_image.pNext = NULL;
			semaphore_create_info_acquire_next_image.flags = 0U;
			VkResult res_create_semaphore_acquire_next_image = pfn_create_semaphore(vulkan_device, &semaphore_create_info_acquire_next_image, vulkan_allocation_callbacks, &vulkan_semaphores_acquire_next_image[frame_throtting_index]);
			assert(VK_SUCCESS == res_create_semaphore_acquire_next_image);

			VkSemaphoreCreateInfo semaphore_create_info_queue_submit;
			semaphore_create_info_queue_submit.sType = VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO;
			semaphore_create_info_queue_submit.pNext = NULL;
			semaphore_create_info_queue_submit.flags = 0U;
			VkResult res_create_semaphore_queue_submit = pfn_create_semaphore(vulkan_device, &semaphore_create_info_queue_submit, vulkan_allocation_callbacks, &vulkan_semaphores_queue_submit[frame_throtting_index]);
			assert(VK_SUCCESS == res_create_semaphore_queue_submit);

			VkFenceCreateInfo fence_create_info;
			fence_create_info.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
			fence_create_info.pNext = NULL;
			fence_create_info.flags = VK_FENCE_CREATE_SIGNALED_BIT;
			VkResult res_create_fence = pfn_create_fence(vulkan_device, &fence_create_info, vulkan_allocation_callbacks, &vulkan_fences[frame_throtting_index]);
			assert(VK_SUCCESS == res_create_fence);
		}
	}

	VkDeviceSize vulkan_upload_ring_buffer_size;
	VkBuffer vulkan_upload_ring_buffer;
	VkDeviceMemory vulkan_upload_ring_buffer_device_memory;
	void *vulkan_upload_ring_buffer_device_memory_pointer;
	uint32_t vulkan_color_input_transient_attachment_memory_index;
	VkFormat vulkan_depth_format;
	uint32_t vulkan_depth_stencil_transient_attachment_memory_index;
	{
		PFN_vkGetPhysicalDeviceMemoryProperties pfn_get_physical_device_memory_properties = reinterpret_cast<PFN_vkGetPhysicalDeviceMemoryProperties>(pfn_get_instance_proc_addr(vulkan_instance, "vkGetPhysicalDeviceMemoryProperties"));
		PFN_vkCreateBuffer pfn_create_buffer = reinterpret_cast<PFN_vkCreateBuffer>(pfn_get_device_proc_addr(vulkan_device, "vkCreateBuffer"));
		PFN_vkGetBufferMemoryRequirements pfn_get_buffer_memory_requirements = reinterpret_cast<PFN_vkGetBufferMemoryRequirements>(pfn_get_device_proc_addr(vulkan_device, "vkGetBufferMemoryRequirements"));
		PFN_vkAllocateMemory pfn_allocate_memory = reinterpret_cast<PFN_vkAllocateMemory>(pfn_get_device_proc_addr(vulkan_device, "vkAllocateMemory"));
		PFN_vkMapMemory pfn_map_memory = reinterpret_cast<PFN_vkMapMemory>(pfn_get_device_proc_addr(vulkan_device, "vkMapMemory"));
		PFN_vkBindBufferMemory pfn_bind_buffer_memory = reinterpret_cast<PFN_vkBindBufferMemory>(pfn_get_device_proc_addr(vulkan_device, "vkBindBufferMemory"));
		PFN_vkGetPhysicalDeviceFormatProperties pfn_get_physical_device_format_properties = reinterpret_cast<PFN_vkGetPhysicalDeviceFormatProperties>(pfn_get_instance_proc_addr(vulkan_instance, "vkGetPhysicalDeviceFormatProperties"));
		PFN_vkCreateImage pfn_create_image = reinterpret_cast<PFN_vkCreateImage>(pfn_get_device_proc_addr(vulkan_device, "vkCreateImage"));
		PFN_vkGetImageMemoryRequirements pfn_get_image_memory_requirements = reinterpret_cast<PFN_vkGetImageMemoryRequirements>(pfn_get_device_proc_addr(vulkan_device, "vkGetImageMemoryRequirements"));
		PFN_vkDestroyImage pfn_destroy_image = reinterpret_cast<PFN_vkDestroyImage>(pfn_get_device_proc_addr(vulkan_device, "vkDestroyImage"));

		VkPhysicalDeviceMemoryProperties physical_device_memory_properties;
		pfn_get_physical_device_memory_properties(vulkan_physical_device, &physical_device_memory_properties);

		// upload ring buffer
		// NVIDIA Driver 128 MB
		// \[Gruen 2015\] [Holger Gruen. "Constant Buffers without Constant Pain." NVIDIA GameWorks Blog 2015.](https://developer.nvidia.com/content/constant-buffers-without-constant-pain-0)
		// AMD Special Pool 256MB
		// \[Sawicki 2018\] [Adam Sawicki. "Memory Management in Vulkan and DX12." GDC 2018.](https://gpuopen.com/events/gdc-2018-presentations)
		vulkan_upload_ring_buffer_size = (224ULL * 1024ULL * 1024ULL); // 224MB
		{
			vulkan_upload_ring_buffer = VK_NULL_HANDLE;
			VkDeviceSize memory_requirements_size = VkDeviceSize(-1);
			uint32_t memory_requirements_memory_type_bits = 0U;
			{
				struct VkBufferCreateInfo buffer_create_info;
				buffer_create_info.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
				buffer_create_info.pNext = NULL;
				buffer_create_info.flags = 0U;
				buffer_create_info.size = vulkan_upload_ring_buffer_size;
				buffer_create_info.usage = VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT;
				buffer_create_info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
				buffer_create_info.queueFamilyIndexCount = 0U;
				buffer_create_info.pQueueFamilyIndices = NULL;

				VkResult res_create_buffer = pfn_create_buffer(vulkan_device, &buffer_create_info, vulkan_allocation_callbacks, &vulkan_upload_ring_buffer);
				assert(VK_SUCCESS == res_create_buffer);

				struct VkMemoryRequirements memory_requirements;
				pfn_get_buffer_memory_requirements(vulkan_device, vulkan_upload_ring_buffer, &memory_requirements);
				memory_requirements_size = memory_requirements.size;
				memory_requirements_memory_type_bits = memory_requirements.memoryTypeBits;
			}

			// VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT
			// We use the "AMD Special Pool" for upload ring buffer
			uint32_t vulkan_upload_ring_buffer_memory_index = __internal_find_lowest_memory_type_index(&physical_device_memory_properties, memory_requirements_size, memory_requirements_memory_type_bits, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
			assert(VK_MAX_MEMORY_TYPES > vulkan_upload_ring_buffer_memory_index);
			assert(physical_device_memory_properties.memoryTypeCount > vulkan_upload_ring_buffer_memory_index);

			vulkan_upload_ring_buffer_device_memory = VK_NULL_HANDLE;
			{
				VkMemoryAllocateInfo memory_allocate_info;
				memory_allocate_info.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
				memory_allocate_info.pNext = NULL;
				memory_allocate_info.allocationSize = memory_requirements_size;
				memory_allocate_info.memoryTypeIndex = vulkan_upload_ring_buffer_memory_index;
				VkResult res_allocate_memory = pfn_allocate_memory(vulkan_device, &memory_allocate_info, vulkan_allocation_callbacks, &vulkan_upload_ring_buffer_device_memory);
				assert(VK_SUCCESS == res_allocate_memory);
			}

			vulkan_upload_ring_buffer_device_memory_pointer = NULL;
			VkResult res_map_memory = pfn_map_memory(vulkan_device, vulkan_upload_ring_buffer_device_memory, 0U, vulkan_upload_ring_buffer_size, 0U, &vulkan_upload_ring_buffer_device_memory_pointer);
			assert(VK_SUCCESS == res_map_memory);

			VkResult res_bind_buffer_memory = pfn_bind_buffer_memory(vulkan_device, vulkan_upload_ring_buffer, vulkan_upload_ring_buffer_device_memory, 0U);
			assert(VK_SUCCESS == res_bind_buffer_memory);
		}

		// color attachment
		vulkan_color_input_transient_attachment_memory_index = VK_MAX_MEMORY_TYPES;
		{
			VkDeviceSize memory_requirements_size = VkDeviceSize(-1);
			uint32_t memory_requirements_memory_type_bits = 0U;
			{
				struct VkImageCreateInfo image_create_info_depth_stencil_transient_tiling_optimal;
				image_create_info_depth_stencil_transient_tiling_optimal.sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;
				image_create_info_depth_stencil_transient_tiling_optimal.pNext = NULL;
				image_create_info_depth_stencil_transient_tiling_optimal.flags = 0U;
				image_create_info_depth_stencil_transient_tiling_optimal.imageType = VK_IMAGE_TYPE_2D;
				image_create_info_depth_stencil_transient_tiling_optimal.format = VK_FORMAT_R8G8B8A8_UNORM;
				image_create_info_depth_stencil_transient_tiling_optimal.extent.width = 8U;
				image_create_info_depth_stencil_transient_tiling_optimal.extent.height = 8U;
				image_create_info_depth_stencil_transient_tiling_optimal.extent.depth = 1U;
				image_create_info_depth_stencil_transient_tiling_optimal.mipLevels = 1U;
				image_create_info_depth_stencil_transient_tiling_optimal.arrayLayers = 1U;
				image_create_info_depth_stencil_transient_tiling_optimal.samples = VK_SAMPLE_COUNT_1_BIT;
				image_create_info_depth_stencil_transient_tiling_optimal.tiling = VK_IMAGE_TILING_OPTIMAL;
				image_create_info_depth_stencil_transient_tiling_optimal.usage = VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_INPUT_ATTACHMENT_BIT | VK_IMAGE_USAGE_TRANSIENT_ATTACHMENT_BIT;
				image_create_info_depth_stencil_transient_tiling_optimal.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
				image_create_info_depth_stencil_transient_tiling_optimal.queueFamilyIndexCount = 0U;
				image_create_info_depth_stencil_transient_tiling_optimal.pQueueFamilyIndices = NULL;
				image_create_info_depth_stencil_transient_tiling_optimal.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;

				VkImage dummy_img;
				VkResult res_create_image = pfn_create_image(vulkan_device, &image_create_info_depth_stencil_transient_tiling_optimal, vulkan_allocation_callbacks, &dummy_img);
				assert(VK_SUCCESS == res_create_image);

				struct VkMemoryRequirements memory_requirements;
				pfn_get_image_memory_requirements(vulkan_device, dummy_img, &memory_requirements);
				memory_requirements_size = memory_requirements.size;
				memory_requirements_memory_type_bits = memory_requirements.memoryTypeBits;

				pfn_destroy_image(vulkan_device, dummy_img, vulkan_allocation_callbacks);
			}

			// The lower index indicates the more performance
			vulkan_color_input_transient_attachment_memory_index = __internal_find_lowest_memory_type_index(&physical_device_memory_properties, memory_requirements_size, memory_requirements_memory_type_bits, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_MEMORY_PROPERTY_LAZILY_ALLOCATED_BIT);
			assert(VK_MAX_MEMORY_TYPES > vulkan_color_input_transient_attachment_memory_index);
			assert(physical_device_memory_properties.memoryTypeCount > vulkan_color_input_transient_attachment_memory_index);
		}

		// depth attachment
		vulkan_depth_format = VK_FORMAT_UNDEFINED;
		vulkan_depth_stencil_transient_attachment_memory_index = VK_MAX_MEMORY_TYPES;
		{
			struct VkFormatProperties format_properties;
			vulkan_depth_format = VK_FORMAT_D32_SFLOAT;
			pfn_get_physical_device_format_properties(vulkan_physical_device, vulkan_depth_format, &format_properties);
			if (0 == (format_properties.optimalTilingFeatures & VK_FORMAT_FEATURE_DEPTH_STENCIL_ATTACHMENT_BIT))
			{
				vulkan_depth_format = VK_FORMAT_X8_D24_UNORM_PACK32;
				pfn_get_physical_device_format_properties(vulkan_physical_device, vulkan_depth_format, &format_properties);
				assert(format_properties.optimalTilingFeatures & VK_FORMAT_FEATURE_DEPTH_STENCIL_ATTACHMENT_BIT);
			}

			VkDeviceSize memory_requirements_size = VkDeviceSize(-1);
			uint32_t memory_requirements_memory_type_bits = 0U;
			{
				struct VkImageCreateInfo image_create_info_depth_stencil_transient_tiling_optimal;
				image_create_info_depth_stencil_transient_tiling_optimal.sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;
				image_create_info_depth_stencil_transient_tiling_optimal.pNext = NULL;
				image_create_info_depth_stencil_transient_tiling_optimal.flags = 0U;
				image_create_info_depth_stencil_transient_tiling_optimal.imageType = VK_IMAGE_TYPE_2D;
				image_create_info_depth_stencil_transient_tiling_optimal.format = vulkan_depth_format;
				image_create_info_depth_stencil_transient_tiling_optimal.extent.width = 8U;
				image_create_info_depth_stencil_transient_tiling_optimal.extent.height = 8U;
				image_create_info_depth_stencil_transient_tiling_optimal.extent.depth = 1U;
				image_create_info_depth_stencil_transient_tiling_optimal.mipLevels = 1U;
				image_create_info_depth_stencil_transient_tiling_optimal.arrayLayers = 1U;
				image_create_info_depth_stencil_transient_tiling_optimal.samples = VK_SAMPLE_COUNT_1_BIT;
				image_create_info_depth_stencil_transient_tiling_optimal.tiling = VK_IMAGE_TILING_OPTIMAL;
				image_create_info_depth_stencil_transient_tiling_optimal.usage = VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT | VK_IMAGE_USAGE_TRANSIENT_ATTACHMENT_BIT;
				image_create_info_depth_stencil_transient_tiling_optimal.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
				image_create_info_depth_stencil_transient_tiling_optimal.queueFamilyIndexCount = 0U;
				image_create_info_depth_stencil_transient_tiling_optimal.pQueueFamilyIndices = NULL;
				image_create_info_depth_stencil_transient_tiling_optimal.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;

				VkImage dummy_img;
				VkResult res_create_image = pfn_create_image(vulkan_device, &image_create_info_depth_stencil_transient_tiling_optimal, vulkan_allocation_callbacks, &dummy_img);
				assert(VK_SUCCESS == res_create_image);

				struct VkMemoryRequirements memory_requirements;
				pfn_get_image_memory_requirements(vulkan_device, dummy_img, &memory_requirements);
				memory_requirements_size = memory_requirements.size;
				memory_requirements_memory_type_bits = memory_requirements.memoryTypeBits;

				pfn_destroy_image(vulkan_device, dummy_img, vulkan_allocation_callbacks);
			}

			// The lower index indicates the more performance
			vulkan_depth_stencil_transient_attachment_memory_index = __internal_find_lowest_memory_type_index(&physical_device_memory_properties, memory_requirements_size, memory_requirements_memory_type_bits, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, VK_MEMORY_PROPERTY_LAZILY_ALLOCATED_BIT);
			assert(VK_MAX_MEMORY_TYPES > vulkan_depth_stencil_transient_attachment_memory_index);
			assert(physical_device_memory_properties.memoryTypeCount > vulkan_depth_stencil_transient_attachment_memory_index);
		}
	}

	uint32_t vulkan_upload_ring_buffer_begin[FRAME_THROTTLING_COUNT];
	uint32_t vulkan_upload_ring_buffer_end[FRAME_THROTTLING_COUNT];
	uint32_t vulkan_upload_ring_buffer_current[FRAME_THROTTLING_COUNT];
	{
		assert(vulkan_upload_ring_buffer_size < static_cast<VkDeviceSize>(UINT32_MAX));

		for (uint32_t frame_throtting_index = 0U; frame_throtting_index < FRAME_THROTTLING_COUNT; ++frame_throtting_index)
		{
			vulkan_upload_ring_buffer_begin[frame_throtting_index] = (static_cast<uint32_t>(vulkan_upload_ring_buffer_size) * frame_throtting_index) / FRAME_THROTTLING_COUNT;
			vulkan_upload_ring_buffer_end[frame_throtting_index] = (static_cast<uint32_t>(vulkan_upload_ring_buffer_size) * (frame_throtting_index + 1U)) / FRAME_THROTTLING_COUNT;
		}
	}

	VmaAllocator vulkan_asset_allocator;
	{
		VmaVulkanFunctions vulkan_functions = {};
		vulkan_functions.vkGetInstanceProcAddr = pfn_get_instance_proc_addr;
		vulkan_functions.vkGetDeviceProcAddr = pfn_get_device_proc_addr;

		VmaAllocatorCreateInfo allocator_create_info = {};
		allocator_create_info.vulkanApiVersion = vulkan_api_version;
		allocator_create_info.physicalDevice = vulkan_physical_device;
		allocator_create_info.device = vulkan_device;
		allocator_create_info.instance = vulkan_instance;
		allocator_create_info.pAllocationCallbacks = vulkan_allocation_callbacks;
		allocator_create_info.pVulkanFunctions = &vulkan_functions;

		vmaCreateAllocator(&allocator_create_info, &vulkan_asset_allocator);
	}

	PFN_vkWaitForFences pfn_wait_for_fences = reinterpret_cast<PFN_vkWaitForFences>(pfn_get_device_proc_addr(vulkan_device, "vkWaitForFences"));
	assert(NULL != pfn_wait_for_fences);
	PFN_vkResetFences pfn_reset_fences = reinterpret_cast<PFN_vkResetFences>(pfn_get_device_proc_addr(vulkan_device, "vkResetFences"));
	assert(NULL != pfn_reset_fences);
	PFN_vkResetCommandPool pfn_reset_command_pool = reinterpret_cast<PFN_vkResetCommandPool>(pfn_get_device_proc_addr(vulkan_device, "vkResetCommandPool"));
	assert(NULL != pfn_reset_command_pool);
	PFN_vkBeginCommandBuffer pfn_begin_command_buffer = reinterpret_cast<PFN_vkBeginCommandBuffer>(pfn_get_device_proc_addr(vulkan_device, "vkBeginCommandBuffer"));
	assert(NULL != pfn_begin_command_buffer);
	PFN_vkEndCommandBuffer pfn_end_command_buffer = reinterpret_cast<PFN_vkEndCommandBuffer>(pfn_get_device_proc_addr(vulkan_device, "vkEndCommandBuffer"));
	assert(NULL != pfn_end_command_buffer);
	PFN_vkAcquireNextImageKHR pfn_acquire_next_image = reinterpret_cast<PFN_vkAcquireNextImageKHR>(pfn_get_device_proc_addr(vulkan_device, "vkAcquireNextImageKHR"));
	assert(NULL != pfn_acquire_next_image);
	PFN_vkQueueSubmit pfn_queue_submit = reinterpret_cast<PFN_vkQueueSubmit>(pfn_get_device_proc_addr(vulkan_device, "vkQueueSubmit"));
	assert(NULL != pfn_queue_submit);
	PFN_vkQueuePresentKHR pfn_queue_present = reinterpret_cast<PFN_vkQueuePresentKHR>(pfn_get_device_proc_addr(vulkan_device, "vkQueuePresentKHR"));
	assert(NULL != pfn_queue_present);

	// Demo
	{
		class Demo demo;

		// Streaming
		{
			VkCommandPool vulkan_streaming_transfer_command_pool = VK_NULL_HANDLE;
			VkCommandBuffer vulkan_streaming_transfer_command_buffer = VK_NULL_HANDLE;
			VkCommandPool vulkan_streaming_graphics_command_pool = VK_NULL_HANDLE;
			VkCommandBuffer vulkan_streaming_graphics_command_buffer = VK_NULL_HANDLE;
			VkSemaphore vulkan_streaming_semaphore = VK_NULL_HANDLE;
			VkFence vulkan_streaming_fence = VK_NULL_HANDLE;
			{
				PFN_vkCreateCommandPool pfn_create_command_pool = reinterpret_cast<PFN_vkCreateCommandPool>(pfn_get_device_proc_addr(vulkan_device, "vkCreateCommandPool"));
				assert(NULL != pfn_create_command_pool);
				PFN_vkAllocateCommandBuffers pfn_allocate_command_buffers = reinterpret_cast<PFN_vkAllocateCommandBuffers>(pfn_get_device_proc_addr(vulkan_device, "vkAllocateCommandBuffers"));
				assert(NULL != pfn_allocate_command_buffers);
				PFN_vkCreateFence pfn_create_fence = reinterpret_cast<PFN_vkCreateFence>(pfn_get_device_proc_addr(vulkan_device, "vkCreateFence"));
				assert(NULL != pfn_create_fence);
				PFN_vkCreateSemaphore pfn_create_semaphore = reinterpret_cast<PFN_vkCreateSemaphore>(pfn_get_device_proc_addr(vulkan_device, "vkCreateSemaphore"));
				assert(NULL != pfn_create_semaphore);

				if (vulkan_has_dedicated_transfer_queue)
				{
					if (vulkan_queue_graphics_family_index != vulkan_queue_transfer_family_index)
					{

						VkCommandPoolCreateInfo transfer_command_pool_create_info;
						transfer_command_pool_create_info.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
						transfer_command_pool_create_info.pNext = NULL;
						transfer_command_pool_create_info.flags = 0U;
						transfer_command_pool_create_info.queueFamilyIndex = vulkan_queue_transfer_family_index;
						VkResult res_transfer_create_command_pool = pfn_create_command_pool(vulkan_device, &transfer_command_pool_create_info, vulkan_allocation_callbacks, &vulkan_streaming_transfer_command_pool);
						assert(VK_SUCCESS == res_transfer_create_command_pool);

						VkCommandBufferAllocateInfo transfer_command_buffer_allocate_info;
						transfer_command_buffer_allocate_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
						transfer_command_buffer_allocate_info.pNext = NULL;
						transfer_command_buffer_allocate_info.commandPool = vulkan_streaming_transfer_command_pool;
						transfer_command_buffer_allocate_info.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
						transfer_command_buffer_allocate_info.commandBufferCount = 1U;
						VkResult res_transfer_allocate_command_buffers = pfn_allocate_command_buffers(vulkan_device, &transfer_command_buffer_allocate_info, &vulkan_streaming_transfer_command_buffer);
						assert(VK_SUCCESS == res_transfer_allocate_command_buffers);

						VkCommandPoolCreateInfo graphics_command_pool_create_info;
						graphics_command_pool_create_info.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
						graphics_command_pool_create_info.pNext = NULL;
						graphics_command_pool_create_info.flags = 0U;
						graphics_command_pool_create_info.queueFamilyIndex = vulkan_queue_graphics_family_index;
						VkResult res_graphics_create_command_pool = pfn_create_command_pool(vulkan_device, &graphics_command_pool_create_info, vulkan_allocation_callbacks, &vulkan_streaming_graphics_command_pool);
						assert(VK_SUCCESS == res_graphics_create_command_pool);

						VkCommandBufferAllocateInfo graphics_command_buffer_allocate_info;
						graphics_command_buffer_allocate_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
						graphics_command_buffer_allocate_info.pNext = NULL;
						graphics_command_buffer_allocate_info.commandPool = vulkan_streaming_graphics_command_pool;
						graphics_command_buffer_allocate_info.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
						graphics_command_buffer_allocate_info.commandBufferCount = 1U;
						VkResult res_graphics_allocate_command_buffers = pfn_allocate_command_buffers(vulkan_device, &graphics_command_buffer_allocate_info, &vulkan_streaming_graphics_command_buffer);
						assert(VK_SUCCESS == res_graphics_allocate_command_buffers);

						VkSemaphoreCreateInfo semaphore_create_info;
						semaphore_create_info.sType = VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO;
						semaphore_create_info.pNext = NULL;
						semaphore_create_info.flags = 0U;
						VkResult res_create_semaphore = pfn_create_semaphore(vulkan_device, &semaphore_create_info, vulkan_allocation_callbacks, &vulkan_streaming_semaphore);
						assert(VK_SUCCESS == res_create_semaphore);
					}
					else
					{

						VkCommandPoolCreateInfo transfer_command_pool_create_info;
						transfer_command_pool_create_info.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
						transfer_command_pool_create_info.pNext = NULL;
						transfer_command_pool_create_info.flags = 0U;
						transfer_command_pool_create_info.queueFamilyIndex = vulkan_queue_transfer_family_index;
						VkResult res_transfer_create_command_pool = pfn_create_command_pool(vulkan_device, &transfer_command_pool_create_info, vulkan_allocation_callbacks, &vulkan_streaming_transfer_command_pool);
						assert(VK_SUCCESS == res_transfer_create_command_pool);

						VkCommandBufferAllocateInfo transfer_command_buffer_allocate_info;
						transfer_command_buffer_allocate_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
						transfer_command_buffer_allocate_info.pNext = NULL;
						transfer_command_buffer_allocate_info.commandPool = vulkan_streaming_transfer_command_pool;
						transfer_command_buffer_allocate_info.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
						transfer_command_buffer_allocate_info.commandBufferCount = 1U;
						VkResult res_transfer_allocate_command_buffers = pfn_allocate_command_buffers(vulkan_device, &transfer_command_buffer_allocate_info, &vulkan_streaming_transfer_command_buffer);
						assert(VK_SUCCESS == res_transfer_allocate_command_buffers);
					}
				}
				else
				{

					VkCommandPoolCreateInfo graphics_command_pool_create_info;
					graphics_command_pool_create_info.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
					graphics_command_pool_create_info.pNext = NULL;
					graphics_command_pool_create_info.flags = 0U;
					graphics_command_pool_create_info.queueFamilyIndex = vulkan_queue_graphics_family_index;
					VkResult res_graphics_create_command_pool = pfn_create_command_pool(vulkan_device, &graphics_command_pool_create_info, vulkan_allocation_callbacks, &vulkan_streaming_graphics_command_pool);
					assert(VK_SUCCESS == res_graphics_create_command_pool);

					VkCommandBufferAllocateInfo graphics_command_buffer_allocate_info;
					graphics_command_buffer_allocate_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
					graphics_command_buffer_allocate_info.pNext = NULL;
					graphics_command_buffer_allocate_info.commandPool = vulkan_streaming_graphics_command_pool;
					graphics_command_buffer_allocate_info.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
					graphics_command_buffer_allocate_info.commandBufferCount = 1U;
					VkResult res_graphics_allocate_command_buffers = pfn_allocate_command_buffers(vulkan_device, &graphics_command_buffer_allocate_info, &vulkan_streaming_graphics_command_buffer);
					assert(VK_SUCCESS == res_graphics_allocate_command_buffers);
				}

				VkFenceCreateInfo fence_create_info;
				fence_create_info.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
				fence_create_info.pNext = NULL;
				fence_create_info.flags = 0U;
				VkResult res_create_fence = pfn_create_fence(vulkan_device, &fence_create_info, vulkan_allocation_callbacks, &vulkan_streaming_fence);
				assert(VK_SUCCESS == res_create_fence);
			}

			VkDeviceSize vulkan_staging_buffer_size;
			VkBuffer vulkan_staging_buffer;
			VkDeviceMemory vulkan_staging_buffer_device_memory;
			void *vulkan_staging_buffer_device_memory_pointer;
			uint32_t vulkan_asset_vertex_buffer_memory_index;
			uint32_t vulkan_asset_index_buffer_memory_index;
			uint32_t vulkan_asset_uniform_buffer_memory_index;
			uint32_t vulkan_asset_image_memory_index;
			{
				PFN_vkGetPhysicalDeviceMemoryProperties pfn_get_physical_device_memory_properties = reinterpret_cast<PFN_vkGetPhysicalDeviceMemoryProperties>(pfn_get_instance_proc_addr(vulkan_instance, "vkGetPhysicalDeviceMemoryProperties"));
				PFN_vkCreateBuffer pfn_create_buffer = reinterpret_cast<PFN_vkCreateBuffer>(pfn_get_device_proc_addr(vulkan_device, "vkCreateBuffer"));
				PFN_vkGetBufferMemoryRequirements pfn_get_buffer_memory_requirements = reinterpret_cast<PFN_vkGetBufferMemoryRequirements>(pfn_get_device_proc_addr(vulkan_device, "vkGetBufferMemoryRequirements"));
				PFN_vkAllocateMemory pfn_allocate_memory = reinterpret_cast<PFN_vkAllocateMemory>(pfn_get_device_proc_addr(vulkan_device, "vkAllocateMemory"));
				PFN_vkMapMemory pfn_map_memory = reinterpret_cast<PFN_vkMapMemory>(pfn_get_device_proc_addr(vulkan_device, "vkMapMemory"));
				PFN_vkBindBufferMemory pfn_bind_buffer_memory = reinterpret_cast<PFN_vkBindBufferMemory>(pfn_get_device_proc_addr(vulkan_device, "vkBindBufferMemory"));
				PFN_vkDestroyBuffer pfn_destroy_buffer = reinterpret_cast<PFN_vkDestroyBuffer>(pfn_get_device_proc_addr(vulkan_device, "vkDestroyBuffer"));
				PFN_vkGetPhysicalDeviceFormatProperties pfn_get_physical_device_format_properties = reinterpret_cast<PFN_vkGetPhysicalDeviceFormatProperties>(pfn_get_instance_proc_addr(vulkan_instance, "vkGetPhysicalDeviceFormatProperties"));
				PFN_vkCreateImage pfn_create_image = reinterpret_cast<PFN_vkCreateImage>(pfn_get_device_proc_addr(vulkan_device, "vkCreateImage"));
				PFN_vkGetImageMemoryRequirements pfn_get_image_memory_requirements = reinterpret_cast<PFN_vkGetImageMemoryRequirements>(pfn_get_device_proc_addr(vulkan_device, "vkGetImageMemoryRequirements"));
				PFN_vkDestroyImage pfn_destroy_image = reinterpret_cast<PFN_vkDestroyImage>(pfn_get_device_proc_addr(vulkan_device, "vkDestroyImage"));

				VkPhysicalDeviceMemoryProperties physical_device_memory_properties;
				pfn_get_physical_device_memory_properties(vulkan_physical_device, &physical_device_memory_properties);

				// staging buffer
				vulkan_staging_buffer_size = 64ULL * 1024ULL * 1024ULL;
				{
					vulkan_staging_buffer = VK_NULL_HANDLE;
					VkDeviceSize memory_requirements_size = VkDeviceSize(-1);
					uint32_t memory_requirements_memory_type_bits = 0U;
					{
						struct VkBufferCreateInfo buffer_create_info;
						buffer_create_info.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
						buffer_create_info.pNext = NULL;
						buffer_create_info.flags = 0U;
						buffer_create_info.size = vulkan_staging_buffer_size;
						buffer_create_info.usage = VK_BUFFER_USAGE_TRANSFER_SRC_BIT;
						buffer_create_info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
						buffer_create_info.queueFamilyIndexCount = 0U;
						buffer_create_info.pQueueFamilyIndices = NULL;

						VkResult res_create_buffer = pfn_create_buffer(vulkan_device, &buffer_create_info, vulkan_allocation_callbacks, &vulkan_staging_buffer);
						assert(VK_SUCCESS == res_create_buffer);

						struct VkMemoryRequirements memory_requirements;
						pfn_get_buffer_memory_requirements(vulkan_device, vulkan_staging_buffer, &memory_requirements);
						memory_requirements_size = memory_requirements.size;
						memory_requirements_memory_type_bits = memory_requirements.memoryTypeBits;
					}

					// Do NOT use "VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT"
					// We leave the "AMD Special Pool" for upload ring buffer
					uint32_t vulkan_staging_buffer_memory_index = __internal_find_lowest_memory_type_index(&physical_device_memory_properties, memory_requirements_size, memory_requirements_memory_type_bits, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);
					assert(VK_MAX_MEMORY_TYPES > vulkan_staging_buffer_memory_index);
					assert(physical_device_memory_properties.memoryTypeCount > vulkan_staging_buffer_memory_index);

					vulkan_staging_buffer_device_memory = VK_NULL_HANDLE;
					{
						VkMemoryAllocateInfo memory_allocate_info;
						memory_allocate_info.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
						memory_allocate_info.pNext = NULL;
						memory_allocate_info.allocationSize = memory_requirements_size;
						memory_allocate_info.memoryTypeIndex = vulkan_staging_buffer_memory_index;
						VkResult res_allocate_memory = pfn_allocate_memory(vulkan_device, &memory_allocate_info, vulkan_allocation_callbacks, &vulkan_staging_buffer_device_memory);
						assert(VK_SUCCESS == res_allocate_memory);
					}

					VkResult res_map_memory = pfn_map_memory(vulkan_device, vulkan_staging_buffer_device_memory, 0U, vulkan_staging_buffer_size, 0U, &vulkan_staging_buffer_device_memory_pointer);
					assert(VK_SUCCESS == res_map_memory);

					VkResult res_bind_buffer_memory = pfn_bind_buffer_memory(vulkan_device, vulkan_staging_buffer, vulkan_staging_buffer_device_memory, 0U);
					assert(VK_SUCCESS == res_bind_buffer_memory);
				}

				// https://www.khronos.org/registry/vulkan/specs/1.0/html/chap13.html#VkMemoryRequirements
				// If buffer is a VkBuffer not created with the VK_BUFFER_CREATE_SPARSE_BINDING_BIT bit set, or if image is linear image,
				// then the memoryTypeBits member always contains at least one bit set corresponding to a VkMemoryType with a
				// propertyFlags that has both the VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT bit and the
				// VK_MEMORY_PROPERTY_HOST_COHERENT_BIT bit set. In other words, mappable coherent memory can always be attached to
				// these objects.

				// https://www.khronos.org/registry/vulkan/specs/1.0/html/chap13.html#VkMemoryRequirements
				// The memoryTypeBits member is identical for all VkBuffer objects created with the same value for the flags and usage
				// members in the VkBufferCreateInfo structure passed to vkCreateBuffer. Further, if usage1 and usage2 of type
				// VkBufferUsageFlags are such that the bits set in usage2 are a subset of the bits set in usage1, and they have the same flags,
				// then the bits set in memoryTypeBits returned for usage1 must be a subset of the bits set in memoryTypeBits returned for
				// usage2, for all values of flags.

				// asset vertex buffer
				vulkan_asset_vertex_buffer_memory_index = VK_MAX_MEMORY_TYPES;
				{
					VkDeviceSize memory_requirements_size = VkDeviceSize(-1);
					uint32_t memory_requirements_memory_type_bits = 0U;
					{
						struct VkBufferCreateInfo buffer_create_info;
						buffer_create_info.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
						buffer_create_info.pNext = NULL;
						buffer_create_info.flags = 0U;
						buffer_create_info.size = 8U;
						buffer_create_info.usage = VK_BUFFER_USAGE_TRANSFER_DST_BIT | VK_BUFFER_USAGE_VERTEX_BUFFER_BIT;
						buffer_create_info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
						buffer_create_info.queueFamilyIndexCount = 0U;
						buffer_create_info.pQueueFamilyIndices = NULL;

						VkBuffer dummy_buf;
						VkResult res_create_buffer = pfn_create_buffer(vulkan_device, &buffer_create_info, vulkan_allocation_callbacks, &dummy_buf);
						assert(VK_SUCCESS == res_create_buffer);

						struct VkMemoryRequirements memory_requirements;
						pfn_get_buffer_memory_requirements(vulkan_device, dummy_buf, &memory_requirements);
						memory_requirements_size = memory_requirements.size;
						memory_requirements_memory_type_bits = memory_requirements.memoryTypeBits;

						pfn_destroy_buffer(vulkan_device, dummy_buf, vulkan_allocation_callbacks);
					}

					// NOT HOST_VISIBLE
					// The UMA driver may compress the buffer/texture to boost performance
					vulkan_asset_vertex_buffer_memory_index = __internal_find_lowest_memory_type_index(&physical_device_memory_properties, memory_requirements_size, memory_requirements_memory_type_bits, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
					assert(VK_MAX_MEMORY_TYPES > vulkan_asset_vertex_buffer_memory_index);
					assert(physical_device_memory_properties.memoryTypeCount > vulkan_asset_vertex_buffer_memory_index);
				}

				// asset index buffer
				vulkan_asset_index_buffer_memory_index = VK_MAX_MEMORY_TYPES;
				{
					VkDeviceSize memory_requirements_size = VkDeviceSize(-1);
					uint32_t memory_requirements_memory_type_bits = 0U;
					{
						struct VkBufferCreateInfo buffer_create_info;
						buffer_create_info.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
						buffer_create_info.pNext = NULL;
						buffer_create_info.flags = 0U;
						buffer_create_info.size = 8U;
						buffer_create_info.usage = VK_BUFFER_USAGE_TRANSFER_DST_BIT | VK_BUFFER_USAGE_INDEX_BUFFER_BIT;
						buffer_create_info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
						buffer_create_info.queueFamilyIndexCount = 0U;
						buffer_create_info.pQueueFamilyIndices = NULL;

						VkBuffer dummy_buf;
						VkResult res_create_buffer = pfn_create_buffer(vulkan_device, &buffer_create_info, vulkan_allocation_callbacks, &dummy_buf);
						assert(VK_SUCCESS == res_create_buffer);

						struct VkMemoryRequirements memory_requirements;
						pfn_get_buffer_memory_requirements(vulkan_device, dummy_buf, &memory_requirements);
						memory_requirements_size = memory_requirements.size;
						memory_requirements_memory_type_bits = memory_requirements.memoryTypeBits;

						pfn_destroy_buffer(vulkan_device, dummy_buf, vulkan_allocation_callbacks);
					}

					// NOT HOST_VISIBLE
					// The UMA driver may compress the buffer/texture to boost performance
					vulkan_asset_index_buffer_memory_index = __internal_find_lowest_memory_type_index(&physical_device_memory_properties, memory_requirements_size, memory_requirements_memory_type_bits, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
					assert(VK_MAX_MEMORY_TYPES > vulkan_asset_index_buffer_memory_index);
					assert(physical_device_memory_properties.memoryTypeCount > vulkan_asset_index_buffer_memory_index);
				}

				// asset uniform buffer
				vulkan_asset_uniform_buffer_memory_index = VK_MAX_MEMORY_TYPES;
				{
					VkDeviceSize memory_requirements_size = VkDeviceSize(-1);
					uint32_t memory_requirements_memory_type_bits = 0U;
					{
						struct VkBufferCreateInfo buffer_create_info;
						buffer_create_info.sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
						buffer_create_info.pNext = NULL;
						buffer_create_info.flags = 0U;
						buffer_create_info.size = 8U;
						buffer_create_info.usage = VK_BUFFER_USAGE_TRANSFER_DST_BIT | VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT;
						buffer_create_info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
						buffer_create_info.queueFamilyIndexCount = 0U;
						buffer_create_info.pQueueFamilyIndices = NULL;

						VkBuffer dummy_buf;
						VkResult res_create_buffer = pfn_create_buffer(vulkan_device, &buffer_create_info, vulkan_allocation_callbacks, &dummy_buf);
						assert(VK_SUCCESS == res_create_buffer);

						struct VkMemoryRequirements memory_requirements;
						pfn_get_buffer_memory_requirements(vulkan_device, dummy_buf, &memory_requirements);
						memory_requirements_size = memory_requirements.size;
						memory_requirements_memory_type_bits = memory_requirements.memoryTypeBits;

						pfn_destroy_buffer(vulkan_device, dummy_buf, vulkan_allocation_callbacks);
					}

					// NOT HOST_VISIBLE
					// The UMA driver may compress the buffer/texture to boost performance
					vulkan_asset_uniform_buffer_memory_index = __internal_find_lowest_memory_type_index(&physical_device_memory_properties, memory_requirements_size, memory_requirements_memory_type_bits, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
					assert(VK_MAX_MEMORY_TYPES > vulkan_asset_uniform_buffer_memory_index);
					assert(physical_device_memory_properties.memoryTypeCount > vulkan_asset_uniform_buffer_memory_index);
				}

				// https://www.khronos.org/registry/vulkan/specs/1.0/html/chap13.html#VkMemoryRequirements
				// For images created with a color format, the memoryTypeBits member is identical for all VkImage objects created with the
				// same combination of values for the tiling member, the VK_IMAGE_CREATE_SPARSE_BINDING_BIT bit of the flags member, and
				// the VK_IMAGE_USAGE_TRANSIENT_ATTACHMENT_BIT of the usage member in the VkImageCreateInfo structure passed to
				// vkCreateImage.

				// vulkaninfo
				// https://github.com/KhronosGroup/Vulkan-Tools/tree/master/vulkaninfo/vulkaninfo/vulkaninfo.h
				// GetImageCreateInfo
				// FillImageTypeSupport
				// https://github.com/KhronosGroup/Vulkan-Tools/tree/master/vulkaninfo/vulkaninfo.cpp
				// GpuDumpMemoryProps //"usable for"

				vulkan_asset_image_memory_index = VK_MAX_MEMORY_TYPES;
				{
					VkDeviceSize memory_requirements_size = VkDeviceSize(-1);
					uint32_t memory_requirements_memory_type_bits = 0U;
					{
						VkFormat color_format = VK_FORMAT_R8G8B8A8_UNORM;

						struct VkFormatProperties format_properties;
						pfn_get_physical_device_format_properties(vulkan_physical_device, color_format, &format_properties);
						assert(format_properties.optimalTilingFeatures & VK_FORMAT_FEATURE_COLOR_ATTACHMENT_BIT);
						assert(format_properties.optimalTilingFeatures & VK_FORMAT_FEATURE_COLOR_ATTACHMENT_BLEND_BIT);
						assert(format_properties.optimalTilingFeatures & VK_FORMAT_FEATURE_SAMPLED_IMAGE_BIT);
						assert(format_properties.optimalTilingFeatures & VK_FORMAT_FEATURE_SAMPLED_IMAGE_FILTER_LINEAR_BIT);

						struct VkImageCreateInfo image_create_info_regular_tiling_optimal;
						image_create_info_regular_tiling_optimal.sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;
						image_create_info_regular_tiling_optimal.pNext = NULL;
						image_create_info_regular_tiling_optimal.flags = 0U;
						image_create_info_regular_tiling_optimal.imageType = VK_IMAGE_TYPE_2D;
						image_create_info_regular_tiling_optimal.format = color_format;
						image_create_info_regular_tiling_optimal.extent.width = 8U;
						image_create_info_regular_tiling_optimal.extent.height = 8U;
						image_create_info_regular_tiling_optimal.extent.depth = 1U;
						image_create_info_regular_tiling_optimal.mipLevels = 1U;
						image_create_info_regular_tiling_optimal.arrayLayers = 1U;
						image_create_info_regular_tiling_optimal.samples = VK_SAMPLE_COUNT_1_BIT;
						image_create_info_regular_tiling_optimal.tiling = VK_IMAGE_TILING_OPTIMAL;
						image_create_info_regular_tiling_optimal.usage = VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_SAMPLED_BIT;
						image_create_info_regular_tiling_optimal.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
						image_create_info_regular_tiling_optimal.queueFamilyIndexCount = 0U;
						image_create_info_regular_tiling_optimal.pQueueFamilyIndices = NULL;
						image_create_info_regular_tiling_optimal.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;

						VkImage dummy_img;
						VkResult res_create_image = pfn_create_image(vulkan_device, &image_create_info_regular_tiling_optimal, vulkan_allocation_callbacks, &dummy_img);
						assert(VK_SUCCESS == res_create_image);

						struct VkMemoryRequirements memory_requirements;
						pfn_get_image_memory_requirements(vulkan_device, dummy_img, &memory_requirements);
						memory_requirements_size = memory_requirements.size;
						memory_requirements_memory_type_bits = memory_requirements.memoryTypeBits;

						pfn_destroy_image(vulkan_device, dummy_img, vulkan_allocation_callbacks);
					}

					vulkan_asset_image_memory_index = __internal_find_lowest_memory_type_index(&physical_device_memory_properties, memory_requirements_size, memory_requirements_memory_type_bits, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
					assert(VK_MAX_MEMORY_TYPES > vulkan_asset_image_memory_index);
					assert(physical_device_memory_properties.memoryTypeCount > vulkan_asset_image_memory_index);
				}
			}

			if (vulkan_has_dedicated_transfer_queue)
			{
				if (vulkan_queue_graphics_family_index != vulkan_queue_transfer_family_index)
				{
					VkResult res_reset_transfer_command_pool = pfn_reset_command_pool(vulkan_device, vulkan_streaming_transfer_command_pool, 0U);
					assert(VK_SUCCESS == res_reset_transfer_command_pool);

					VkCommandBufferBeginInfo vulkan_transfer_command_buffer_begin_info = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO, NULL, VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT, NULL};
					VkResult res_begin_transfer_command_buffer = pfn_begin_command_buffer(vulkan_streaming_transfer_command_buffer, &vulkan_transfer_command_buffer_begin_info);
					assert(VK_SUCCESS == res_begin_transfer_command_buffer);

					VkResult res_reset_graphics_command_pool = pfn_reset_command_pool(vulkan_device, vulkan_streaming_graphics_command_pool, 0U);
					assert(VK_SUCCESS == res_reset_graphics_command_pool);

					VkCommandBufferBeginInfo vulkan_graphics_command_buffer_begin_info = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO, NULL, VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT, NULL};
					VkResult res_begin_graphics_command_buffer = pfn_begin_command_buffer(vulkan_streaming_graphics_command_buffer, &vulkan_graphics_command_buffer_begin_info);
					assert(VK_SUCCESS == res_begin_graphics_command_buffer);
				}
				else
				{
					VkResult res_reset_transfer_command_pool = pfn_reset_command_pool(vulkan_device, vulkan_streaming_transfer_command_pool, 0U);
					assert(VK_SUCCESS == res_reset_transfer_command_pool);

					VkCommandBufferBeginInfo vulkan_transfer_command_buffer_begin_info = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO, NULL, VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT, NULL};
					VkResult res_begin_transfer_command_buffer = pfn_begin_command_buffer(vulkan_streaming_transfer_command_buffer, &vulkan_transfer_command_buffer_begin_info);
					assert(VK_SUCCESS == res_begin_transfer_command_buffer);
				}
			}
			else
			{
				VkResult res_reset_graphics_command_pool = pfn_reset_command_pool(vulkan_device, vulkan_streaming_graphics_command_pool, 0U);
				assert(VK_SUCCESS == res_reset_graphics_command_pool);

				VkCommandBufferBeginInfo vulkan_graphics_command_buffer_begin_info = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO, NULL, VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT, NULL};
				VkResult res_begin_graphics_command_buffer = pfn_begin_command_buffer(vulkan_streaming_graphics_command_buffer, &vulkan_graphics_command_buffer_begin_info);
				assert(VK_SUCCESS == res_begin_graphics_command_buffer);
			}

			assert(vulkan_staging_buffer_size < static_cast<VkDeviceSize>(UINT32_MAX));
			demo.init(
				vulkan_instance, pfn_get_instance_proc_addr, vulkan_device, pfn_get_device_proc_addr, vulkan_allocation_callbacks,
				vulkan_depth_format, vulkan_swapchain_image_format,
				vulkan_asset_allocator, 0U, static_cast<uint32_t>(vulkan_staging_buffer_size), vulkan_staging_buffer_device_memory_pointer, vulkan_staging_buffer,
				vulkan_asset_vertex_buffer_memory_index, vulkan_asset_index_buffer_memory_index, vulkan_asset_uniform_buffer_memory_index, vulkan_asset_image_memory_index,
				vulkan_optimal_buffer_copy_offset_alignment, vulkan_optimal_buffer_copy_row_pitch_alignment,
				vulkan_has_dedicated_transfer_queue, vulkan_queue_transfer_family_index, vulkan_queue_graphics_family_index, vulkan_streaming_transfer_command_buffer, vulkan_streaming_graphics_command_buffer,
				vulkan_upload_ring_buffer_size, vulkan_upload_ring_buffer);

			if (vulkan_has_dedicated_transfer_queue)
			{
				if (vulkan_queue_graphics_family_index != vulkan_queue_transfer_family_index)
				{
					VkResult res_end_transfer_command_buffer = pfn_end_command_buffer(vulkan_streaming_transfer_command_buffer);
					assert(VK_SUCCESS == res_end_transfer_command_buffer);

					VkResult res_end_graphics_command_buffer = pfn_end_command_buffer(vulkan_streaming_graphics_command_buffer);
					assert(VK_SUCCESS == res_end_graphics_command_buffer);

					VkSubmitInfo submit_info_transfer;
					submit_info_transfer.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
					submit_info_transfer.pNext = NULL;
					submit_info_transfer.waitSemaphoreCount = 0U;
					submit_info_transfer.pWaitSemaphores = NULL;
					submit_info_transfer.pWaitDstStageMask = NULL;
					submit_info_transfer.commandBufferCount = 1U;
					submit_info_transfer.pCommandBuffers = &vulkan_streaming_transfer_command_buffer;
					submit_info_transfer.signalSemaphoreCount = 1U;
					submit_info_transfer.pSignalSemaphores = &vulkan_streaming_semaphore;
					VkResult res_transfer_queue_submit = pfn_queue_submit(vulkan_queue_transfer, 1U, &submit_info_transfer, VK_NULL_HANDLE);
					assert(VK_SUCCESS == res_transfer_queue_submit);

					VkSubmitInfo submit_info_graphics;
					submit_info_graphics.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
					submit_info_graphics.pNext = NULL;
					submit_info_graphics.waitSemaphoreCount = 1U;
					submit_info_graphics.pWaitSemaphores = &vulkan_streaming_semaphore;
					VkPipelineStageFlags wait_dst_stage_mask[1] = {VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT};
					submit_info_graphics.pWaitDstStageMask = wait_dst_stage_mask;
					submit_info_graphics.commandBufferCount = 1U;
					submit_info_graphics.pCommandBuffers = &vulkan_streaming_graphics_command_buffer;
					submit_info_graphics.signalSemaphoreCount = 0U;
					submit_info_graphics.pSignalSemaphores = NULL;
					// queue family ownership transfer - acquire operation
					VkResult res_graphics_queue_submit = pfn_queue_submit(vulkan_queue_graphics, 1U, &submit_info_graphics, vulkan_streaming_fence);
					assert(VK_SUCCESS == res_graphics_queue_submit);
				}
				else
				{
					VkResult res_end_transfer_command_buffer = pfn_end_command_buffer(vulkan_streaming_transfer_command_buffer);
					assert(VK_SUCCESS == res_end_transfer_command_buffer);

					VkSubmitInfo submit_info_transfer;
					submit_info_transfer.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
					submit_info_transfer.pNext = NULL;
					submit_info_transfer.waitSemaphoreCount = 0U;
					submit_info_transfer.pWaitSemaphores = NULL;
					submit_info_transfer.pWaitDstStageMask = NULL;
					submit_info_transfer.commandBufferCount = 1U;
					submit_info_transfer.pCommandBuffers = &vulkan_streaming_transfer_command_buffer;
					submit_info_transfer.signalSemaphoreCount = 0U;
					submit_info_transfer.pSignalSemaphores = NULL;
					VkResult res_transfer_queue_submit = pfn_queue_submit(vulkan_queue_transfer, 1U, &submit_info_transfer, vulkan_streaming_fence);
					assert(VK_SUCCESS == res_transfer_queue_submit);
				}
			}
			else
			{
				VkResult res_end_graphics_command_buffer = pfn_end_command_buffer(vulkan_streaming_graphics_command_buffer);
				assert(VK_SUCCESS == res_end_graphics_command_buffer);

				VkSubmitInfo submit_info_graphics;
				submit_info_graphics.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
				submit_info_graphics.pNext = NULL;
				submit_info_graphics.waitSemaphoreCount = 0U;
				submit_info_graphics.pWaitSemaphores = NULL;
				submit_info_graphics.pWaitDstStageMask = NULL;
				submit_info_graphics.commandBufferCount = 1U;
				submit_info_graphics.pCommandBuffers = &vulkan_streaming_graphics_command_buffer;
				submit_info_graphics.signalSemaphoreCount = 0U;
				submit_info_graphics.pSignalSemaphores = NULL;
				VkResult res_graphics_queue_submit = pfn_queue_submit(vulkan_queue_graphics, 1U, &submit_info_graphics, vulkan_streaming_fence);
				assert(VK_SUCCESS == res_graphics_queue_submit);
			}

			VkResult res_wait_for_fences = pfn_wait_for_fences(vulkan_device, 1U, &vulkan_streaming_fence, VK_TRUE, UINT64_MAX);
			assert(VK_SUCCESS == res_wait_for_fences);

			// free staging buffer
			{
				PFN_vkUnmapMemory pfn_unmap_memory = reinterpret_cast<PFN_vkUnmapMemory>(pfn_get_device_proc_addr(vulkan_device, "vkUnmapMemory"));
				PFN_vkDestroyBuffer pfn_destroy_buffer = reinterpret_cast<PFN_vkDestroyBuffer>(pfn_get_device_proc_addr(vulkan_device, "vkDestroyBuffer"));
				PFN_vkFreeMemory pfn_free_memory = reinterpret_cast<PFN_vkFreeMemory>(pfn_get_device_proc_addr(vulkan_device, "vkFreeMemory"));

				pfn_unmap_memory(vulkan_device, vulkan_staging_buffer_device_memory);
				pfn_destroy_buffer(vulkan_device, vulkan_staging_buffer, vulkan_allocation_callbacks);
				pfn_free_memory(vulkan_device, vulkan_staging_buffer_device_memory, vulkan_allocation_callbacks);
			}

			// TODO:
			{
				// destory fence/semaphore/commandbuffer/commandpool
			}
		}

		// FrameBuffer
		uint32_t vulkan_framebuffer_width = 0U;
		uint32_t vulkan_framebuffer_height = 0U;
		VkSwapchainKHR vulkan_swapchain = VK_NULL_HANDLE;
		uint32_t vulkan_swapchain_image_count = 0U;
		std::vector<VkImage> vulkan_swapchain_images;
		std::vector<VkImageView> vulkan_swapchain_image_views;

		// Rendering
		uint32_t frame_throtting_index = 0U;
		while (!g_window_quit)
		{
			if (0U == vulkan_framebuffer_width || 0U == vulkan_framebuffer_height)
			{
				vulkan_query_frame_buffer_extent(
					vulkan_instance, pfn_get_instance_proc_addr, vulkan_physical_device,
					vulkan_surface,
					vulkan_framebuffer_width, vulkan_framebuffer_height);

				if (0U != vulkan_framebuffer_width && 0U != vulkan_framebuffer_height)
				{
					vulkan_create_swapchain(
						vulkan_device, pfn_get_device_proc_addr, vulkan_allocation_callbacks,
						vulkan_surface, vulkan_swapchain_image_format,
						vulkan_framebuffer_width, vulkan_framebuffer_height,
						vulkan_swapchain,
						vulkan_swapchain_image_count,
						vulkan_swapchain_image_views);

					demo.create_frame_buffer(
						vulkan_instance, pfn_get_instance_proc_addr, vulkan_physical_device, vulkan_device, pfn_get_device_proc_addr, vulkan_allocation_callbacks,
						vulkan_color_input_transient_attachment_memory_index, vulkan_depth_format, vulkan_depth_stencil_transient_attachment_memory_index,
						vulkan_framebuffer_width, vulkan_framebuffer_height,
						vulkan_swapchain,
						vulkan_swapchain_image_count,
						vulkan_swapchain_image_views);
				}

				// skip current frame
				continue;
			}

			VkCommandPool vulkan_command_pool = vulkan_command_pools[frame_throtting_index];
			VkCommandBuffer vulkan_command_buffer = vulkan_command_buffers[frame_throtting_index];
			VkFence vulkan_fence = vulkan_fences[frame_throtting_index];
			VkSemaphore vulkan_semaphore_acquire_next_image = vulkan_semaphores_acquire_next_image[frame_throtting_index];
			VkSemaphore vulkan_semaphore_queue_submit = vulkan_semaphores_queue_submit[frame_throtting_index];

			uint32_t vulkan_swapchain_image_index = uint32_t(-1);
			VkResult res_acquire_next_image = pfn_acquire_next_image(vulkan_device, vulkan_swapchain, UINT64_MAX, vulkan_semaphore_acquire_next_image, VK_NULL_HANDLE, &vulkan_swapchain_image_index);
			if (VK_SUCCESS != res_acquire_next_image)
			{
				assert(VK_SUBOPTIMAL_KHR == res_acquire_next_image || VK_ERROR_OUT_OF_DATE_KHR == res_acquire_next_image);

				VkResult res_wait_for_fences = pfn_wait_for_fences(vulkan_device, FRAME_THROTTLING_COUNT, vulkan_fences, VK_TRUE, UINT64_MAX);
				assert(VK_SUCCESS == res_wait_for_fences);

				demo.destroy_frame_buffer(
					vulkan_device, pfn_get_device_proc_addr, vulkan_allocation_callbacks,
					vulkan_swapchain_image_count);

				vulkan_destroy_swapchain(
					vulkan_device, pfn_get_device_proc_addr, vulkan_allocation_callbacks,
					vulkan_framebuffer_width, vulkan_framebuffer_height,
					vulkan_swapchain,
					vulkan_swapchain_image_count,
					vulkan_swapchain_image_views);

				vulkan_query_frame_buffer_extent(
					vulkan_instance, pfn_get_instance_proc_addr, vulkan_physical_device,
					vulkan_surface,
					vulkan_framebuffer_width, vulkan_framebuffer_height);

				if (0U != vulkan_framebuffer_width && 0U != vulkan_framebuffer_height)
				{

					vulkan_create_swapchain(
						vulkan_device, pfn_get_device_proc_addr, vulkan_allocation_callbacks,
						vulkan_surface, vulkan_swapchain_image_format,
						vulkan_framebuffer_width, vulkan_framebuffer_height,
						vulkan_swapchain,
						vulkan_swapchain_image_count,
						vulkan_swapchain_image_views);

					demo.create_frame_buffer(
						vulkan_instance, pfn_get_instance_proc_addr, vulkan_physical_device, vulkan_device, pfn_get_device_proc_addr, vulkan_allocation_callbacks,
						vulkan_color_input_transient_attachment_memory_index, vulkan_depth_format, vulkan_depth_stencil_transient_attachment_memory_index,
						vulkan_framebuffer_width, vulkan_framebuffer_height,
						vulkan_swapchain,
						vulkan_swapchain_image_count,
						vulkan_swapchain_image_views);
				}

				// submission (without any commandbuffers) to unsignal semaphore
				VkSubmitInfo submit_info;
				submit_info.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
				submit_info.pNext = NULL;
				submit_info.waitSemaphoreCount = 1U;
				submit_info.pWaitSemaphores = &vulkan_semaphore_acquire_next_image;
				VkPipelineStageFlags wait_dst_stage_mask = VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT;
				submit_info.pWaitDstStageMask = &wait_dst_stage_mask;
				submit_info.commandBufferCount = 0U;
				submit_info.pCommandBuffers = NULL;
				submit_info.signalSemaphoreCount = 0U;
				submit_info.pSignalSemaphores = NULL;
				VkResult res_queue_submit = pfn_queue_submit(vulkan_queue_graphics, 1U, &submit_info, VK_NULL_HANDLE);
				assert(VK_SUCCESS == res_queue_submit);

				// skip current frame
				continue;
			}

			VkResult res_wait_for_fences = pfn_wait_for_fences(vulkan_device, 1U, &vulkan_fence, VK_TRUE, UINT64_MAX);
			assert(VK_SUCCESS == res_wait_for_fences);

			VkResult res_reset_command_pool = pfn_reset_command_pool(vulkan_device, vulkan_command_pool, 0U);
			assert(VK_SUCCESS == res_reset_command_pool);

			VkCommandBufferBeginInfo vulkan_command_buffer_begin_info = {VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO, NULL, VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT, NULL};
			VkResult res_begin_command_buffer = pfn_begin_command_buffer(vulkan_command_buffer, &vulkan_command_buffer_begin_info);
			assert(VK_SUCCESS == res_begin_command_buffer);

			vulkan_upload_ring_buffer_current[frame_throtting_index] = vulkan_upload_ring_buffer_begin[frame_throtting_index];

			demo.tick(vulkan_command_buffer, vulkan_swapchain_image_index, vulkan_framebuffer_width, vulkan_framebuffer_height,
					  vulkan_upload_ring_buffer_device_memory_pointer, vulkan_upload_ring_buffer_current[frame_throtting_index], vulkan_upload_ring_buffer_end[frame_throtting_index],
					  vulkan_min_uniform_buffer_offset_alignment);

			VkResult res_end_command_buffer = pfn_end_command_buffer(vulkan_command_buffer);
			assert(VK_SUCCESS == res_end_command_buffer);

			VkResult res_reset_fences = pfn_reset_fences(vulkan_device, 1U, &vulkan_fence);
			assert(VK_SUCCESS == res_reset_fences);

			VkSubmitInfo submit_info;
			submit_info.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
			submit_info.pNext = NULL;
			submit_info.waitSemaphoreCount = 1U;
			submit_info.pWaitSemaphores = &vulkan_semaphore_acquire_next_image;
			VkPipelineStageFlags wait_dst_stage_mask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
			submit_info.pWaitDstStageMask = &wait_dst_stage_mask;
			submit_info.commandBufferCount = 1U;
			submit_info.pCommandBuffers = &vulkan_command_buffer;
			submit_info.signalSemaphoreCount = 1U;
			submit_info.pSignalSemaphores = &vulkan_semaphore_queue_submit;
			VkResult res_queue_submit = pfn_queue_submit(vulkan_queue_graphics, 1U, &submit_info, vulkan_fence);
			assert(VK_SUCCESS == res_queue_submit);

			VkPresentInfoKHR present_info;
			present_info.sType = VK_STRUCTURE_TYPE_PRESENT_INFO_KHR;
			present_info.pNext = NULL;
			present_info.waitSemaphoreCount = 1U;
			present_info.pWaitSemaphores = &vulkan_semaphore_queue_submit;
			present_info.swapchainCount = 1U;
			present_info.pSwapchains = &vulkan_swapchain;
			present_info.pImageIndices = &vulkan_swapchain_image_index;
			present_info.pResults = NULL;
			VkResult res_queue_present = pfn_queue_present(vulkan_queue_graphics, &present_info);
			if (VK_SUCCESS != res_acquire_next_image)
			{
				assert(VK_SUBOPTIMAL_KHR == res_acquire_next_image || VK_ERROR_OUT_OF_DATE_KHR == res_acquire_next_image);

				demo.destroy_frame_buffer(
					vulkan_device, pfn_get_device_proc_addr, vulkan_allocation_callbacks,
					vulkan_swapchain_image_count);

				vulkan_destroy_swapchain(
					vulkan_device, pfn_get_device_proc_addr, vulkan_allocation_callbacks,
					vulkan_framebuffer_width, vulkan_framebuffer_height,
					vulkan_swapchain,
					vulkan_swapchain_image_count,
					vulkan_swapchain_image_views);

				vulkan_query_frame_buffer_extent(
					vulkan_instance, pfn_get_instance_proc_addr, vulkan_physical_device,
					vulkan_surface,
					vulkan_framebuffer_width, vulkan_framebuffer_height);

				if (0U != vulkan_framebuffer_width && 0U != vulkan_framebuffer_height)
				{
					vulkan_create_swapchain(
						vulkan_device, pfn_get_device_proc_addr, vulkan_allocation_callbacks,
						vulkan_surface, vulkan_swapchain_image_format,
						vulkan_framebuffer_width, vulkan_framebuffer_height,
						vulkan_swapchain,
						vulkan_swapchain_image_count,
						vulkan_swapchain_image_views);

					demo.create_frame_buffer(
						vulkan_instance, pfn_get_instance_proc_addr, vulkan_physical_device, vulkan_device, pfn_get_device_proc_addr, vulkan_allocation_callbacks,
						vulkan_color_input_transient_attachment_memory_index, vulkan_depth_format, vulkan_depth_stencil_transient_attachment_memory_index,
						vulkan_framebuffer_width, vulkan_framebuffer_height,
						vulkan_swapchain,
						vulkan_swapchain_image_count,
						vulkan_swapchain_image_views);
				}
			}

			++frame_throtting_index;
			frame_throtting_index %= FRAME_THROTTLING_COUNT;
		}

		// Destroy
		{
			VkResult res_wait_for_fences = pfn_wait_for_fences(vulkan_device, FRAME_THROTTLING_COUNT, vulkan_fences, VK_TRUE, UINT64_MAX);
			assert(VK_SUCCESS == res_wait_for_fences);

			demo.destroy(vulkan_device, pfn_get_device_proc_addr, vulkan_allocation_callbacks, &vulkan_asset_allocator);
		}
	}

	// TODO:
	// destory

	return 0U;
}

static VkBool32 VKAPI_CALL vulkan_debug_utils_messenger_callback(VkDebugUtilsMessageSeverityFlagBitsEXT messageSeverity, VkDebugUtilsMessageTypeFlagsEXT messageTypes, const VkDebugUtilsMessengerCallbackDataEXT *pCallbackData, void *pUserData)
{
	OutputDebugStringA(pCallbackData->pMessage);
	OutputDebugStringA("\n");
	return VK_FALSE;
}

static inline uint32_t __internal_find_lowest_memory_type_index(struct VkPhysicalDeviceMemoryProperties const *physical_device_memory_properties, VkDeviceSize memory_requirements_size, uint32_t memory_requirements_memory_type_bits, VkMemoryPropertyFlags required_property_flags)
{
	uint32_t memory_type_count = physical_device_memory_properties->memoryTypeCount;
	assert(VK_MAX_MEMORY_TYPES >= memory_type_count);

	// The lower memory_type_index indicates the more performance
	for (uint32_t memory_type_index = 0; memory_type_index < memory_type_count; ++memory_type_index)
	{
		uint32_t memory_type_bits = (1U << memory_type_index);
		bool is_required_memory_type = ((memory_requirements_memory_type_bits & memory_type_bits) != 0) ? true : false;

		VkMemoryPropertyFlags property_flags = physical_device_memory_properties->memoryTypes[memory_type_index].propertyFlags;
		bool has_required_property_flags = ((property_flags & required_property_flags) == required_property_flags) ? true : false;

		uint32_t heap_index = physical_device_memory_properties->memoryTypes[memory_type_index].heapIndex;
		VkDeviceSize heap_budget = physical_device_memory_properties->memoryHeaps[heap_index].size;
		// The application is not alone and there may be other applications which interact with the Vulkan as well.
		// The allocation may success even if the budget has been exceeded. However, this may result in performance issue.
		bool is_within_budget = (memory_requirements_size <= heap_budget) ? true : false;

		if (is_required_memory_type && has_required_property_flags && is_within_budget)
		{
			return memory_type_index;
		}
	}

	return VK_MAX_MEMORY_TYPES;
}

static inline uint32_t __internal_find_lowest_memory_type_index(struct VkPhysicalDeviceMemoryProperties const *physical_device_memory_properties, VkDeviceSize memory_requirements_size, uint32_t memory_requirements_memory_type_bits, VkMemoryPropertyFlags required_property_flags, VkMemoryPropertyFlags preferred_property_flags)
{
	VkMemoryPropertyFlags optimal_property_flags = (required_property_flags | preferred_property_flags);
	uint32_t memory_type_index = __internal_find_lowest_memory_type_index(physical_device_memory_properties, memory_requirements_size, memory_requirements_memory_type_bits, optimal_property_flags);
	if (VK_MAX_MEMORY_TYPES != memory_type_index)
	{
		assert(VK_MAX_MEMORY_TYPES > memory_type_index);
		return memory_type_index;
	}
	else
	{
		return __internal_find_lowest_memory_type_index(physical_device_memory_properties, memory_requirements_size, memory_requirements_memory_type_bits, required_property_flags);
	}
}

static inline void vulkan_query_frame_buffer_extent(
	VkInstance vulkan_instance, PFN_vkGetInstanceProcAddr pfn_get_instance_proc_addr, VkPhysicalDevice vulkan_physical_device,
	VkSurfaceKHR vulkan_surface,
	uint32_t &vulkan_framebuffer_width, uint32_t &vulkan_framebuffer_height)
{
	assert(0U == vulkan_framebuffer_width);
	assert(0U == vulkan_framebuffer_height);
	{
		PFN_vkGetPhysicalDeviceSurfaceCapabilitiesKHR pfn_get_physical_device_surface_capabilities = reinterpret_cast<PFN_vkGetPhysicalDeviceSurfaceCapabilitiesKHR>(pfn_get_instance_proc_addr(vulkan_instance, "vkGetPhysicalDeviceSurfaceCapabilitiesKHR"));
		assert(NULL != pfn_get_physical_device_surface_capabilities);

		VkSurfaceCapabilitiesKHR surface_capabilities;
		VkResult res_get_physical_device_surface_capablilities = pfn_get_physical_device_surface_capabilities(vulkan_physical_device, vulkan_surface, &surface_capabilities);
		assert(VK_SUCCESS == res_get_physical_device_surface_capablilities);

		if (surface_capabilities.currentExtent.width != 0XFFFFFFFFU)
		{
			vulkan_framebuffer_width = surface_capabilities.currentExtent.width;
		}
		else
		{
			vulkan_framebuffer_width = g_resolution_width;
		}
		vulkan_framebuffer_width = (vulkan_framebuffer_width < surface_capabilities.minImageExtent.width) ? surface_capabilities.minImageExtent.width : (vulkan_framebuffer_width > surface_capabilities.maxImageExtent.width) ? surface_capabilities.maxImageExtent.width
																																																							   : vulkan_framebuffer_width;

		if (surface_capabilities.currentExtent.height != 0XFFFFFFFFU)
		{
			vulkan_framebuffer_height = surface_capabilities.currentExtent.height;
		}
		else
		{
			vulkan_framebuffer_height = g_resolution_height;
		}
		vulkan_framebuffer_height = (vulkan_framebuffer_height < surface_capabilities.minImageExtent.height) ? surface_capabilities.minImageExtent.height : (vulkan_framebuffer_height > surface_capabilities.maxImageExtent.height) ? surface_capabilities.maxImageExtent.height
																																																									 : vulkan_framebuffer_height;
		assert(VK_SURFACE_TRANSFORM_IDENTITY_BIT_KHR == (VK_SURFACE_TRANSFORM_IDENTITY_BIT_KHR & surface_capabilities.supportedTransforms));
		assert(VK_COMPOSITE_ALPHA_OPAQUE_BIT_KHR == (VK_COMPOSITE_ALPHA_OPAQUE_BIT_KHR & surface_capabilities.supportedCompositeAlpha));
	}
}

static inline void vulkan_create_swapchain(
	VkDevice vulkan_device, PFN_vkGetDeviceProcAddr pfn_get_device_proc_addr, VkAllocationCallbacks *vulkan_allocation_callbacks,
	VkSurfaceKHR vulkan_surface, VkFormat vulkan_swapchain_image_format,
	uint32_t vulkan_framebuffer_width, uint32_t vulkan_framebuffer_height,
	VkSwapchainKHR &vulkan_swapchain,
	uint32_t &vulkan_swapchain_image_count,
	std::vector<VkImageView> &vulkan_swapchain_image_views)
{
	assert(0U != vulkan_framebuffer_width);
	assert(0U != vulkan_framebuffer_height);

	assert(VK_NULL_HANDLE == vulkan_swapchain);
	{
		PFN_vkCreateSwapchainKHR pfn_create_swapchain = reinterpret_cast<PFN_vkCreateSwapchainKHR>(pfn_get_device_proc_addr(vulkan_device, "vkCreateSwapchainKHR"));
		assert(NULL != pfn_create_swapchain);

		VkSwapchainCreateInfoKHR swapchain_create_info;
		swapchain_create_info.sType = VK_STRUCTURE_TYPE_SWAPCHAIN_CREATE_INFO_KHR;
		swapchain_create_info.pNext = NULL;
		swapchain_create_info.flags = 0U;
		swapchain_create_info.surface = vulkan_surface;
		swapchain_create_info.minImageCount = FRAME_THROTTLING_COUNT;
		swapchain_create_info.imageFormat = vulkan_swapchain_image_format;
		swapchain_create_info.imageColorSpace = VK_COLOR_SPACE_SRGB_NONLINEAR_KHR;
		swapchain_create_info.imageExtent.width = vulkan_framebuffer_width;
		swapchain_create_info.imageExtent.height = vulkan_framebuffer_height;
		swapchain_create_info.imageArrayLayers = 1U;
		swapchain_create_info.imageUsage = VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT;
		swapchain_create_info.imageSharingMode = VK_SHARING_MODE_EXCLUSIVE;
		swapchain_create_info.queueFamilyIndexCount = 0U;
		swapchain_create_info.pQueueFamilyIndices = NULL;
		swapchain_create_info.preTransform = VK_SURFACE_TRANSFORM_IDENTITY_BIT_KHR;
		swapchain_create_info.compositeAlpha = VK_COMPOSITE_ALPHA_OPAQUE_BIT_KHR;
		swapchain_create_info.presentMode = VK_PRESENT_MODE_FIFO_KHR;
		swapchain_create_info.clipped = VK_FALSE;
		swapchain_create_info.oldSwapchain = VK_NULL_HANDLE;

		VkResult res_create_swapchain = pfn_create_swapchain(vulkan_device, &swapchain_create_info, vulkan_allocation_callbacks, &vulkan_swapchain);
		assert(VK_SUCCESS == res_create_swapchain);
	}

	assert(0U == vulkan_swapchain_image_count);
	std::vector<VkImage> vulkan_swapchain_images;
	{
		PFN_vkGetSwapchainImagesKHR pfn_get_swapchain_images = reinterpret_cast<PFN_vkGetSwapchainImagesKHR>(pfn_get_device_proc_addr(vulkan_device, "vkGetSwapchainImagesKHR"));
		assert(NULL != pfn_get_swapchain_images);

		uint32_t swapchain_image_count_1 = uint32_t(-1);
		VkResult res_get_swapchain_images_1 = pfn_get_swapchain_images(vulkan_device, vulkan_swapchain, &swapchain_image_count_1, NULL);
		assert(VK_SUCCESS == res_get_swapchain_images_1);

		VkImage *swapchain_images = static_cast<VkImage *>(malloc(sizeof(VkImage) * swapchain_image_count_1));
		assert(NULL != swapchain_images);

		uint32_t swapchain_image_count_2 = swapchain_image_count_1;
		VkResult res_get_swapchain_images_2 = pfn_get_swapchain_images(vulkan_device, vulkan_swapchain, &swapchain_image_count_2, swapchain_images);
		assert(VK_SUCCESS == res_get_swapchain_images_2 && swapchain_image_count_2 == swapchain_image_count_1);

		vulkan_swapchain_image_count = swapchain_image_count_2;
		vulkan_swapchain_images.assign(swapchain_images, swapchain_images + swapchain_image_count_2);

		free(swapchain_images);
	}

	assert(0U == vulkan_swapchain_image_views.size());
	{
		vulkan_swapchain_image_views.resize(vulkan_swapchain_image_count);

		PFN_vkCreateImageView pfn_create_image_view = reinterpret_cast<PFN_vkCreateImageView>(pfn_get_device_proc_addr(vulkan_device, "vkCreateImageView"));
		assert(NULL != pfn_create_image_view);

		for (uint32_t vulkan_swapchain_image_index = 0U; vulkan_swapchain_image_index < vulkan_swapchain_image_count; ++vulkan_swapchain_image_index)
		{
			VkImageViewCreateInfo image_view_create_info = {
				VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO,
				NULL,
				0U,
				vulkan_swapchain_images[vulkan_swapchain_image_index],
				VK_IMAGE_VIEW_TYPE_2D,
				vulkan_swapchain_image_format,
				{VK_COMPONENT_SWIZZLE_R, VK_COMPONENT_SWIZZLE_G, VK_COMPONENT_SWIZZLE_B, VK_COMPONENT_SWIZZLE_A},
				{VK_IMAGE_ASPECT_COLOR_BIT, 0U, 1U, 0U, 1U}};

			VkResult res_create_image_view = pfn_create_image_view(vulkan_device, &image_view_create_info, vulkan_allocation_callbacks, &vulkan_swapchain_image_views[vulkan_swapchain_image_index]);
			assert(VK_SUCCESS == res_create_image_view);
		}
	}
}

static inline void vulkan_destroy_swapchain(
	VkDevice vulkan_device, PFN_vkGetDeviceProcAddr pfn_get_device_proc_addr, VkAllocationCallbacks *vulkan_allocation_callbacks,
	uint32_t &vulkan_framebuffer_width, uint32_t &vulkan_framebuffer_height,
	VkSwapchainKHR &vulkan_swapchain,
	uint32_t &vulkan_swapchain_image_count,
	std::vector<VkImageView> &vulkan_swapchain_image_views)
{
	// Swapchain ImageView
	{
		PFN_vkDestroyImageView pfn_destroy_image_view = reinterpret_cast<PFN_vkDestroyImageView>(pfn_get_device_proc_addr(vulkan_device, "vkDestroyImageView"));

		assert(static_cast<uint32_t>(vulkan_swapchain_image_views.size()) == vulkan_swapchain_image_count);

		for (uint32_t vulkan_swapchain_image_index = 0U; vulkan_swapchain_image_index < vulkan_swapchain_image_count; ++vulkan_swapchain_image_index)
		{
			pfn_destroy_image_view(vulkan_device, vulkan_swapchain_image_views[vulkan_swapchain_image_index], vulkan_allocation_callbacks);
		}
		vulkan_swapchain_image_views.clear();
	}

	// SwapChain
	{
		PFN_vkDestroySwapchainKHR pfn_destroy_swapchain = reinterpret_cast<PFN_vkDestroySwapchainKHR>(pfn_get_device_proc_addr(vulkan_device, "vkDestroySwapchainKHR"));
		pfn_destroy_swapchain(vulkan_device, vulkan_swapchain, vulkan_allocation_callbacks);
		vulkan_swapchain = VK_NULL_HANDLE;
	}

	// Count and Extent
	{
		vulkan_swapchain_image_count = 0U;
		vulkan_framebuffer_width = 0U;
		vulkan_framebuffer_height = 0U;
	}
}