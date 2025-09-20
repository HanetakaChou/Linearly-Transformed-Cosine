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

#ifndef _DEMO_H_
#define _DEMO_H_ 1

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN 1
#endif
#ifndef NOMINMAX
#define NOMINMAX 1
#endif
#ifndef VK_USE_PLATFORM_WIN32_KHR
#define VK_USE_PLATFORM_WIN32_KHR 1
#endif
#include <sdkddkver.h>
#include <windows.h>
#include "../thirdparty/Vulkan-Headers/include/vulkan/vulkan.h"
#include "../thirdparty/VulkanMemoryAllocator/include/vk_mem_alloc.h"
#include <vector>

class Demo
{
	VkRenderPass m_vulkan_render_pass;
	VkImage m_vulkan_depth_image;
	VkDeviceMemory m_vulkan_depth_device_memory;
	VkImageView m_vulkan_depth_image_view;
	VkImage m_vulkan_backup_image;
	VkDeviceMemory m_vulkan_backup_device_memory;
	VkImageView m_vulkan_backup_image_view;
	std::vector<VkFramebuffer> m_vulkan_framebuffers;
	VkDescriptorSet m_vulkan_framebuffer_set;

	VkSampler m_vulkan_linear_clamp_sampler;

	VkDescriptorSetLayout m_vulkan_global_set_layout;
	VkDescriptorSetLayout m_vulkan_material_set_layout;
	VkDescriptorSetLayout m_vulkan_framebuffer_set_layout;
	VkPipelineLayout m_vulkan_pipeline_layout;
	VkPipeline m_vulkan_forward_shading_pipeline;
	VkPipeline m_vulkan_post_process_pipeline;

	VkDescriptorPool m_vulkan_descriptor_pool;

	VmaAllocation m_ltc_matrix_lut_allocation;
	VkImage m_ltc_matrix_lut;
	VkImageView m_ltc_matrix_lut_view;
	VmaAllocation m_ltc_norm_lut_allocation;
	VkImage m_ltc_norm_lut;
	VkImageView m_ltc_norm_lut_view;
	VkDescriptorSet m_vulkan_global_set;

	uint32_t m_plane_vertex_count;
	VmaAllocation m_plane_vertex_position_allocation;
	VkBuffer m_plane_vertex_position_buffer;
	VmaAllocation m_plane_vertex_varying_allocation;
	VkBuffer m_plane_vertex_varying_buffer;
	VmaAllocation m_plane_material_buffer_allocation;
	VkBuffer m_plane_material_buffer;
	VkDescriptorSet m_plane_material_set;

	PFN_vkCmdBeginRenderPass m_pfn_cmd_begin_render_pass;
	PFN_vkCmdNextSubpass m_pfn_cmd_next_subpass;
	PFN_vkCmdEndRenderPass m_pfn_cmd_end_render_pass;
	PFN_vkCmdBindPipeline m_pfn_cmd_bind_pipeline;
	PFN_vkCmdSetViewport m_pfn_cmd_set_viewport;
	PFN_vkCmdSetScissor m_pfn_cmd_set_scissor;
	PFN_vkCmdBindDescriptorSets m_pfn_cmd_bind_descriptor_sets;
	PFN_vkCmdBindVertexBuffers m_pfn_cmd_bind_vertex_buffers;
	PFN_vkCmdDraw m_pfn_cmd_draw;
#ifndef NDEBUG
	PFN_vkCmdBeginDebugUtilsLabelEXT m_pfn_cmd_begin_debug_utils_label;
	PFN_vkCmdEndDebugUtilsLabelEXT m_pfn_cmd_end_debug_utils_label;
#endif

public:
	Demo();

	void init(
		VkInstance vulkan_instance, PFN_vkGetInstanceProcAddr pfn_get_instance_proc_addr, VkDevice vulkan_device, PFN_vkGetDeviceProcAddr pfn_get_device_proc_addr, VkAllocationCallbacks *vulkan_allocation_callbacks,
		VkFormat vulkan_depth_format, VkFormat vulkan_swapchain_image_format,
		VmaAllocator vulkan_asset_allocator, uint32_t staging_buffer_current, uint32_t vulkan_staging_buffer_end, void *vulkan_staging_buffer_device_memory_pointer, VkBuffer vulkan_staging_buffer,
		uint32_t vulkan_asset_vertex_buffer_memory_index, uint32_t vulkan_asset_index_buffer_memory_index, uint32_t vulkan_asset_uniform_buffer_memory_index, uint32_t vulkan_asset_image_memory_index,
		uint32_t vulkan_optimal_buffer_copy_offset_alignment, uint32_t vulkan_optimal_buffer_copy_row_pitch_alignment,
		bool vulkan_has_dedicated_transfer_queue, uint32_t vulkan_queue_transfer_family_index, uint32_t vulkan_queue_graphics_family_index, VkCommandBuffer vulkan_streaming_transfer_command_buffer, VkCommandBuffer vulkan_streaming_graphics_command_buffer,
		VkDeviceSize vulkan_upload_ring_buffer_size, VkBuffer vulkan_upload_ring_buffer);

	void create_frame_buffer(
		VkInstance vulkan_instance, PFN_vkGetInstanceProcAddr pfn_get_instance_proc_addr, VkPhysicalDevice vulkan_physical_device, VkDevice vulkan_device, PFN_vkGetDeviceProcAddr pfn_get_device_proc_addr, VkAllocationCallbacks *vulkan_allocation_callbacks,
		uint32_t vulkan_color_input_transient_attachment_memory_index, VkFormat vulkan_depth_format, uint32_t vulkan_depth_stencil_transient_attachment_memory_index,
		uint32_t vulkan_framebuffer_width, uint32_t vulkan_framebuffer_height,
		VkSwapchainKHR vulkan_swapchain,
		uint32_t vulkan_swapchain_image_count,
		std::vector<VkImageView> const &vulkan_swapchain_image_views);

	void destroy_frame_buffer(
		VkDevice vulkan_device, PFN_vkGetDeviceProcAddr pfn_get_device_proc_addr, VkAllocationCallbacks *vulkan_allocation_callbacks, uint32_t vulkan_swapchain_image_count);

	void tick(
		VkCommandBuffer vulkan_command_buffer, uint32_t vulkan_swapchain_image_index, uint32_t vulkan_framebuffer_width, uint32_t vulkan_framebuffer_height,
		void *vulkan_upload_ring_buffer_device_memory_pointer, uint32_t vulkan_upload_ring_buffer_current, uint32_t vulkan_upload_ring_buffer_end,
		uint32_t vulkan_min_uniform_buffer_offset_alignment);

	void destroy(VkDevice vulkan_device, PFN_vkGetDeviceProcAddr pfn_get_device_proc_addr, VkAllocationCallbacks *vulkan_allocation_callbacks, VmaAllocator *vulkan_asset_allocator);
};

#endif