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

#include "streaming.h"
#include <assert.h>

void streaming_staging_to_asset_buffer(
	PFN_vkCmdPipelineBarrier pfn_cmd_pipeline_barrier, PFN_vkCmdCopyBuffer pfn_cmd_copy_buffer,
	bool vulkan_has_dedicated_transfer_queue, uint32_t vulkan_queue_transfer_family_index, uint32_t vulkan_queue_graphics_family_index, VkCommandBuffer vulkan_streaming_transfer_command_buffer, VkCommandBuffer vulkan_streaming_graphics_command_buffer,
	VkBuffer staging_buffer, VkBuffer asset_buffer, uint32_t region_count, VkBufferCopy *const regions, enum STREAMING_ASSET_BUFFER_TYPE asset_buffer_type)
{
	if (vulkan_has_dedicated_transfer_queue)
	{
		if (vulkan_queue_graphics_family_index != vulkan_queue_transfer_family_index)
		{
			// barrier undefine - transfer_dst
			VkBufferMemoryBarrier buffer_memory_barrier_undefine_to_transfer_dst[1];
			buffer_memory_barrier_undefine_to_transfer_dst[0].sType = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER;
			buffer_memory_barrier_undefine_to_transfer_dst[0].pNext = NULL;
			buffer_memory_barrier_undefine_to_transfer_dst[0].srcAccessMask = 0U;
			buffer_memory_barrier_undefine_to_transfer_dst[0].dstAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
			buffer_memory_barrier_undefine_to_transfer_dst[0].srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
			buffer_memory_barrier_undefine_to_transfer_dst[0].dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
			buffer_memory_barrier_undefine_to_transfer_dst[0].buffer = asset_buffer;
			buffer_memory_barrier_undefine_to_transfer_dst[0].offset = 0U;
			buffer_memory_barrier_undefine_to_transfer_dst[0].size = VK_WHOLE_SIZE;
			pfn_cmd_pipeline_barrier(vulkan_streaming_transfer_command_buffer, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, VK_PIPELINE_STAGE_TRANSFER_BIT, 0U, 0U, NULL, 1U, buffer_memory_barrier_undefine_to_transfer_dst, 0U, NULL);

			// copy buffer operation
			pfn_cmd_copy_buffer(vulkan_streaming_transfer_command_buffer, staging_buffer, asset_buffer, region_count, regions);

			// queue family ownership transfer - release operation
			VkBufferMemoryBarrier command_buffer_release_ownership[1];
			command_buffer_release_ownership[0].sType = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER;
			command_buffer_release_ownership[0].pNext = NULL;
			command_buffer_release_ownership[0].srcAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
			command_buffer_release_ownership[0].dstAccessMask = 0U;
			command_buffer_release_ownership[0].srcQueueFamilyIndex = vulkan_queue_transfer_family_index;
			command_buffer_release_ownership[0].dstQueueFamilyIndex = vulkan_queue_graphics_family_index;
			command_buffer_release_ownership[0].buffer = asset_buffer;
			command_buffer_release_ownership[0].offset = 0U;
			command_buffer_release_ownership[0].size = VK_WHOLE_SIZE;
			pfn_cmd_pipeline_barrier(vulkan_streaming_transfer_command_buffer, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT, 0U, 0U, NULL, 1U, command_buffer_release_ownership, 0U, NULL);

			// queue family ownership transfer - acquire operation
			VkBufferMemoryBarrier buffer_memory_barrier_acquire_ownership[1];
			buffer_memory_barrier_acquire_ownership[0].sType = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER;
			buffer_memory_barrier_acquire_ownership[0].pNext = NULL;
			buffer_memory_barrier_acquire_ownership[0].srcAccessMask = 0U;
			buffer_memory_barrier_acquire_ownership[0].dstAccessMask = 0U;
			buffer_memory_barrier_acquire_ownership[0].srcQueueFamilyIndex = vulkan_queue_transfer_family_index;
			buffer_memory_barrier_acquire_ownership[0].dstQueueFamilyIndex = vulkan_queue_graphics_family_index;
			buffer_memory_barrier_acquire_ownership[0].buffer = asset_buffer;
			buffer_memory_barrier_acquire_ownership[0].offset = 0U;
			buffer_memory_barrier_acquire_ownership[0].size = VK_WHOLE_SIZE;
			pfn_cmd_pipeline_barrier(vulkan_streaming_graphics_command_buffer, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT, 0U, 0U, NULL, 1U, buffer_memory_barrier_acquire_ownership, 0U, NULL);
		}
		else
		{
			// barrier undefine - transfer_dst
			VkBufferMemoryBarrier buffer_memory_barrier_undefine_to_transfer_dst[1];
			buffer_memory_barrier_undefine_to_transfer_dst[0].sType = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER;
			buffer_memory_barrier_undefine_to_transfer_dst[0].pNext = NULL;
			buffer_memory_barrier_undefine_to_transfer_dst[0].srcAccessMask = 0U;
			buffer_memory_barrier_undefine_to_transfer_dst[0].dstAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
			buffer_memory_barrier_undefine_to_transfer_dst[0].srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
			buffer_memory_barrier_undefine_to_transfer_dst[0].dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
			buffer_memory_barrier_undefine_to_transfer_dst[0].buffer = asset_buffer;
			buffer_memory_barrier_undefine_to_transfer_dst[0].offset = 0U;
			buffer_memory_barrier_undefine_to_transfer_dst[0].size = VK_WHOLE_SIZE;
			pfn_cmd_pipeline_barrier(vulkan_streaming_transfer_command_buffer, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, VK_PIPELINE_STAGE_TRANSFER_BIT, 0U, 0U, NULL, 1U, buffer_memory_barrier_undefine_to_transfer_dst, 0U, NULL);

			// copy buffer operation
			pfn_cmd_copy_buffer(vulkan_streaming_transfer_command_buffer, staging_buffer, asset_buffer, region_count, regions);

			// barrier transfer_dst - asset_read
			VkAccessFlags buffer_memory_barrier_transfer_dst_to_asset_read_dst_access_mask;
			VkPipelineStageFlags buffer_memory_barrier_transfer_dst_to_asset_read_dst_dst_stage_mask;
			switch (asset_buffer_type)
			{
			case ASSET_VERTEX_BUFFER:
				buffer_memory_barrier_transfer_dst_to_asset_read_dst_access_mask = VK_ACCESS_VERTEX_ATTRIBUTE_READ_BIT;
				buffer_memory_barrier_transfer_dst_to_asset_read_dst_dst_stage_mask = VK_PIPELINE_STAGE_VERTEX_INPUT_BIT;
				break;
			case ASSET_INDEX_BUFFER:
				buffer_memory_barrier_transfer_dst_to_asset_read_dst_access_mask = VK_ACCESS_INDEX_READ_BIT;
				buffer_memory_barrier_transfer_dst_to_asset_read_dst_dst_stage_mask = VK_PIPELINE_STAGE_VERTEX_INPUT_BIT;
				break;
			default:
				assert(ASSET_UNIFORM_BUFFER == asset_buffer_type);
				buffer_memory_barrier_transfer_dst_to_asset_read_dst_access_mask = VK_ACCESS_UNIFORM_READ_BIT;
				buffer_memory_barrier_transfer_dst_to_asset_read_dst_dst_stage_mask = VK_PIPELINE_STAGE_VERTEX_SHADER_BIT | VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT;
				break;
			}
			VkBufferMemoryBarrier buffer_memory_barrier_transfer_dst_to_asset_read[1];
			buffer_memory_barrier_transfer_dst_to_asset_read[0].sType = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER;
			buffer_memory_barrier_transfer_dst_to_asset_read[0].pNext = NULL;
			buffer_memory_barrier_transfer_dst_to_asset_read[0].srcAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
			buffer_memory_barrier_transfer_dst_to_asset_read[0].dstAccessMask = buffer_memory_barrier_transfer_dst_to_asset_read_dst_access_mask;
			buffer_memory_barrier_transfer_dst_to_asset_read[0].srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
			buffer_memory_barrier_transfer_dst_to_asset_read[0].dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
			buffer_memory_barrier_transfer_dst_to_asset_read[0].buffer = asset_buffer;
			buffer_memory_barrier_transfer_dst_to_asset_read[0].offset = 0U;
			buffer_memory_barrier_transfer_dst_to_asset_read[0].size = VK_WHOLE_SIZE;
			pfn_cmd_pipeline_barrier(vulkan_streaming_transfer_command_buffer, VK_PIPELINE_STAGE_TRANSFER_BIT, buffer_memory_barrier_transfer_dst_to_asset_read_dst_dst_stage_mask, 0U, 0U, NULL, 1U, buffer_memory_barrier_transfer_dst_to_asset_read, 0U, NULL);
		}
	}
	else
	{
		// barrier undefine - transfer_dst
		VkBufferMemoryBarrier buffer_memory_barrier_undefine_to_transfer_dst[1];
		buffer_memory_barrier_undefine_to_transfer_dst[0].sType = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER;
		buffer_memory_barrier_undefine_to_transfer_dst[0].pNext = NULL;
		buffer_memory_barrier_undefine_to_transfer_dst[0].srcAccessMask = 0U;
		buffer_memory_barrier_undefine_to_transfer_dst[0].dstAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
		buffer_memory_barrier_undefine_to_transfer_dst[0].srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
		buffer_memory_barrier_undefine_to_transfer_dst[0].dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
		buffer_memory_barrier_undefine_to_transfer_dst[0].buffer = asset_buffer;
		buffer_memory_barrier_undefine_to_transfer_dst[0].offset = 0U;
		buffer_memory_barrier_undefine_to_transfer_dst[0].size = VK_WHOLE_SIZE;
		pfn_cmd_pipeline_barrier(vulkan_streaming_graphics_command_buffer, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, VK_PIPELINE_STAGE_TRANSFER_BIT, 0U, 0U, NULL, 1U, buffer_memory_barrier_undefine_to_transfer_dst, 0U, NULL);

		// copy buffer operation
		pfn_cmd_copy_buffer(vulkan_streaming_graphics_command_buffer, staging_buffer, asset_buffer, region_count, regions);

		// barrier transfer_dst - asset_read
		VkAccessFlags buffer_memory_barrier_transfer_dst_to_asset_read_dst_access_mask;
		VkPipelineStageFlags buffer_memory_barrier_transfer_dst_to_asset_read_dst_dst_stage_mask;
		switch (asset_buffer_type)
		{
		case ASSET_VERTEX_BUFFER:
			buffer_memory_barrier_transfer_dst_to_asset_read_dst_access_mask = VK_ACCESS_VERTEX_ATTRIBUTE_READ_BIT;
			buffer_memory_barrier_transfer_dst_to_asset_read_dst_dst_stage_mask = VK_PIPELINE_STAGE_VERTEX_INPUT_BIT;
			break;
		case ASSET_INDEX_BUFFER:
			buffer_memory_barrier_transfer_dst_to_asset_read_dst_access_mask = VK_ACCESS_INDEX_READ_BIT;
			buffer_memory_barrier_transfer_dst_to_asset_read_dst_dst_stage_mask = VK_PIPELINE_STAGE_VERTEX_INPUT_BIT;
			break;
		default:
			assert(ASSET_UNIFORM_BUFFER == asset_buffer_type);
			buffer_memory_barrier_transfer_dst_to_asset_read_dst_access_mask = VK_ACCESS_UNIFORM_READ_BIT;
			buffer_memory_barrier_transfer_dst_to_asset_read_dst_dst_stage_mask = VK_PIPELINE_STAGE_VERTEX_SHADER_BIT | VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT;
			break;
		}
		VkBufferMemoryBarrier buffer_memory_barrier_transfer_dst_to_asset_read[1];
		buffer_memory_barrier_transfer_dst_to_asset_read[0].sType = VK_STRUCTURE_TYPE_BUFFER_MEMORY_BARRIER;
		buffer_memory_barrier_transfer_dst_to_asset_read[0].pNext = NULL;
		buffer_memory_barrier_transfer_dst_to_asset_read[0].srcAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
		buffer_memory_barrier_transfer_dst_to_asset_read[0].dstAccessMask = buffer_memory_barrier_transfer_dst_to_asset_read_dst_access_mask;
		buffer_memory_barrier_transfer_dst_to_asset_read[0].srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
		buffer_memory_barrier_transfer_dst_to_asset_read[0].dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
		buffer_memory_barrier_transfer_dst_to_asset_read[0].buffer = asset_buffer;
		buffer_memory_barrier_transfer_dst_to_asset_read[0].offset = 0U;
		buffer_memory_barrier_transfer_dst_to_asset_read[0].size = VK_WHOLE_SIZE;
		pfn_cmd_pipeline_barrier(vulkan_streaming_graphics_command_buffer, VK_PIPELINE_STAGE_TRANSFER_BIT, buffer_memory_barrier_transfer_dst_to_asset_read_dst_dst_stage_mask, 0U, 0U, NULL, 1U, buffer_memory_barrier_transfer_dst_to_asset_read, 0U, NULL);
	}
}

void streaming_staging_to_asset_image(
	PFN_vkCmdPipelineBarrier pfn_cmd_pipeline_barrier, PFN_vkCmdCopyBufferToImage pfn_cmd_copy_buffer_to_image,
	bool vulkan_has_dedicated_transfer_queue, uint32_t vulkan_queue_transfer_family_index, uint32_t vulkan_queue_graphics_family_index, VkCommandBuffer vulkan_streaming_transfer_command_buffer, VkCommandBuffer vulkan_streaming_graphics_command_buffer,
	VkBuffer staging_buffer, VkImage asset_image, uint32_t region_count, VkBufferImageCopy *const regions, VkImageSubresourceRange const &asset_image_subresource_range)
{
	if (vulkan_has_dedicated_transfer_queue)
	{
		if (vulkan_queue_graphics_family_index != vulkan_queue_transfer_family_index)
		{
			// barrier undefine - transfer_dst
			VkImageMemoryBarrier image_memory_barrier_undefine_to_transfer_dst[1] = {
				{VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER,
				 NULL,
				 0U,
				 VK_ACCESS_TRANSFER_WRITE_BIT,
				 VK_IMAGE_LAYOUT_UNDEFINED,
				 VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
				 VK_QUEUE_FAMILY_IGNORED,
				 VK_QUEUE_FAMILY_IGNORED,
				 asset_image,
				 asset_image_subresource_range}};
			pfn_cmd_pipeline_barrier(vulkan_streaming_transfer_command_buffer, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, VK_PIPELINE_STAGE_TRANSFER_BIT, 0U, 0U, NULL, 0U, NULL, 1U, image_memory_barrier_undefine_to_transfer_dst);

			// copy buffer to image operation
			pfn_cmd_copy_buffer_to_image(vulkan_streaming_transfer_command_buffer, staging_buffer, asset_image, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, region_count, regions);

			// queue family ownership transfer - release operation
			VkImageMemoryBarrier command_buffer_release_ownership[1] = {
				{VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER,
				 NULL,
				 VK_ACCESS_TRANSFER_WRITE_BIT,
				 0,
				 VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
				 VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL,
				 vulkan_queue_transfer_family_index,
				 vulkan_queue_graphics_family_index,
				 asset_image,
				 asset_image_subresource_range}};
			pfn_cmd_pipeline_barrier(vulkan_streaming_transfer_command_buffer, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT, 0U, 0U, NULL, 0U, NULL, 1U, command_buffer_release_ownership);

			// queue family ownership transfer - acquire operation
			VkImageMemoryBarrier image_memory_barrier_acquire_ownership[1] = {
				{VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER,
				 NULL,
				 0,
				 0,
				 VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
				 VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL,
				 vulkan_queue_transfer_family_index,
				 vulkan_queue_graphics_family_index,
				 asset_image,
				 asset_image_subresource_range}};
			pfn_cmd_pipeline_barrier(vulkan_streaming_graphics_command_buffer, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, VK_PIPELINE_STAGE_BOTTOM_OF_PIPE_BIT, 0U, 0U, NULL, 0U, NULL, 1U, image_memory_barrier_acquire_ownership);
		}
		else
		{
			// barrier undefine - transfer_dst
			VkImageMemoryBarrier image_memory_barrier_undefine_to_transfer_dst[1] = {
				{VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER,
				 NULL,
				 0,
				 VK_ACCESS_TRANSFER_WRITE_BIT,
				 VK_IMAGE_LAYOUT_UNDEFINED,
				 VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
				 VK_QUEUE_FAMILY_IGNORED,
				 VK_QUEUE_FAMILY_IGNORED,
				 asset_image,
				 asset_image_subresource_range}};
			pfn_cmd_pipeline_barrier(vulkan_streaming_transfer_command_buffer, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, VK_PIPELINE_STAGE_TRANSFER_BIT, 0U, 0U, NULL, 0U, NULL, 1U, image_memory_barrier_undefine_to_transfer_dst);

			// copy buffer to image operation
			pfn_cmd_copy_buffer_to_image(vulkan_streaming_transfer_command_buffer, staging_buffer, asset_image, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, region_count, regions);

			// barrier transfer_dst - shader_read_only
			VkImageMemoryBarrier image_memory_barrier_transfer_dst_to_shader_read_only[1] = {
				{VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER,
				 NULL,
				 VK_ACCESS_TRANSFER_WRITE_BIT,
				 VK_ACCESS_SHADER_READ_BIT,
				 VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
				 VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL,
				 VK_QUEUE_FAMILY_IGNORED,
				 VK_QUEUE_FAMILY_IGNORED,
				 asset_image,
				 asset_image_subresource_range}};
			pfn_cmd_pipeline_barrier(vulkan_streaming_transfer_command_buffer, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT, 0U, 0U, NULL, 0U, NULL, 1U, image_memory_barrier_transfer_dst_to_shader_read_only);
		}
	}
	else
	{

		// barrier undefine - transfer_dst
		VkImageMemoryBarrier image_memory_barrier_undefine_to_transfer_dst[1] = {
			{VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER,
			 NULL,
			 0,
			 VK_ACCESS_TRANSFER_WRITE_BIT,
			 VK_IMAGE_LAYOUT_UNDEFINED,
			 VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
			 VK_QUEUE_FAMILY_IGNORED,
			 VK_QUEUE_FAMILY_IGNORED,
			 asset_image,
			 asset_image_subresource_range}};

		pfn_cmd_pipeline_barrier(vulkan_streaming_graphics_command_buffer, VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT, VK_PIPELINE_STAGE_TRANSFER_BIT, 0U, 0U, NULL, 0U, NULL, 1U, image_memory_barrier_undefine_to_transfer_dst);

		// copy buffer to image operation
		pfn_cmd_copy_buffer_to_image(vulkan_streaming_graphics_command_buffer, staging_buffer, asset_image, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, region_count, regions);

		// barrier transfer_dst - shader_read_only
		VkImageMemoryBarrier image_memory_barrier_transfer_dst_to_shader_read_only[1] = {
			{VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER,
			 NULL,
			 VK_ACCESS_TRANSFER_WRITE_BIT,
			 VK_ACCESS_SHADER_READ_BIT,
			 VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
			 VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL,
			 VK_QUEUE_FAMILY_IGNORED,
			 VK_QUEUE_FAMILY_IGNORED,
			 asset_image,
			 asset_image_subresource_range}};
		pfn_cmd_pipeline_barrier(vulkan_streaming_graphics_command_buffer, VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT, 0U, 0U, NULL, 0U, NULL, 1U, image_memory_barrier_transfer_dst_to_shader_read_only);
	}
}