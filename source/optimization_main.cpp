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

#include "../thirdparty/Linearly-Transformed-Cosine/include/brx_linearly_transformed_cosine_optimization.h"

#include <DirectXPackedVector.h>
#include <cinttypes>

static inline uintptr_t internal_brx_open(const char *path);

static inline void internal_brx_write(uintptr_t fildes, void const *buf, size_t nbyte);

static inline void internal_brx_close(uintptr_t fildes);

int main(int argc, char *argv[])
{
    constexpr uint32_t const ltc_lut_width = 128U;
    constexpr uint32_t const ltc_lut_height = 128U;

    mcrt_vector<DirectX::XMFLOAT4> ltc_lut_matrices(static_cast<size_t>(ltc_lut_width * ltc_lut_height));
    brx_linearly_transformed_cosine_optimization_compute_matrices(ltc_lut_matrices.data(), ltc_lut_width, ltc_lut_height);

    // write DDS
    {
        uintptr_t file = internal_brx_open("brx_linearly_transformed_cosine_look_up_table_matrices.dds");

        {
            uint32_t dds_metadata[] =
                {
                    // DDS_MAGIC
                    0X20534444U,
                    // sizeof(DDS_HEADER)
                    124U,
                    // DDSD_CAPS | DDSD_HEIGHT | DDSD_WIDTH | DDSD_PIXELFORMAT
                    0X1007U,
                    // Height
                    ltc_lut_height,
                    // Width
                    ltc_lut_width,
                    // PitchOrLinearSize
                    0U,
                    // Depth,
                    0U,
                    // MipMapCount
                    0U,
                    // Reserved1[11]
                    0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U, 0U,
                    // sizeof(DDS_PIXELFORMAT)
                    32U,
                    // DDPF_FOURCC
                    0X4U,
                    // D3DFMT_A16B16G16R16F
                    113U,
                    // RGBBitCount
                    0U,
                    // RBitMask
                    0U,
                    // GBitMask
                    0U,
                    // BBitMask
                    0U,
                    // ABitMask
                    0U,
                    // DDSCAPS_TEXTURE
                    0X1000,
                    // Caps2
                    0U,
                    // Caps3
                    0U,
                    // Caps4
                    0U,
                    // Reserved2
                    0U};

            internal_brx_write(file, dds_metadata, sizeof(dds_metadata));
        }

        {
            mcrt_vector<uint16_t> half(static_cast<size_t>(4U * ltc_lut_width * ltc_lut_height));

            for (uint32_t ltc_lut_index = 0U; ltc_lut_index < (ltc_lut_width * ltc_lut_height); ++ltc_lut_index)
            {
                half[4U * ltc_lut_index] = DirectX::PackedVector::XMConvertFloatToHalf(ltc_lut_matrices[ltc_lut_index].x);
                half[4U * ltc_lut_index + 1U] = DirectX::PackedVector::XMConvertFloatToHalf(ltc_lut_matrices[ltc_lut_index].y);
                half[4U * ltc_lut_index + 2U] = DirectX::PackedVector::XMConvertFloatToHalf(ltc_lut_matrices[ltc_lut_index].z);
                half[4U * ltc_lut_index + 3U] = DirectX::PackedVector::XMConvertFloatToHalf(ltc_lut_matrices[ltc_lut_index].w);
            }

            assert((sizeof(uint16_t) * half.size()) == (sizeof(uint16_t) * 4U * ltc_lut_width * ltc_lut_height));

            internal_brx_write(file, half.data(), sizeof(uint16_t) * half.size());
        }

        internal_brx_close(file);
    }

    // write C header
    {
        uintptr_t file = internal_brx_open("../thirdparty/Linearly-Transformed-Cosine/include/brx_linearly_transformed_cosine_look_up_table_matrices.h");

        {
            constexpr char const string[] = {"//\r\n// Copyright (C) YuqiaoZhang(HanetakaChou)\r\n//\r\n// This program is free software: you can redistribute it and/or modify\r\n// it under the terms of the GNU Lesser General Public License as published\r\n// by the Free Software Foundation, either version 3 of the License, or\r\n// (at your option) any later version.\r\n//\r\n// This program is distributed in the hope that it will be useful,\r\n// but WITHOUT ANY WARRANTY; without even the implied warranty of\r\n// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\r\n// GNU Lesser General Public License for more details.\r\n//\r\n// You should have received a copy of the GNU Lesser General Public License\r\n// along with this program.  If not, see <https://www.gnu.org/licenses/>.\r\n//\r\n\r\n#ifndef _BRX_LINEARLY_TRANSFORMED_COSINE_LOOK_UP_TABLE_MATRICES_H_\r\n#define _BRX_LINEARLY_TRANSFORMED_COSINE_LOOK_UP_TABLE_MATRICES_H_ 1\r\n\r\n"};

            internal_brx_write(file, string, (sizeof(string) / sizeof(string[0]) - 1U));
        }

        {
            constexpr char const format[] = {"static constexpr uint32_t const g_brx_ltc_lut_width = %dU;\r\n"};

            char string[256];
            int const nchar_written = std::snprintf(string, (sizeof(string) / sizeof(string[0])), format, static_cast<int>(ltc_lut_width));
            assert(nchar_written < (sizeof(string) / sizeof(string[0])));

            internal_brx_write(file, string, nchar_written);
        }

        {
            constexpr char const format[] = {"static constexpr uint32_t const g_brx_ltc_lut_height = %dU;\r\n"};

            char string[256];
            int const nchar_written = std::snprintf(string, (sizeof(string) / sizeof(string[0])), format, static_cast<int>(ltc_lut_height));
            assert(nchar_written < (sizeof(string) / sizeof(string[0])));

            internal_brx_write(file, string, nchar_written);
        }

        {
            constexpr char const string[] = {"static constexpr uint16_t const g_brx_ltc_lut_matrices[4U * g_brx_ltc_lut_width * g_brx_ltc_lut_height] = {\r\n"};

            internal_brx_write(file, string, (sizeof(string) / sizeof(string[0]) - 1U));
        }

        for (uint32_t ltc_lut_index = 0U; ltc_lut_index < (ltc_lut_width * ltc_lut_height); ++ltc_lut_index)
        {
            uint16_t const half2_x = DirectX::PackedVector::XMConvertFloatToHalf(ltc_lut_matrices[ltc_lut_index].x);
            uint16_t const half2_y = DirectX::PackedVector::XMConvertFloatToHalf(ltc_lut_matrices[ltc_lut_index].y);
            uint16_t const half2_z = DirectX::PackedVector::XMConvertFloatToHalf(ltc_lut_matrices[ltc_lut_index].z);
            uint16_t const half2_w = DirectX::PackedVector::XMConvertFloatToHalf(ltc_lut_matrices[ltc_lut_index].w);

            constexpr char const format[] = {"    0X%04" PRIX16 "U, 0X%04" PRIX16 "U, 0X%04" PRIX16 "U, 0X%04" PRIX16 "U,\r\n"};

            char string[256];
            int const nchar_written = std::snprintf(string, (sizeof(string) / sizeof(string[0])), format, half2_x, half2_y, half2_z, half2_w);
            assert(nchar_written < (sizeof(string) / sizeof(string[0])));

            internal_brx_write(file, string, nchar_written);
        }

        {
            constexpr char const string[] = {"};\r\n\r\n#endif\r\n"};
            internal_brx_write(file, string, (sizeof(string) / sizeof(string[0]) - 1U));
        }

        internal_brx_close(file);
    }

    return 0;
}

#define NOMINMAX 1
#define WIN32_LEAN_AND_MEAN 1
#include <sdkddkver.h>
#include <Windows.h>

static inline uintptr_t internal_brx_open(const char *path)
{
    HANDLE const file = CreateFileA(path, FILE_GENERIC_WRITE, 0U, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    assert(INVALID_HANDLE_VALUE != file);
    return reinterpret_cast<uintptr_t>(file);
}

static inline void internal_brx_write(uintptr_t fildes, void const *buf, size_t nbyte)
{
    HANDLE const file = reinterpret_cast<HANDLE>(fildes);
    assert(INVALID_HANDLE_VALUE != file);

    DWORD nbyte_written;
    BOOL result_write_file = WriteFile(file, buf, static_cast<DWORD>(nbyte), &nbyte_written, NULL);
    assert(FALSE != result_write_file);
    assert(nbyte == nbyte_written);
}

static inline void internal_brx_close(uintptr_t fildes)
{
    HANDLE const file = reinterpret_cast<HANDLE>(fildes);
    assert(INVALID_HANDLE_VALUE != file);

    BOOL result_close_handle = CloseHandle(file);
    assert(FALSE != result_close_handle);
}
