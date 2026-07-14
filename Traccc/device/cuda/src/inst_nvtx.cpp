/* Copyright (c) 2014, NVIDIA CORPORATION. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#if __has_include(<cxxabi.h>)
#define HAVE_CXXABI 1
#include <cxxabi.h>
#endif

#include <dlfcn.h>

#include <cstring>
#include <memory>
#include <optional>
#include <string>

#include "nvtx3/nvToolsExt.h"

namespace {
uint32_t __attribute__((no_instrument_function)) djb2(const std::string &str) {
    uint32_t hash = 5381u;

    for (auto &chr : str) {
        hash = 33u * hash + static_cast<uint32_t>(chr);
    }

    return hash;
}

void __attribute__((no_instrument_function)) nvtxRangePushWrapper(
    const std::optional<std::string> &name) {
    nvtxEventAttributes_t eventAttrib;
    std::memset(&eventAttrib, 0, sizeof(nvtxEventAttributes_t));

    eventAttrib.version = NVTX_VERSION;
    eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;

    if (name) {
        eventAttrib.colorType = NVTX_COLOR_ARGB;
        eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
        eventAttrib.color = 0xFF000000 | (djb2(*name) & 0x00FFFFFF);
        eventAttrib.message.ascii = name->c_str();
    } else {
        eventAttrib.colorType = NVTX_COLOR_UNKNOWN;
        eventAttrib.messageType = NVTX_MESSAGE_UNKNOWN;
    }

    nvtxRangePushEx(&eventAttrib);
}
}  // namespace

extern "C" {
void __attribute__((no_instrument_function)) __cyg_profile_func_enter(
    [[maybe_unused]] void *this_fn, void *) {
    Dl_info this_fn_info;

    if (dladdr(this_fn, &this_fn_info)) {
#ifdef HAVE_CXXABI
        int status = 0;
        std::unique_ptr<char[]> fname(abi::__cxa_demangle(
            this_fn_info.dli_sname, nullptr, nullptr, &status));

        if (status == 0 && fname) {
            nvtxRangePushWrapper(fname.get());
        } else {
            nvtxRangePushWrapper({});
        }
#else
        nvtxRangePushWrapper(this_fn_info.dli_sname);
#endif

    } else {
        nvtxRangePushWrapper({});
    }
}

void __attribute__((no_instrument_function)) __cyg_profile_func_exit(void *,
                                                                     void *) {
    nvtxRangePop();
}
}
