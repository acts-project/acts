// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cfenv>
#include <cstddef>
#include <cstdint>
#include <optional>

#include <boost/stacktrace/frame.hpp>

#include "FpeMonitorPlatform.hpp"

namespace ActsPlugins::detail::darwin {

struct RegisterState {
  std::uintptr_t sp;
  std::uintptr_t fp;
  std::uintptr_t pc;
};

inline int exceptMaskForType(FpeType type) {
  using enum FpeType;
  switch (type) {
    case INTDIV:
    case FLTDIV:
      return FE_DIVBYZERO;
    case INTOVF:
    case FLTOVF:
      return FE_OVERFLOW;
    case FLTUND:
      return FE_UNDERFLOW;
    case FLTRES:
      return FE_INEXACT;
    case FLTINV:
    case FLTSUB:
      return FE_INVALID;
    default:
      return 0;
  }
}

inline std::optional<FpeType> fpeTypeFromSiCode(int siCode) {
  using enum FpeType;
  switch (siCode) {
    case FPE_INTDIV:
      return INTDIV;
    case FPE_INTOVF:
      return INTOVF;
    case FPE_FLTDIV:
      return FLTDIV;
    case FPE_FLTOVF:
      return FLTOVF;
    case FPE_FLTUND:
      return FLTUND;
    case FPE_FLTRES:
      return FLTRES;
    case FPE_FLTINV:
      return FLTINV;
    case FPE_FLTSUB:
      return FLTSUB;
    default:
      return std::nullopt;
  }
}

inline void clearPendingExceptions(int excepts) {
  std::feclearexcept(excepts);
}

template <typename RegisterStateExtractor>
std::size_t captureStackFromSignalContext(void* ctx, void* buffer,
                                          std::size_t bufferBytes,
                                          RegisterStateExtractor&& extractor) {
  using NativeFramePtr = boost::stacktrace::frame::native_frame_ptr_t;
  auto* frames = static_cast<NativeFramePtr*>(buffer);
  const std::size_t maxFrames = bufferBytes / sizeof(NativeFramePtr);
  std::size_t count = 0;

  if (ctx == nullptr || maxFrames == 0) {
    return 0;
  }

  const RegisterState state = extractor(ctx);
  const std::uintptr_t sp = state.sp;
  std::uintptr_t pc = state.pc;
  std::uintptr_t fp = state.fp;

  auto push = [&](std::uintptr_t address) {
    if (address == 0 || count >= maxFrames) {
      return;
    }
    frames[count++] = reinterpret_cast<NativeFramePtr>(address);
  };

  push(pc);

  constexpr std::uintptr_t kMaxStackWindow = 16 * 1024 * 1024;
  auto inStackWindow = [&](std::uintptr_t address) {
    if (address < sp || address > sp + kMaxStackWindow) {
      return false;
    }
    return (address % alignof(std::uintptr_t)) == 0;
  };

  struct FrameRecord {
    std::uintptr_t prevFp;
    std::uintptr_t returnAddress;
  };

  while (count < maxFrames && inStackWindow(fp) &&
         fp + sizeof(FrameRecord) <= sp + kMaxStackWindow) {
    const auto* record = reinterpret_cast<const FrameRecord*>(fp);
    const std::uintptr_t prevFp = record->prevFp;
    const std::uintptr_t returnAddress = record->returnAddress;
    push(returnAddress);

    if (prevFp <= fp || !inStackWindow(prevFp)) {
      break;
    }
    fp = prevFp;
  }

  return count * sizeof(NativeFramePtr);
}

}  // namespace ActsPlugins::detail::darwin
