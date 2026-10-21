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

inline void clearPendingExceptions(int excepts) {
  std::feclearexcept(excepts);
}

template <typename RegisterStateExtractor>
std::size_t captureStackFromSignalContext(void* ctx, void* buffer,
                                          std::size_t bufferBytes,
                                          RegisterStateExtractor&& extractor) {
  // Why this helper exists:
  // In a signal handler we need the stack of the interrupted faulting context
  // (PC/SP/FP from ucontext), not the stack of the handler itself.
  // Generic signal-safe dumps start at the current handler frame and cannot be
  // seeded with those saved registers, so they may miss the real fault site or
  // include mostly signal trampoline frames. Walking from saved FP/PC gives a
  // deterministic trace rooted at the trapping instruction on Darwin.
  using NativeFramePtr = boost::stacktrace::frame::native_frame_ptr_t;
  auto* frames = static_cast<NativeFramePtr*>(buffer);
  const std::size_t maxFrames = bufferBytes / sizeof(NativeFramePtr);
  std::size_t count = 0;

  if (ctx == nullptr || maxFrames == 0) {
    return 0;
  }

  // The platform TU provides arch-specific extraction of SP/FP/PC from the raw
  // signal context while this helper keeps the frame-walk logic shared.
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

  // Keep unwinding constrained to a finite stack window above SP to avoid
  // dereferencing arbitrary memory if the frame chain is corrupted.
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

  // Standard frame-pointer chain walk:
  //   fp -> {prev_fp, return_address}
  // Stops when the chain is non-monotonic, leaves the allowed stack window,
  // or we exhaust caller-provided buffer capacity.
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
