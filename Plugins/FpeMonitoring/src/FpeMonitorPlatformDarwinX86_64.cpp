// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/FpeMonitoring/FpeMonitor.hpp"

#include <cfenv>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <optional>

#include <boost/stacktrace/frame.hpp>

#include "FpeMonitorPlatform.hpp"

#if !defined(__APPLE__) || !defined(__x86_64__)
#error "This translation unit is only valid for macOS x86_64"
#endif

namespace ActsPlugins::detail {
namespace {

int exceptMaskForType(FpeType type) {
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

std::optional<FpeType> fpeTypeFromSiCode(int siCode) {
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

}  // namespace

bool isRuntimeSupported() {
  return true;
}

std::optional<FpeType> decodeFpeType(int signal, siginfo_t* si, void* ctx) {
  static_cast<void>(ctx);
  if (signal != SIGFPE || si == nullptr) {
    return std::nullopt;
  }
  return fpeTypeFromSiCode(si->si_code);
}

void clearPendingExceptions(int excepts) {
  std::feclearexcept(excepts);
}

void enableExceptions(int excepts) {
  fenv_t env{};
  if (fegetenv(&env) != 0) {
    return;
  }
  env.__control &= ~static_cast<unsigned short>(excepts);
  env.__mxcsr &= ~(static_cast<unsigned int>(excepts) << 7u);
  env.__status &= ~static_cast<unsigned short>(FE_ALL_EXCEPT);
  env.__mxcsr &= ~static_cast<unsigned int>(FE_ALL_EXCEPT);
  fesetenv(&env);
}

void disableExceptions(int excepts) {
  fenv_t env{};
  if (fegetenv(&env) != 0) {
    return;
  }
  env.__control |= static_cast<unsigned short>(excepts);
  env.__mxcsr |= (static_cast<unsigned int>(excepts) << 7u);
  fesetenv(&env);
}

void maskTrapsInSignalContext(void* ctx, FpeType type) {
  const int excepts = exceptMaskForType(type);
  auto* uc = static_cast<ucontext_t*>(ctx);
  uc->uc_mcontext->__fs.__fpu_fcw |= static_cast<unsigned short>(excepts);
  uc->uc_mcontext->__fs.__fpu_fsw &=
      ~static_cast<unsigned short>(FE_ALL_EXCEPT);
  uc->uc_mcontext->__fs.__fpu_mxcsr |=
      (static_cast<unsigned int>(excepts) << 7u);
  uc->uc_mcontext->__fs.__fpu_mxcsr &=
      ~static_cast<unsigned int>(FE_ALL_EXCEPT);
}

std::size_t captureStackFromSignalContext(void* ctx, void* buffer,
                                          std::size_t bufferBytes) {
  using NativeFramePtr = boost::stacktrace::frame::native_frame_ptr_t;
  auto* frames = static_cast<NativeFramePtr*>(buffer);
  const std::size_t maxFrames = bufferBytes / sizeof(NativeFramePtr);
  std::size_t count = 0;

  if (ctx == nullptr || maxFrames == 0) {
    return 0;
  }

  auto* uc = static_cast<ucontext_t*>(ctx);
  const std::uintptr_t sp = uc->uc_mcontext->__ss.__rsp;
  std::uintptr_t fp = uc->uc_mcontext->__ss.__rbp;
  const std::uintptr_t pc = uc->uc_mcontext->__ss.__rip;

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
    const std::uintptr_t ra = record->returnAddress;
    push(ra);

    if (prevFp <= fp || !inStackWindow(prevFp)) {
      break;
    }
    fp = prevFp;
  }

  return count * sizeof(NativeFramePtr);
}

std::size_t safeDumpSkipFrames() {
  return 1;
}

bool shouldFailFastOnUnknownSignal() {
  return false;
}

void installSignalHandlers(void (*handler)(int, siginfo_t*, void*)) {
  struct sigaction action {};
  action.sa_sigaction = handler;
  action.sa_flags = SA_SIGINFO;
  sigaction(SIGFPE, &action, nullptr);
}

}  // namespace ActsPlugins::detail
