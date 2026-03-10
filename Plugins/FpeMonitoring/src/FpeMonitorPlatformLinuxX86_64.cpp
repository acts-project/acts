// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <cfenv>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <optional>

#include <ucontext.h>

#include "FpeMonitorPlatform.hpp"

#if !defined(__linux__) || !defined(__x86_64__)
#error "This translation unit is only valid for Linux x86_64"
#endif

namespace ActsPlugins::detail {
namespace {

constexpr std::uint16_t kFpuExceptionMask = 0x3f;
constexpr std::uint16_t kFpuStatusFlags = 0xff;
constexpr std::uint32_t kSseStatusFlags = kFpuExceptionMask;

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
  feenableexcept(excepts);
}

void disableExceptions(int excepts) {
  fedisableexcept(excepts);
}

void maskTrapsInSignalContext(void* ctx, FpeType type) {
  static_cast<void>(type);
  auto* uc = static_cast<ucontext_t*>(ctx);
  __uint16_t* cw = &uc->uc_mcontext.fpregs->cwd;
  *cw |= kFpuExceptionMask;

  __uint16_t* sw = &uc->uc_mcontext.fpregs->swd;
  *sw &= ~kFpuStatusFlags;

  __uint32_t* mxcsr = &uc->uc_mcontext.fpregs->mxcsr;
  *mxcsr |= ((*mxcsr & kSseStatusFlags) << 7);
  *mxcsr &= ~kSseStatusFlags;
}

std::size_t captureStackFromSignalContext(void* ctx, void* buffer,
                                          std::size_t bufferBytes) {
  static_cast<void>(ctx);
  static_cast<void>(buffer);
  static_cast<void>(bufferBytes);
  return 0;
}

std::size_t safeDumpSkipFrames() {
  return 2;
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
