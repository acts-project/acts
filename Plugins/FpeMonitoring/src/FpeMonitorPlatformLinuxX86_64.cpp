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

// Linux x86_64 signal context exposes raw x87/SSE control and status words.
// These masks are used when disabling the current trap and clearing sticky
// exception bits in the interrupted context.
constexpr std::uint16_t kFpuExceptionMask = 0x3f;
constexpr std::uint16_t kFpuStatusFlags = 0xff;
constexpr std::uint32_t kSseStatusFlags = kFpuExceptionMask;

}  // namespace

bool isRuntimeSupported() {
  // Linux x86_64 supports feenableexcept/fedisableexcept and SIGFPE si_code
  // decoding, so trapping mode is fully available.
  return true;
}

std::optional<FpeType> decodeFpeType(int signal, siginfo_t* si, void* ctx) {
  // This backend only installs SIGFPE handlers. si_code carries the exception
  // category we expose in FpeType.
  static_cast<void>(ctx);
  if (signal != SIGFPE || si == nullptr) {
    return std::nullopt;
  }
  return fpeTypeFromSiCode(si->si_code);
}

void clearPendingExceptions(int excepts) {
  // Clear stale sticky flags before changing trap state.
  std::feclearexcept(excepts);
}

void enableExceptions(int excepts) {
  // glibc helper enables hardware traps for requested FE_* bits.
  feenableexcept(excepts);
}

void disableExceptions(int excepts) {
  // glibc helper disables hardware traps for requested FE_* bits.
  fedisableexcept(excepts);
}

void maskTrapsInSignalContext(void* ctx, FpeType type) {
  // Linux reports enough detail in si_code, so "type" is currently unused.
  // We mask all x87 trap bits and clear pending x87/SSE status flags to allow
  // safe unwinding past the faulting instruction.
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
  // Linux x86_64 currently falls back to boost::stacktrace::safe_dump_to in
  // the signal handler, so no context-based frame extraction is attempted here.
  static_cast<void>(ctx);
  static_cast<void>(buffer);
  static_cast<void>(bufferBytes);
  return 0;
}

std::size_t safeDumpSkipFrames() {
  // Skip two frames to hide signal-handler and dump-helper internals.
  return 2;
}

bool shouldFailFastOnUnknownSignal() {
  // Unknown/unsupported signal payloads are tolerated and ignored on Linux.
  return false;
}

void installSignalHandlers(void (*handler)(int, siginfo_t*, void*)) {
  // Linux only needs SIGFPE for floating-point trap monitoring.
  struct sigaction action {};
  action.sa_sigaction = handler;
  action.sa_flags = SA_SIGINFO;
  sigaction(SIGFPE, &action, nullptr);
}

}  // namespace ActsPlugins::detail
