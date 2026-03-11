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
#include <optional>

#include "FpeMonitorPlatform.hpp"
#include "FpeMonitorPlatformDarwinCommon.hpp"

#if !defined(__APPLE__) || !defined(__x86_64__)
#error "This translation unit is only valid for macOS x86_64"
#endif

namespace ActsPlugins::detail {

bool isRuntimeSupported() {
  return true;
}

std::optional<FpeType> decodeFpeType(int signal, siginfo_t* si, void* ctx) {
  static_cast<void>(ctx);
  if (signal != SIGFPE || si == nullptr) {
    return std::nullopt;
  }
  return darwin::fpeTypeFromSiCode(si->si_code);
}

void clearPendingExceptions(int excepts) {
  darwin::clearPendingExceptions(excepts);
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
  const int excepts = darwin::exceptMaskForType(type);
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
  return darwin::captureStackFromSignalContext(
      ctx, buffer, bufferBytes, [](void* rawCtx) {
        const auto& uc = *static_cast<ucontext_t*>(rawCtx);
        return darwin::RegisterState{
            .sp = uc.uc_mcontext->__ss.__rsp,
            .fp = uc.uc_mcontext->__ss.__rbp,
            .pc = uc.uc_mcontext->__ss.__rip,
        };
      });
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
