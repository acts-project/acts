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

// Darwin x86_64 has full trap support via the legacy x87 control word plus
// SSE MXCSR, so runtime support is always available on this build target.
bool isRuntimeSupported() {
  return true;
}

std::optional<FpeType> decodeFpeType(int signal, const siginfo_t* si,
                                     void* ctx) {
  // On this platform we only install a SIGFPE handler; si_code is enough to
  // classify all FPE kinds we track.
  static_cast<void>(ctx);
  if (signal != SIGFPE || si == nullptr) {
    return std::nullopt;
  }
  return fpeTypeFromSiCode(si->si_code);
}

void clearPendingExceptions(int excepts) {
  // Clear sticky exception state before enabling traps to avoid immediate
  // retriggering from stale flags.
  darwin::clearPendingExceptions(excepts);
}

void enableExceptions(int excepts) {
  // Darwin x86_64 mirrors Linux semantics:
  // - x87 control word: trap enabled when corresponding mask bit is cleared
  // - MXCSR: trap enabled when corresponding mask bit [7:12] is cleared
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
  // Restore masking for requested exceptions in both x87 and SSE domains.
  fenv_t env{};
  if (fegetenv(&env) != 0) {
    return;
  }
  env.__control |= static_cast<unsigned short>(excepts);
  env.__mxcsr |= (static_cast<unsigned int>(excepts) << 7u);
  fesetenv(&env);
}

void maskTrapsInSignalContext(void* ctx, FpeType type) {
  // We mask only the trap that fired in the interrupted context so the faulting
  // instruction can be unwound safely, then clear status flags in both units.
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
  // Use the shared Darwin frame-walker and provide x86_64 register extraction
  // from the interrupted thread state.
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
  // Skip one frame to hide the monitor's own dump helper from final traces.
  return 1;
}

bool shouldFailFastOnUnknownSignal() {
  // x86_64 should only see SIGFPE with known si_code values. Unknowns are
  // treated as non-fatal to preserve backward compatibility.
  return false;
}

void installSignalHandlers(void (*handler)(int, siginfo_t*, void*)) {
  // A single SIGFPE handler is sufficient on Darwin x86_64.
  struct sigaction action {};
  action.sa_sigaction = handler;
  action.sa_flags = SA_SIGINFO;
  sigaction(SIGFPE, &action, nullptr);
}

}  // namespace ActsPlugins::detail
