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

#include "FpeMonitorPlatform.hpp"
#include "FpeMonitorPlatformDarwinCommon.hpp"

#if !defined(__APPLE__) || !defined(__arm64__)
#error "This translation unit is only valid for macOS arm64"
#endif

namespace ActsPlugins::detail {
namespace {

std::uint32_t darwinArm64TrapMask(int excepts) {
  std::uint32_t mask = 0;
  if ((excepts & FE_INVALID) != 0) {
    mask |= __fpcr_trap_invalid;
  }
  if ((excepts & FE_DIVBYZERO) != 0) {
    mask |= __fpcr_trap_divbyzero;
  }
  if ((excepts & FE_OVERFLOW) != 0) {
    mask |= __fpcr_trap_overflow;
  }
  if ((excepts & FE_UNDERFLOW) != 0) {
    mask |= __fpcr_trap_underflow;
  }
  if ((excepts & FE_INEXACT) != 0) {
    mask |= __fpcr_trap_inexact;
  }
  return mask;
}

std::optional<FpeType> fpeTypeFromDarwinArm64Esr(std::uint32_t esr) {
  constexpr std::uint32_t kEsrExceptionClassShift = 26u;
  constexpr std::uint32_t kEsrExceptionClassMask = 0x3fu;
  constexpr std::uint32_t kFpExceptionClass = 0x2cu;
  const std::uint32_t exceptionClass =
      (esr >> kEsrExceptionClassShift) & kEsrExceptionClassMask;
  if (exceptionClass != kFpExceptionClass) {
    return std::nullopt;
  }

  const std::uint32_t flags = esr & static_cast<std::uint32_t>(FE_ALL_EXCEPT);
  if ((flags & FE_INVALID) != 0) {
    return FpeType::FLTINV;
  }
  if ((flags & FE_DIVBYZERO) != 0) {
    return FpeType::FLTDIV;
  }
  if ((flags & FE_OVERFLOW) != 0) {
    return FpeType::FLTOVF;
  }
  if ((flags & FE_UNDERFLOW) != 0) {
    return FpeType::FLTUND;
  }
  if ((flags & FE_INEXACT) != 0) {
    return FpeType::FLTRES;
  }
  return std::nullopt;
}

}  // namespace

bool isRuntimeSupported() {
  return true;
}

std::optional<FpeType> decodeFpeType(int signal, siginfo_t* si, void* ctx) {
  if (signal == SIGFPE && si != nullptr) {
    return darwin::fpeTypeFromSiCode(si->si_code);
  }

  if (signal == SIGILL && ctx != nullptr) {
    auto* uc = static_cast<ucontext_t*>(ctx);
    return fpeTypeFromDarwinArm64Esr(uc->uc_mcontext->__es.__esr);
  }

  return std::nullopt;
}

void clearPendingExceptions(int excepts) {
  darwin::clearPendingExceptions(excepts);
}

void enableExceptions(int excepts) {
  fenv_t env{};
  if (fegetenv(&env) != 0) {
    return;
  }
  env.__fpcr |= static_cast<unsigned long long>(darwinArm64TrapMask(excepts));
  env.__fpsr &= ~static_cast<unsigned long long>(FE_ALL_EXCEPT);
  fesetenv(&env);
}

void disableExceptions(int excepts) {
  fenv_t env{};
  if (fegetenv(&env) != 0) {
    return;
  }
  env.__fpcr &= ~static_cast<unsigned long long>(darwinArm64TrapMask(excepts));
  fesetenv(&env);
}

void maskTrapsInSignalContext(void* ctx, FpeType type) {
  const int excepts = darwin::exceptMaskForType(type);
  auto* uc = static_cast<ucontext_t*>(ctx);
  uc->uc_mcontext->__ns.__fpcr &=
      ~static_cast<std::uint32_t>(darwinArm64TrapMask(excepts));
  uc->uc_mcontext->__ns.__fpsr &= ~static_cast<std::uint32_t>(FE_ALL_EXCEPT);
}

std::size_t captureStackFromSignalContext(void* ctx, void* buffer,
                                          std::size_t bufferBytes) {
  return darwin::captureStackFromSignalContext(
      ctx, buffer, bufferBytes, [](void* rawCtx) {
        const auto& uc = *static_cast<ucontext_t*>(rawCtx);
        return darwin::RegisterState{
            .sp = __darwin_arm_thread_state64_get_sp(uc.uc_mcontext->__ss),
            .fp = __darwin_arm_thread_state64_get_fp(uc.uc_mcontext->__ss),
            .pc = __darwin_arm_thread_state64_get_pc(uc.uc_mcontext->__ss),
        };
      });
}

std::size_t safeDumpSkipFrames() {
  return 1;
}

bool shouldFailFastOnUnknownSignal() {
  return true;
}

void installSignalHandlers(void (*handler)(int, siginfo_t*, void*)) {
  struct sigaction action {};
  action.sa_sigaction = handler;
  action.sa_flags = SA_SIGINFO;
  sigaction(SIGFPE, &action, nullptr);
  sigaction(SIGILL, &action, nullptr);
}

}  // namespace ActsPlugins::detail
