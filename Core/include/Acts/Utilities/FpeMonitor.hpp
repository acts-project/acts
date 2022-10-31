// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Helpers.hpp"

// #include <backward.hpp>
#include <cfenv>
#include <iostream>

#include <fenv.h>
#include <signal.h>

namespace Acts {

static void ill_signal_handler(int /*sig*/, siginfo_t* /*sip*/, void* /*scp*/) {
  std::cout << "Caught SIGILL (likely via FPE signal handler)" << std::endl;

  // using namespace backward;
  // StackTrace st;
  // st.load_here(32);
  // Printer p;
  // p.print(st, stderr);
  // std::abort();
  std::feclearexcept(FE_ALL_EXCEPT);
}

static void fpe_signal_handler(int /*sig*/, siginfo_t* sip, void* /*scp*/) {
  int fe_code = sip->si_code;

  switch (fe_code) {
    case FPE_INTDIV:
      std::cout << "Integer divide by zero" << std::endl;
      break;
    case FPE_INTOVF:
      std::cout << "Integer overflow" << std::endl;
      break;
    case FPE_FLTDIV:
      std::cout << "Floating point divide by zero" << std::endl;
      break;
    case FPE_FLTOVF:
      std::cout << "Floating point overflow" << std::endl;
      break;
    case FPE_FLTUND:
      std::cout << "Floating point underflow" << std::endl;
      break;
    case FPE_FLTRES:
      std::cout << "Floating point inexact result" << std::endl;
      break;
    case FPE_FLTINV:
      std::cout << "Floating point invalid operation" << std::endl;
      break;
    case FPE_FLTSUB:
      std::cout << "Floating point subscript out of range" << std::endl;
      break;
    default:
      std::cerr << "Unknown signal caught:" << fe_code << std::endl;
      std::abort();
  }

  // using namespace backward;
  // StackTrace st;
  // st.load_here(32);
  // Printer p;
  // p.print(st, stderr);

  std::abort();
}

#if defined(__APPLE__)
inline int feenableexcept(unsigned int excepts) {
#if defined(__aarch64__)
  fenv_t env;
  fegetenv(&env);

  if (ACTS_CHECK_BIT(excepts, FE_DIVBYZERO)) {
    env.__fpcr |= __fpcr_trap_divbyzero;
  }

  if (ACTS_CHECK_BIT(excepts, FE_INVALID)) {
    env.__fpcr |= __fpcr_trap_invalid;
  }

  if (ACTS_CHECK_BIT(excepts, FE_OVERFLOW)) {
    env.__fpcr |= __fpcr_trap_overflow;
  }

  return fesetenv(&env);
#else
  ACTS_LOCAL_LOGGER(
      Acts::getDefaultLogger("feenableexcept", Acts::Logging::WARNING));
  ACTS_WARNING(
      "FPE exception raising currently not implemented for macOS x86_64");
#endif

  // unmask
  // fenv.__control &= ~new_excepts;
  // fenv.__mxcsr &= ~(new_excepts << 7);
}
#endif

void enable_floating_point_exceptions() {
  // macOS
  // fegetenv(&env);
  // env.__fpcr = env.__fpcr | __fpcr_trap_invalid | __fpcr_trap_divbyzero |
  // __fpcr_trap_overflow;
  // fesetenv(&env);

  // std::fenv_t env;
  // std::fegetenv(&env);
  // std::fesetenv(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  {
    struct sigaction act {};
    act.sa_sigaction = fpe_signal_handler;
    sigemptyset(&act.sa_mask);
    act.sa_flags = SA_SIGINFO;
    sigaction(SIGFPE, &act, nullptr);
  }

  {
    struct sigaction act {};
    act.sa_sigaction = ill_signal_handler;
    sigemptyset(&act.sa_mask);
    act.sa_flags = SA_SIGINFO;
    sigaction(SIGILL, &act, nullptr);
  }
}
}  // namespace Acts
