// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/FpeMonitor.hpp"

#include "Acts/Utilities/Helpers.hpp"

#include <bitset>
#include <cfenv>
#include <csignal>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string_view>

#include <boost/stacktrace.hpp>
#include <fenv.h>
#include <signal.h>

#define FPU_EXCEPTION_MASK 0x3f
#define FPU_STATUS_FLAGS 0xff
#define SSE_STATUS_FLAGS FPU_EXCEPTION_MASK
#define SSE_EXCEPTION_MASK (FPU_EXCEPTION_MASK << 7)

namespace Acts {

FpeMonitor::FpeMonitor()
    : FpeMonitor{FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW} {}

FpeMonitor::FpeMonitor(int excepts) : m_excepts(excepts) {
  enable();
}

FpeMonitor::~FpeMonitor() {
  disable();
}

void FpeMonitor::signalHandler(int /*signal*/, siginfo_t *si, void *ctx) {
  std::cout << "SIGACTION" << std::endl;
  switch (si->si_code) {
    case FPE_INTDIV:
      std::cerr << "integer divide by zero";
      break;
    case FPE_INTOVF:
      std::cerr << "integer overflow";
      break;
    case FPE_FLTDIV:
      std::cerr << "floating point divide by zero";
      break;
    case FPE_FLTOVF:
      std::cerr << "floating point overflow";
      break;
    case FPE_FLTUND:
      std::cerr << "floating point underflow";
      break;
    case FPE_FLTRES:
      std::cerr << "floating point inexact result";
      break;
    case FPE_FLTINV:
      std::cerr << "floating point invalid operation";
      break;
    case FPE_FLTSUB:
      std::cerr << "subscript out of range";
      break;
  }

  // std::cout << boost::stacktrace::stacktrace();
  std::cout << std::endl;

  FpeMonitor &fpe = *stack().top();
  fpe.m_encountered |= 1 << si->si_code;

  __uint16_t *cw = &((ucontext_t *)ctx)->uc_mcontext.fpregs->cwd;
  *cw |= FPU_EXCEPTION_MASK;

  __uint16_t *sw = &((ucontext_t *)ctx)->uc_mcontext.fpregs->swd;
  *sw &= ~FPU_STATUS_FLAGS;

  __uint32_t *mxcsr = &((ucontext_t *)ctx)->uc_mcontext.fpregs->mxcsr;
  // *mxcsr |= SSE_EXCEPTION_MASK; [> disable all SSE exceptions <]
  *mxcsr |= ((*mxcsr & SSE_STATUS_FLAGS) << 7);
  *mxcsr &= ~SSE_STATUS_FLAGS; /* clear all pending SSE exceptions */
}

void FpeMonitor::enable() {
#if defined(__APPLE__)
  std::cerr << "FPE monitoring currently not supported on Apple" << std::endl;
#else
  std::feclearexcept(m_excepts);
  feenableexcept(m_excepts);

  ensureSignalHandlerInstalled();

  stack().push(this);
#endif
}

void FpeMonitor::ensureSignalHandlerInstalled() {
  auto &state = globalState();
  if (state.isSignalHandlerInstalled) {
    return;
  }

  std::lock_guard lock{state.mutex};

  struct sigaction action;
  action.sa_sigaction = &signalHandler;
  action.sa_flags = SA_SIGINFO;
  sigaction(SIGFPE, &action, nullptr);

  state.isSignalHandlerInstalled = true;
}

void FpeMonitor::disable() {
#if defined(__APPLE__)
  std::cerr << "FPE monitoring currently not supported on Apple" << std::endl;
#else
  std::feclearexcept(m_excepts);
  fedisableexcept(m_excepts);
  stack().pop();
#endif

  std::cout << std::bitset<32>(m_encountered) << std::endl;
}

std::stack<FpeMonitor *> &FpeMonitor::stack() {
  static thread_local std::stack<FpeMonitor *> monitors;
  return monitors;
}

FpeMonitor::GlobalState &FpeMonitor::globalState() {
  static GlobalState state{};
  return state;
}

bool FpeMonitor::encountered(FpeType type) const {
  // std::cout << "check 0x" << std::bitset<32>(1 <<
  // static_cast<uint32_t>(type))
  // << std::endl;
  // std::cout << "store 0x" << std::bitset<32>(m_encountered) << std::endl;
  uint32_t mask = 1 << static_cast<uint32_t>(type);
  return ACTS_CHECK_BIT(m_encountered, mask);
}
}  // namespace Acts
