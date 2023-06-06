// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/FpeMonitor.hpp"

#include "Acts/Utilities/Helpers.hpp"

#include <algorithm>
#include <bitset>
#include <cfenv>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <memory>
#include <memory_resource>
#include <mutex>
#include <stdexcept>
#include <string_view>
#include <vector>

#include <fenv.h>
#include <signal.h>

#define FPU_EXCEPTION_MASK 0x3f
#define FPU_STATUS_FLAGS 0xff
#define SSE_STATUS_FLAGS FPU_EXCEPTION_MASK
#define SSE_EXCEPTION_MASK (FPU_EXCEPTION_MASK << 7)

namespace Acts {

FpeMonitor::Result::Result(std::pmr::memory_resource &mem)
    : m_stracktraces{&mem} {}

FpeMonitor::Result FpeMonitor::Result::merged(const Result &with) const {
  Result result{*m_stracktraces.get_allocator().resource()};

  for (unsigned int i = 0; i < m_counts.size(); i++) {
    result.m_counts[i] = m_counts[i] + with.m_counts[i];
  }

  std::copy(with.m_stracktraces.begin(), with.m_stracktraces.end(),
            std::back_inserter(result.m_stracktraces));
  std::copy(m_stracktraces.begin(), m_stracktraces.end(),
            std::back_inserter(result.m_stracktraces));

  result.deduplicate();

  return result;
}

void FpeMonitor::Result::merge(const Result &with) {
  for (unsigned int i = 0; i < m_counts.size(); i++) {
    m_counts[i] = m_counts[i] + with.m_counts[i];
  }

  std::copy(with.m_stracktraces.begin(), with.m_stracktraces.end(),
            std::back_inserter(m_stracktraces));

  deduplicate();
}

const FpeMonitor::Result &FpeMonitor::result() const {
  return m_result;
}

FpeMonitor::Result &FpeMonitor::result() {
  return m_result;
}

unsigned int FpeMonitor::Result::count(FpeType type) const {
  return m_counts.at(static_cast<uint32_t>(type));
}

unsigned int FpeMonitor::Result::numStackTraces() const {
  return m_stracktraces.size();
}

const std::pmr::vector<FpeMonitor::Result::FpeInfo>
    &FpeMonitor::Result::stackTraces() const {
  return m_stracktraces;
}

bool FpeMonitor::Result::encountered(FpeType type) const {
  return count(type) > 0;
}

void FpeMonitor::Result::summary(std::ostream &os, std::size_t depth) const {
  os << "FPE result summary:\n";
  static const std::vector<FpeType> types = {
      FpeType::INTDIV, FpeType::INTOVF, FpeType::FLTDIV, FpeType::FLTOVF,
      FpeType::FLTUND, FpeType::FLTRES, FpeType::FLTINV, FpeType::FLTSUB};

  for (auto type : types) {
    os << "- " << type << ": " << count(type) << "\n";
  }

  os << "\nStack traces:\n";
  for (const auto &[count, type, st] : stackTraces()) {
    os << "- " << type << ": (" << count << " times)\n" << st.toString(depth);
  }
}

void FpeMonitor::Result::deduplicate() {
  std::pmr::vector<FpeInfo> copy{m_stracktraces.get_allocator()};
  copy = std::move(m_stracktraces);
  m_stracktraces.clear();

  for (auto &info : copy) {
    // auto type = it.first;
    // auto &st = it.second;
    auto it = std::find_if(m_stracktraces.begin(), m_stracktraces.end(),
                           [&info](const FpeInfo &el) {
                             return el.type == info.type && el.st == info.st;
                           });
    if (it != m_stracktraces.end()) {
      it->count += info.count;
      continue;
    }
    m_stracktraces.push_back({info.count, info.type, std::move(info.st)});
  }
}

FpeMonitor::FpeMonitor(std::pmr::memory_resource &mem)
    : m_excepts{FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW},
      m_result{mem} {
  enable();
}

FpeMonitor::FpeMonitor(int excepts, std::pmr::memory_resource &mem)
    : m_excepts(excepts), m_result{mem} {
  enable();
}

FpeMonitor::~FpeMonitor() {
  disable();
}

void FpeMonitor::signalHandler(int /*signal*/, siginfo_t *si, void *ctx) {
  if (stack().empty()) {
    return;
  }

  FpeMonitor &fpe = *stack().top();
  fpe.m_result.m_counts.at(si->si_code)++;

  std::size_t maxDepth = static_cast<std::size_t>(-1);

  try {
    // collect stack trace skipping 2 frames, which should be the signal handler
    // and the calling facility. This might be platform specific, not sure
    if (fpe.m_result.m_stracktraces.size() < fpe.stackLimit()) {
      fpe.m_result.m_stracktraces.push_back(
          {1, static_cast<FpeType>(si->si_code),
           StackTrace(
               2, maxDepth,
               *fpe.m_result.m_stracktraces.get_allocator().resource())});
    }

  } catch (const std::bad_alloc &e) {
    std::cout << "Unable to collect stack trace due to memory limit"
              << std::endl;
  }

  // @TODO: Disable this on non-x86_64
  __uint16_t *cw = &((ucontext_t *)ctx)->uc_mcontext.fpregs->cwd;
  *cw |= FPU_EXCEPTION_MASK;

  __uint16_t *sw = &((ucontext_t *)ctx)->uc_mcontext.fpregs->swd;
  *sw &= ~FPU_STATUS_FLAGS;

  __uint32_t *mxcsr = &((ucontext_t *)ctx)->uc_mcontext.fpregs->mxcsr;
  // *mxcsr |= SSE_EXCEPTION_MASK;  // disable all SSE exceptions
  *mxcsr |= ((*mxcsr & SSE_STATUS_FLAGS) << 7);
  *mxcsr &= ~SSE_STATUS_FLAGS;  // clear all pending SSE exceptions
}

void FpeMonitor::enable() {
  // @TODO: Disable this on non-x86_64
#if defined(__APPLE__)
  std::cerr << "FPE monitoring currently not supported on Apple" << std::endl;
#else
  std::feclearexcept(m_excepts);
  feenableexcept(m_excepts);

  ensureSignalHandlerInstalled();

  stack().push(this);
#endif
}

void FpeMonitor::rearm() const {
  feenableexcept(m_excepts);
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
  // @TODO: Disable this on non-x86_64
#if defined(__APPLE__)
  std::cerr << "FPE monitoring currently not supported on Apple" << std::endl;
#else
  std::feclearexcept(m_excepts);
  // fedisableexcept(m_excepts);
  stack().pop();
#endif

  // std::cout << std::bitset<32>(m_impl->m_encountered) << std::endl;
}

std::stack<FpeMonitor *> &FpeMonitor::stack() {
  static thread_local std::stack<FpeMonitor *> monitors;
  return monitors;
}

FpeMonitor::GlobalState &FpeMonitor::globalState() {
  static GlobalState state{};
  return state;
}

std::ostream &operator<<(std::ostream &os, FpeType type) {
#define CASE(x)    \
  case FpeType::x: \
    os << #x;      \
    break;

  switch (type) {
    CASE(INTDIV)
    CASE(INTOVF)
    CASE(FLTDIV)
    CASE(FLTOVF)
    CASE(FLTUND)
    CASE(FLTRES)
    CASE(FLTINV)
    CASE(FLTSUB)
  }
#undef CASE

  return os;
}

}  // namespace Acts
