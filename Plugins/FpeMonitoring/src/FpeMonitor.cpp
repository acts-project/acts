// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/FpeMonitoring/FpeMonitor.hpp"

#include <algorithm>
#include <cfenv>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <memory>
#include <mutex>
#include <vector>

#include <boost/stacktrace/frame.hpp>
#include <boost/stacktrace/safe_dump_to.hpp>
#include <boost/stacktrace/stacktrace.hpp>
#include <boost/stacktrace/stacktrace_fwd.hpp>

#include "FpeMonitorPlatform.hpp"

namespace ActsPlugins {

namespace {

bool canMergeFpeInfo(const FpeMonitor::Result::FpeInfo &existing, FpeType type,
                     std::uintptr_t location,
                     const boost::stacktrace::stacktrace &st) {
  if (existing.type != type) {
    return false;
  }

  if (location != 0 && existing.location != 0) {
    return location == existing.location;
  }

  const auto &existingFrame = *existing.st->begin();
  const auto &candidateFrame = *st.begin();
  return boost::stacktrace::hash_value(existingFrame) ==
         boost::stacktrace::hash_value(candidateFrame);
}
}  // namespace

FpeMonitor::Result::FpeInfo::~FpeInfo() = default;

FpeMonitor::Result::FpeInfo::FpeInfo(
    std::size_t countIn, FpeType typeIn,
    std::shared_ptr<const boost::stacktrace::stacktrace> stIn,
    std::uintptr_t locationIn)
    : count{countIn}, type{typeIn}, st{std::move(stIn)}, location{locationIn} {}

FpeMonitor::Result FpeMonitor::Result::merged(const Result &with) const {
  Result result{};

  for (unsigned int i = 0; i < m_counts.size(); i++) {
    result.m_counts[i] = m_counts[i] + with.m_counts[i];
  }

  std::ranges::copy(with.m_stackTraces,
                    std::back_inserter(result.m_stackTraces));
  std::ranges::copy(m_stackTraces, std::back_inserter(result.m_stackTraces));

  result.deduplicate();

  return result;
}

void FpeMonitor::Result::merge(const Result &with) {
  for (unsigned int i = 0; i < m_counts.size(); i++) {
    m_counts[i] = m_counts[i] + with.m_counts[i];
  }

  std::ranges::copy(with.m_stackTraces, std::back_inserter(m_stackTraces));

  deduplicate();
}

void FpeMonitor::Result::add(FpeType type, void *stackPtr,
                             std::size_t bufferSize, std::uintptr_t location) {
  auto st = std::make_unique<boost::stacktrace::stacktrace>(
      boost::stacktrace::stacktrace::from_dump(stackPtr, bufferSize));

  for (auto &el : m_stackTraces) {
    if (canMergeFpeInfo(el, type, location, *st)) {
      el.count += 1;
      return;
    }
  }

  m_stackTraces.emplace_back(1, type, std::move(st), location);
}

bool FpeMonitor::Result::contains(const FpeInfo &info) const {
  return std::ranges::any_of(m_stackTraces, [&](const FpeInfo &el) {
    return canMergeFpeInfo(el, info.type, info.location, *info.st);
  });
}

FpeMonitor::Result &FpeMonitor::result() {
  consumeRecorded();
  return m_result;
}

void FpeMonitor::consumeRecorded() {
  if (m_recorded.empty()) {
    return;
  }

  for (auto [type, stackPtr, remaining, location] : m_recorded) {
    m_result.add(type, stackPtr, remaining, location);
  }

  m_buffer.reset();
  m_recorded.clear();
}

unsigned int FpeMonitor::Result::count(FpeType type) const {
  return m_counts.at(static_cast<std::uint32_t>(type));
}

unsigned int FpeMonitor::Result::numStackTraces() const {
  return m_stackTraces.size();
}

const std::vector<FpeMonitor::Result::FpeInfo> &
FpeMonitor::Result::stackTraces() const {
  return m_stackTraces;
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
  for (const auto &info : stackTraces()) {
    os << "- " << info.type << ": (" << info.count << " times)\n";

    os << stackTraceToString(*info.st, depth);
  }
  os << std::endl;
}

void FpeMonitor::Result::deduplicate() {
  std::vector<FpeInfo> copy = std::move(m_stackTraces);
  m_stackTraces.clear();
  m_stackTraces.reserve(copy.size());

  for (auto &info : copy) {
    const auto mergeTarget =
        std::ranges::find_if(m_stackTraces, [&](const FpeInfo &existing) {
          return canMergeFpeInfo(existing, info.type, info.location, *info.st);
        });

    if (mergeTarget != m_stackTraces.end()) {
      mergeTarget->count += info.count;
    } else {
      m_stackTraces.push_back(info);
    }
  }
}

FpeMonitor::FpeMonitor()
    : m_excepts{FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW} {
  enable();
}

FpeMonitor::FpeMonitor(int excepts) : m_excepts(excepts) {
  enable();
}

FpeMonitor::~FpeMonitor() {
  disable();
}

void FpeMonitor::signalHandler(int signal, siginfo_t *si, void *ctx) {
  if (stack().empty()) {
    return;
  }

  FpeMonitor &fpe = *stack().top();
  auto type = detail::decodeFpeType(signal, si, ctx);
  if (!type.has_value()) {
    if (detail::shouldFailFastOnUnknownSignal()) {
      // Must use _Exit: async-signal-safe. std::terminate is not (calls
      // terminate handler). std::abort raises SIGABRT from within a handler.
      std::_Exit(EXIT_FAILURE);
    }
    return;
  }
  fpe.m_result.m_counts.at(static_cast<std::uint32_t>(*type))++;
  std::uintptr_t location =
      si != nullptr ? reinterpret_cast<std::uintptr_t>(si->si_addr) : 0;

  try {
    auto [buffer, remaining] = fpe.m_buffer.next();
    using NativeFramePtr = boost::stacktrace::frame::native_frame_ptr_t;
    std::size_t stored =
        detail::captureStackFromSignalContext(ctx, buffer, remaining);
    if (stored == 0) {
      std::size_t depth = boost::stacktrace::safe_dump_to(
          detail::safeDumpSkipFrames(), buffer, remaining);
      stored = depth * sizeof(NativeFramePtr);
    }
    if (stored > 0) {
      fpe.m_buffer.pushOffset(stored);  // record how much storage was consumed
      fpe.m_recorded.emplace_back(
          *type, buffer, stored,
          location);  // record consumed stack dump and trap location
    }

  } catch (const std::bad_alloc &e) {
    std::cout << "Unable to collect stack trace due to memory limit"
              << std::endl;
  }

  detail::maskTrapsInSignalContext(ctx, *type);
}

void FpeMonitor::enable() {
  if (!detail::isRuntimeSupported()) {
    return;
  }
  ensureSignalHandlerInstalled();

  // clear pending exceptions so they don't immediately fire
  detail::clearPendingExceptions(m_excepts);

  if (!stack().empty()) {
    // unset previous except state
    detail::disableExceptions(stack().top()->m_excepts);
  }
  // apply this stack
  detail::enableExceptions(m_excepts);

  stack().push(this);
}

void FpeMonitor::rearm() {
  consumeRecorded();
  if (!detail::isRuntimeSupported()) {
    return;
  }
  detail::clearPendingExceptions(m_excepts);
  detail::enableExceptions(m_excepts);
}

void FpeMonitor::ensureSignalHandlerInstalled() {
  auto &state = globalState();
  if (state.isSignalHandlerInstalled) {
    return;
  }

  std::lock_guard lock{state.mutex};
  if (state.isSignalHandlerInstalled) {
    return;
  }

  detail::installSignalHandlers(&signalHandler);

  state.isSignalHandlerInstalled = true;
}

void FpeMonitor::disable() {
  if (!detail::isRuntimeSupported()) {
    return;
  }
  detail::clearPendingExceptions(m_excepts);
  assert(!stack().empty() && "FPE stack shouldn't be empty at this point");
  stack().pop();
  // disable excepts we enabled here
  detail::disableExceptions(m_excepts);
  if (!stack().empty()) {
    // restore excepts from next stack element
    detail::clearPendingExceptions(stack().top()->m_excepts);
    detail::enableExceptions(stack().top()->m_excepts);
  }
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

std::string FpeMonitor::stackTraceToString(
    const boost::stacktrace::stacktrace &st, std::size_t depth) {
  return boost::stacktrace::detail::to_string(st.as_vector().data(),
                                              std::min(depth, st.size()));
}

std::string FpeMonitor::getSourceLocation(
    const boost::stacktrace::frame &frame) {
  return frame.source_file() + ":" + std::to_string(frame.source_line());
}

bool FpeMonitor::canSymbolize() {
#if defined(BOOST_STACKTRACE_USE_NOOP)
  return false;
#else
  return true;
#endif
}

bool FpeMonitor::isSupported() {
  return detail::isRuntimeSupported();
}

}  // namespace ActsPlugins
