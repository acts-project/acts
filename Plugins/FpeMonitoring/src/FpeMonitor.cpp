// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/FpeMonitoring/FpeMonitor.hpp"

#include "Acts/Utilities/Helpers.hpp"

#include <algorithm>
#include <bitset>
#include <cfenv>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <vector>

#include <boost/stacktrace/frame.hpp>
#include <boost/stacktrace/safe_dump_to.hpp>
#include <boost/stacktrace/stacktrace.hpp>
#include <boost/stacktrace/stacktrace_fwd.hpp>

#define FPU_EXCEPTION_MASK 0x3f
#define FPU_STATUS_FLAGS 0xff
#define SSE_STATUS_FLAGS FPU_EXCEPTION_MASK
#define SSE_EXCEPTION_MASK (FPU_EXCEPTION_MASK << 7)

namespace ActsPlugins {

namespace {

#if (defined(__linux__) && defined(__x86_64__)) || \
    (defined(__APPLE__) && (defined(__x86_64__) || defined(__arm64__)))
constexpr bool kFpeRuntimeSupported = true;
#else
// Keep helper boundaries architecture-oriented so Linux aarch64 can be added
// by implementing decode + context masking paths, without touching call sites.
constexpr bool kFpeRuntimeSupported = false;
#endif

bool areFpesEquivalent(
    std::pair<FpeType, const boost::stacktrace::stacktrace &> lhs,
    std::pair<FpeType, const boost::stacktrace::stacktrace &> rhs) {
  const auto &fl = *lhs.second.begin();
  const auto &fr = *rhs.second.begin();
  return lhs.first == rhs.first && (boost::stacktrace::hash_value(fl) ==
                                    boost::stacktrace::hash_value(fr));
}

std::optional<FpeType> fpeTypeFromSiCode(int siCode) {
  switch (siCode) {
    case FPE_INTDIV:
      return FpeType::INTDIV;
    case FPE_INTOVF:
      return FpeType::INTOVF;
    case FPE_FLTDIV:
      return FpeType::FLTDIV;
    case FPE_FLTOVF:
      return FpeType::FLTOVF;
    case FPE_FLTUND:
      return FpeType::FLTUND;
    case FPE_FLTRES:
      return FpeType::FLTRES;
    case FPE_FLTINV:
      return FpeType::FLTINV;
    case FPE_FLTSUB:
      return FpeType::FLTSUB;
    default:
      return std::nullopt;
  }
}

#if defined(__APPLE__) && defined(__arm64__)

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

  // The low ESR bits encode IEEE FP exception classes on Darwin arm64.
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

#endif

std::optional<FpeType> decodeFpeType(int signal, siginfo_t *si, void *ctx) {
  if (signal == SIGFPE && si != nullptr) {
    return fpeTypeFromSiCode(si->si_code);
  }

#if defined(__APPLE__) && defined(__arm64__)
  if (signal == SIGILL && ctx != nullptr) {
    auto *uc = static_cast<ucontext_t *>(ctx);
    return fpeTypeFromDarwinArm64Esr(uc->uc_mcontext->__es.__esr);
  }
#endif

  return std::nullopt;
}

int exceptMaskForType(FpeType type) {
  switch (type) {
    case FpeType::INTDIV:
    case FpeType::FLTDIV:
      return FE_DIVBYZERO;
    case FpeType::INTOVF:
    case FpeType::FLTOVF:
      return FE_OVERFLOW;
    case FpeType::FLTUND:
      return FE_UNDERFLOW;
    case FpeType::FLTRES:
      return FE_INEXACT;
    case FpeType::FLTINV:
    case FpeType::FLTSUB:
      return FE_INVALID;
    default:
      return 0;
  }
}

void clearPendingExceptions(int excepts) { std::feclearexcept(excepts); }

void enableExceptions(int excepts) {
#if defined(__linux__) && defined(__x86_64__)
  feenableexcept(excepts);
#elif defined(__APPLE__) && defined(__x86_64__)
  fenv_t env{};
  if (fegetenv(&env) != 0) {
    return;
  }
  env.__control &= ~static_cast<unsigned short>(excepts);
  env.__mxcsr &= ~(static_cast<unsigned int>(excepts) << 7u);
  env.__status &= ~static_cast<unsigned short>(FE_ALL_EXCEPT);
  env.__mxcsr &= ~static_cast<unsigned int>(FE_ALL_EXCEPT);
  fesetenv(&env);
#elif defined(__APPLE__) && defined(__arm64__)
  fenv_t env{};
  if (fegetenv(&env) != 0) {
    return;
  }
  env.__fpcr |= static_cast<unsigned long long>(darwinArm64TrapMask(excepts));
  env.__fpsr &= ~static_cast<unsigned long long>(FE_ALL_EXCEPT);
  fesetenv(&env);
#else
  static_cast<void>(excepts);
#endif
}

void disableExceptions(int excepts) {
#if defined(__linux__) && defined(__x86_64__)
  fedisableexcept(excepts);
#elif defined(__APPLE__) && defined(__x86_64__)
  fenv_t env{};
  if (fegetenv(&env) != 0) {
    return;
  }
  env.__control |= static_cast<unsigned short>(excepts);
  env.__mxcsr |= (static_cast<unsigned int>(excepts) << 7u);
  fesetenv(&env);
#elif defined(__APPLE__) && defined(__arm64__)
  fenv_t env{};
  if (fegetenv(&env) != 0) {
    return;
  }
  env.__fpcr &= ~static_cast<unsigned long long>(darwinArm64TrapMask(excepts));
  fesetenv(&env);
#else
  static_cast<void>(excepts);
#endif
}

void maskTrapsInSignalContext(void *ctx, FpeType type) {
  const int excepts = exceptMaskForType(type);
#if defined(__linux__) && defined(__x86_64__)
  auto *uc = static_cast<ucontext_t *>(ctx);
  __uint16_t *cw = &uc->uc_mcontext.fpregs->cwd;
  *cw |= FPU_EXCEPTION_MASK;

  __uint16_t *sw = &uc->uc_mcontext.fpregs->swd;
  *sw &= ~FPU_STATUS_FLAGS;

  __uint32_t *mxcsr = &uc->uc_mcontext.fpregs->mxcsr;
  *mxcsr |= ((*mxcsr & SSE_STATUS_FLAGS) << 7);
  *mxcsr &= ~SSE_STATUS_FLAGS;
#elif defined(__APPLE__) && defined(__x86_64__)
  auto *uc = static_cast<ucontext_t *>(ctx);
  uc->uc_mcontext->__fs.__fpu_fcw |=
      static_cast<unsigned short>(excepts);
  uc->uc_mcontext->__fs.__fpu_fsw &= ~static_cast<unsigned short>(FE_ALL_EXCEPT);
  uc->uc_mcontext->__fs.__fpu_mxcsr |=
      (static_cast<unsigned int>(excepts) << 7u);
  uc->uc_mcontext->__fs.__fpu_mxcsr &= ~static_cast<unsigned int>(FE_ALL_EXCEPT);
#elif defined(__APPLE__) && defined(__arm64__)
  auto *uc = static_cast<ucontext_t *>(ctx);
  uc->uc_mcontext->__ns.__fpcr &=
      ~static_cast<std::uint32_t>(darwinArm64TrapMask(excepts));
  uc->uc_mcontext->__ns.__fpsr &= ~static_cast<std::uint32_t>(FE_ALL_EXCEPT);
#else
  static_cast<void>(ctx);
  static_cast<void>(excepts);
#endif
}

std::size_t captureStackFromSignalContext(void *ctx, void *buffer,
                                          std::size_t bufferBytes) {
  using NativeFramePtr = boost::stacktrace::frame::native_frame_ptr_t;
  auto *frames = static_cast<NativeFramePtr *>(buffer);
  const std::size_t maxFrames = bufferBytes / sizeof(NativeFramePtr);
  std::size_t count = 0;

  if (ctx == nullptr || maxFrames == 0) {
    return 0;
  }

#if defined(__APPLE__) && defined(__arm64__)
  auto *uc = static_cast<ucontext_t *>(ctx);
  const std::uintptr_t sp =
      __darwin_arm_thread_state64_get_sp(uc->uc_mcontext->__ss);
  std::uintptr_t fp = __darwin_arm_thread_state64_get_fp(uc->uc_mcontext->__ss);
  const std::uintptr_t pc =
      __darwin_arm_thread_state64_get_pc(uc->uc_mcontext->__ss);

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
    const auto *record = reinterpret_cast<const FrameRecord *>(fp);
    const std::uintptr_t prevFp = record->prevFp;
    const std::uintptr_t lr = record->returnAddress;
    push(lr);

    if (prevFp <= fp || !inStackWindow(prevFp)) {
      break;
    }
    fp = prevFp;
  }
#elif defined(__APPLE__) && defined(__x86_64__)
  auto *uc = static_cast<ucontext_t *>(ctx);
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
    const auto *record = reinterpret_cast<const FrameRecord *>(fp);
    const std::uintptr_t prevFp = record->prevFp;
    const std::uintptr_t ra = record->returnAddress;
    push(ra);

    if (prevFp <= fp || !inStackWindow(prevFp)) {
      break;
    }
    fp = prevFp;
  }
#else
  static_cast<void>(ctx);
  static_cast<void>(frames);
  static_cast<void>(maxFrames);
#endif

  return count * sizeof(NativeFramePtr);
}
}  // namespace

FpeMonitor::Result::FpeInfo::~FpeInfo() = default;

FpeMonitor::Result::FpeInfo::FpeInfo(
    std::size_t countIn, FpeType typeIn,
    std::shared_ptr<const boost::stacktrace::stacktrace> stIn)
    : count{countIn}, type{typeIn}, st{std::move(stIn)} {}

FpeMonitor::Result FpeMonitor::Result::merged(const Result &with) const {
  Result result{};

  for (unsigned int i = 0; i < m_counts.size(); i++) {
    result.m_counts[i] = m_counts[i] + with.m_counts[i];
  }

  std::copy(with.m_stackTraces.begin(), with.m_stackTraces.end(),
            std::back_inserter(result.m_stackTraces));
  std::copy(with.m_locations.begin(), with.m_locations.end(),
            std::back_inserter(result.m_locations));
  std::copy(m_stackTraces.begin(), m_stackTraces.end(),
            std::back_inserter(result.m_stackTraces));
  std::copy(m_locations.begin(), m_locations.end(),
            std::back_inserter(result.m_locations));

  result.deduplicate();

  return result;
}

void FpeMonitor::Result::merge(const Result &with) {
  for (unsigned int i = 0; i < m_counts.size(); i++) {
    m_counts[i] = m_counts[i] + with.m_counts[i];
  }

  std::copy(with.m_stackTraces.begin(), with.m_stackTraces.end(),
            std::back_inserter(m_stackTraces));
  std::copy(with.m_locations.begin(), with.m_locations.end(),
            std::back_inserter(m_locations));

  deduplicate();
}

void FpeMonitor::Result::add(FpeType type, void *stackPtr,
                             std::size_t bufferSize, std::uintptr_t location) {
  auto st = std::make_unique<boost::stacktrace::stacktrace>(
      boost::stacktrace::stacktrace::from_dump(stackPtr, bufferSize));

  for (std::size_t i = 0; i < m_stackTraces.size(); ++i) {
    auto &el = m_stackTraces[i];
    if (el.type != type) {
      continue;
    }

    if (location != 0 && m_locations[i] != 0) {
      if (location == m_locations[i]) {
        el.count += 1;
        return;
      }
      continue;
    }

    if (areFpesEquivalent({el.type, *el.st}, {type, *st})) {
      el.count += 1;
      return;
    }
  }

  m_stackTraces.push_back({1, type, std::move(st)});
  m_locations.push_back(location);
}

bool FpeMonitor::Result::contains(
    FpeType type, const boost::stacktrace::stacktrace &st) const {
  return std::ranges::any_of(m_stackTraces, [&](const FpeInfo &el) {
    return areFpesEquivalent({el.type, *el.st}, {type, st});
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
  for (const auto &[count, type, st] : stackTraces()) {
    os << "- " << type << ": (" << count << " times)\n";

    os << stackTraceToString(*st, depth);
  }
  os << std::endl;
}

void FpeMonitor::Result::deduplicate() {
  std::vector<FpeInfo> copy = std::move(m_stackTraces);
  std::vector<std::uintptr_t> copyLocations = std::move(m_locations);
  m_stackTraces.clear();
  m_locations.clear();

  for (std::size_t i = 0; i < copy.size(); ++i) {
    auto &info = copy[i];
    const std::uintptr_t location = copyLocations[i];

    bool merged = false;
    for (std::size_t j = 0; j < m_stackTraces.size(); ++j) {
      auto &existing = m_stackTraces[j];
      if (existing.type != info.type) {
        continue;
      }

      if (location != 0 && m_locations[j] != 0) {
        if (location == m_locations[j]) {
          existing.count += info.count;
          merged = true;
          break;
        }
        continue;
      }

      if (areFpesEquivalent({existing.type, *existing.st},
                            {info.type, *info.st})) {
        existing.count += info.count;
        merged = true;
        break;
      }
    }

    if (merged) {
      continue;
    }
    m_stackTraces.push_back({info.count, info.type, std::move(info.st)});
    m_locations.push_back(location);
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
  auto type = decodeFpeType(signal, si, ctx);
  if (!type.has_value()) {
    // For unknown SIGILL causes on Darwin arm64, returning without changing the
    // context can re-trigger the signal forever. Fail fast instead.
#if defined(__APPLE__) && defined(__arm64__)
    std::_Exit(EXIT_FAILURE);
#else
    return;
#endif
  }
  fpe.m_result.m_counts.at(static_cast<std::uint32_t>(*type))++;
  std::uintptr_t location =
      si != nullptr ? reinterpret_cast<std::uintptr_t>(si->si_addr) : 0;

  try {
    auto [buffer, remaining] = fpe.m_buffer.next();
    using NativeFramePtr = boost::stacktrace::frame::native_frame_ptr_t;
    std::size_t stored = 0;
#if defined(__APPLE__) && (defined(__arm64__) || defined(__x86_64__))
    // On Darwin, unwind from the interrupted context so masks can match frames
    // between the fault site and callers.
    stored = captureStackFromSignalContext(ctx, buffer, remaining);
    if (stored == 0) {
      std::size_t depth = boost::stacktrace::safe_dump_to(1, buffer, remaining);
      stored = depth * sizeof(NativeFramePtr);
    }
#else
    std::size_t depth = boost::stacktrace::safe_dump_to(2, buffer, remaining);
    stored = depth * sizeof(NativeFramePtr);
#endif
    if (stored > 0) {
      fpe.m_buffer.pushOffset(stored);  // record how much storage was consumed
      fpe.m_recorded.emplace_back(
          *type, buffer,
          stored, location);  // record consumed stack dump and trap location
    }

  } catch (const std::bad_alloc &e) {
    std::cout << "Unable to collect stack trace due to memory limit"
              << std::endl;
  }

  maskTrapsInSignalContext(ctx, *type);
}

void FpeMonitor::enable() {
#if (defined(__linux__) && defined(__x86_64__)) || \
    (defined(__APPLE__) && (defined(__x86_64__) || defined(__arm64__)))
  ensureSignalHandlerInstalled();

  // clear pending exceptions so they don't immediately fire
  clearPendingExceptions(m_excepts);

  if (!stack().empty()) {
    // unset previous except state
    disableExceptions(stack().top()->m_excepts);
  }
  // apply this stack
  enableExceptions(m_excepts);

  stack().push(this);
#else
  static_cast<void>(m_excepts);
#endif
}

void FpeMonitor::rearm() {
  consumeRecorded();
#if (defined(__linux__) && defined(__x86_64__)) || \
    (defined(__APPLE__) && (defined(__x86_64__) || defined(__arm64__)))
  clearPendingExceptions(m_excepts);
  enableExceptions(m_excepts);
#endif
}

void FpeMonitor::ensureSignalHandlerInstalled() {
  auto &state = globalState();
  if (state.isSignalHandlerInstalled) {
    return;
  }

  std::lock_guard lock{state.mutex};

  struct sigaction action{};
  action.sa_sigaction = &signalHandler;
  action.sa_flags = SA_SIGINFO;
  sigaction(SIGFPE, &action, nullptr);
#if defined(__APPLE__) && defined(__arm64__)
  sigaction(SIGILL, &action, nullptr);
#endif

  state.isSignalHandlerInstalled = true;
}

void FpeMonitor::disable() {
#if (defined(__linux__) && defined(__x86_64__)) || \
    (defined(__APPLE__) && (defined(__x86_64__) || defined(__arm64__)))
  clearPendingExceptions(m_excepts);
  assert(!stack().empty() && "FPE stack shouldn't be empty at this point");
  stack().pop();
  // disable excepts we enabled here
  disableExceptions(m_excepts);
  if (!stack().empty()) {
    // restore excepts from next stack element
    clearPendingExceptions(stack().top()->m_excepts);
    enableExceptions(stack().top()->m_excepts);
  }
#endif
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

bool FpeMonitor::isSupported() { return kFpeRuntimeSupported; }

}  // namespace ActsPlugins
