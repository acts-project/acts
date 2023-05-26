// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/StackTrace.hpp"

#include <atomic>
#include <csignal>
#include <memory>
#include <mutex>
#include <stack>

namespace Acts {

enum class FpeType : uint32_t {
  INTDIV = FPE_INTDIV,
  INTOVF = FPE_INTOVF,
  FLTDIV = FPE_FLTDIV,
  FLTOVF = FPE_FLTOVF,
  FLTUND = FPE_FLTUND,
  FLTRES = FPE_FLTRES,
  FLTINV = FPE_FLTINV,
  FLTSUB = FPE_FLTSUB,
};

std::ostream &operator<<(std::ostream &os, FpeType type);

class FpeMonitor {
 public:
  struct Result {
    struct FpeInfo {
      std::size_t count;
      FpeType type;
      StackTrace st;
    };

    Result merge(const Result &with) const;
    bool encountered(FpeType type) const;
    unsigned int count(FpeType type) const;

    const std::vector<FpeInfo> &stackTraces() const;
    unsigned int numStackTraces() const;

    void deduplicate();

    void summary(std::ostream &os) const;

    Result();

   private:
    std::vector<FpeInfo> m_stracktraces;
    std::array<unsigned int, 32> m_counts{};

    friend FpeMonitor;
  };

  FpeMonitor();
  explicit FpeMonitor(int excepts);
  ~FpeMonitor();

  const Result &result() const;
  Result &result();

  void rearm() const;

  std::size_t stackLimit() const { return m_stackLimit; }
  void setStackLimit(std::size_t limit) { m_stackLimit = limit; }

 private:
  void enable();
  void disable();

  static void ensureSignalHandlerInstalled();
  static void signalHandler(int signal, siginfo_t *si, void *ctx);

  struct GlobalState {
    std::atomic_bool isSignalHandlerInstalled{false};
    std::mutex mutex{};
  };

  static std::stack<FpeMonitor *> &stack();
  static GlobalState &globalState();

  int m_excepts;
  std::size_t m_stackLimit = 1000;

  Result m_result;
};

}  // namespace Acts
