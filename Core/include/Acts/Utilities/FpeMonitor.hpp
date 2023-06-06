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
#include <cstddef>
#include <limits>
#include <memory>
#include <memory_resource>
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

    Result merged(const Result &with) const;
    void merge(const Result &with);

    bool encountered(FpeType type) const;
    unsigned int count(FpeType type) const;

    const std::pmr::vector<FpeInfo> &stackTraces() const;
    unsigned int numStackTraces() const;

    void deduplicate();

    void summary(
        std::ostream &os,
        std::size_t depth = std::numeric_limits<std::size_t>::max()) const;

    Result(std::pmr::memory_resource &mem = *std::pmr::new_delete_resource());

    operator bool() const { return !m_stracktraces.empty(); }

   private:
    std::pmr::vector<FpeInfo> m_stracktraces;
    std::array<unsigned int, 32> m_counts{};

    friend FpeMonitor;
  };

  FpeMonitor(std::pmr::memory_resource &mem = *std::pmr::new_delete_resource());
  explicit FpeMonitor(int excepts, std::pmr::memory_resource &mem =
                                       *std::pmr::new_delete_resource());
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

  int m_excepts = 0;
  std::size_t m_stackLimit = 1000;

  Result m_result;
};

}  // namespace Acts
