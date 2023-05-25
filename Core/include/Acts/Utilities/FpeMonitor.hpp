// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Helpers.hpp"

#include <atomic>
#include <csignal>
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

class FpeMonitor {
 public:
  FpeMonitor();
  explicit FpeMonitor(int excepts);
  ~FpeMonitor();

  bool encountered(FpeType type) const;

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

  uint32_t m_encountered = 0;
};

}  // namespace Acts
