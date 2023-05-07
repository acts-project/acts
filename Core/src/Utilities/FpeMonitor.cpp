// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/FpeMonitor.hpp"

#include "Acts/Utilities/Helpers.hpp"

#include <cfenv>
#include <csignal>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string_view>

#include <boost/stacktrace.hpp>
#include <fenv.h>
#include <signal.h>

namespace Acts {
namespace {

void handle_fpe(int except) {
  auto check = [except](int mask) { return ACTS_CHECK_BIT(except, mask); };

  if (check(FE_OVERFLOW)) {
    std::cout << "FE_OVERFLOW" << std::endl;
  }

  if (check(FE_DIVBYZERO)) {
    std::cout << "FE_DIVBYZERO" << std::endl;
  }

  if (check(FE_INVALID)) {
    std::cout << "FE_INVALID" << std::endl;
  }

  if (check(FE_UNDERFLOW)) {
    std::cout << "FE_UNDERFLOW" << std::endl;
  }

  if (check(FE_INEXACT)) {
    std::cout << "FE_INEXACT" << std::endl;
  }

  std::cout << boost::stacktrace::stacktrace();

  std::abort();
}

}  // namespace

FpeMonitor::FpeMonitor()
    : FpeMonitor{FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW} {}

FpeMonitor::FpeMonitor(int excepts) : m_excepts(excepts) {
  enable(m_excepts);
}

FpeMonitor::~FpeMonitor() {
  disable(m_excepts);
}

void FpeMonitor::enable(int excepts) {
#if defined(__APPLE__)
  std::cerr << "FPE monitoring currently not supported on Apple" << std::endl;
  (void)excepts;
  (void)handle_fpe;
#else
  std::feclearexcept(FE_ALL_EXCEPT);
  feenableexcept(excepts);

  std::signal(SIGFPE, handle_fpe);
#endif
}
void FpeMonitor::disable(int excepts) {
#if defined(__APPLE__)
  std::cerr << "FPE monitoring currently not supported on Apple" << std::endl;
  (void)excepts;
#else
  fedisableexcept(excepts);
  std::signal(SIGFPE, SIG_DFL);
#endif
}

}  // namespace Acts
