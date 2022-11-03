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
#include <cxxabi.h>    // for demangling
#include <execinfo.h>  // for backtrace
#include <fenv.h>
#include <link.h>  // for following code in shared libraries
#include <signal.h>

namespace Acts {
namespace {

void handle_fpe(int except) {
  auto check = [except](int mask) { return ACTS_CHECK_BIT(except, mask); };

  if (check(FE_OVERFLOW)) {
    std::cout << "Floating point overflow" << std::endl;
  }

  if (check(FE_DIVBYZERO)) {
    std::cout << "Floating point divide by zero" << std::endl;
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
  std::feclearexcept(FE_ALL_EXCEPT);
  feenableexcept(m_excepts);

  std::signal(SIGFPE, handle_fpe);
}

FpeMonitor::~FpeMonitor() {
  fedisableexcept(m_excepts);
}

}  // namespace Acts
