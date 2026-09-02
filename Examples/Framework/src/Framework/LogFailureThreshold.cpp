// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Framework/LogFailureThreshold.hpp"

#include <atomic>

namespace ActsExamples {

namespace {

// constinit so a read during static initialisation reports MAX. Atomic because
// a write may race a read, although the intended use is a single write at
// startup.
constinit std::atomic<Acts::Logging::Level> s_logFailureThreshold{
    Acts::Logging::Level::MAX};

}  // namespace

Acts::Logging::Level getLogFailureThreshold() {
  return s_logFailureThreshold.load(std::memory_order_relaxed);
}

void setLogFailureThreshold(Acts::Logging::Level level) {
  s_logFailureThreshold.store(level, std::memory_order_relaxed);
}

}  // namespace ActsExamples
