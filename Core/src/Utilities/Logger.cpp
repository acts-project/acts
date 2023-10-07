// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Logger.hpp"

#include <cassert>

namespace Acts {

namespace Logging {

#if defined(ACTS_ENABLE_LOG_FAILURE_THRESHOLD) and \
    not defined(ACTS_LOG_FAILURE_THRESHOLD)
namespace {
Level& getFailureThresholdMutable() {
  static Level _level = []() {
    Level level = Level::MAX;

    const char* envvar = std::getenv("ACTS_LOG_FAILURE_THRESHOLD");
    if (envvar == nullptr) {
      return level;
    }
    std::string slevel = envvar;
    if (slevel == "VERBOSE") {
      level = std::min(level, Level::VERBOSE);
    } else if (slevel == "DEBUG") {
      level = std::min(level, Level::DEBUG);
    } else if (slevel == "INFO") {
      level = std::min(level, Level::INFO);
    } else if (slevel == "WARNING") {
      level = std::min(level, Level::WARNING);
    } else if (slevel == "ERROR") {
      level = std::min(level, Level::ERROR);
    } else if (slevel == "FATAL") {
      level = std::min(level, Level::FATAL);
    } else {
      std::cerr << "ACTS_LOG_FAILURE_THRESHOLD environment variable is set to "
                   "unknown value: "
                << slevel << std::endl;
    }
    return level;
  }();

  return _level;
}
}  // namespace

Level getFailureThreshold() {
  return getFailureThresholdMutable();
}

void setFailureThreshold(Level level) {
  getFailureThresholdMutable() = level;
}

#else

void setFailureThreshold(Level /*lvl*/) {
  throw std::logic_error{
      "Compile-time log failure threshold defined (ACTS_LOG_FAILURE_THRESHOLD "
      "is set or ACTS_ENABLE_LOG_FAILURE_THRESHOLD is OFF), unable to "
      "override. See "
      "https://acts.readthedocs.io/en/latest/core/"
      "logging.html#logging-thresholds"};
}

#endif

namespace {
class NeverFilterPolicy final : public OutputFilterPolicy {
 public:
  ~NeverFilterPolicy() override = default;

  bool doPrint(const Level& /*lvl*/) const override { return false; }

  Level level() const override { return Level::MAX; }

  std::unique_ptr<OutputFilterPolicy> clone(Level /*level*/) const override {
    return std::make_unique<NeverFilterPolicy>();
  }
};

std::unique_ptr<const Logger> makeDummyLogger() {
  using namespace Logging;
  auto output = std::make_unique<DefaultPrintPolicy>(&std::cout);
  auto print = std::make_unique<NeverFilterPolicy>();
  return std::make_unique<const Logger>(std::move(output), std::move(print));
}

}  // namespace
}  // namespace Logging

std::unique_ptr<const Logger> getDefaultLogger(const std::string& name,
                                               const Logging::Level& lvl,
                                               std::ostream* log_stream) {
  using namespace Logging;
  auto output = std::make_unique<LevelOutputDecorator>(
      std::make_unique<NamedOutputDecorator>(
          std::make_unique<TimedOutputDecorator>(
              std::make_unique<DefaultPrintPolicy>(log_stream)),
          name));
  auto print = std::make_unique<DefaultFilterPolicy>(lvl);
  return std::make_unique<const Logger>(std::move(output), std::move(print));
}

const Logger& getDummyLogger() {
  static const std::unique_ptr<const Logger> logger =
      Logging::makeDummyLogger();

  return *logger;
}
}  // namespace Acts
