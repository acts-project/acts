// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Logger.hpp"

#include "Acts/Utilities/Diagnostics.hpp"

#include <cstdlib>
#include <optional>
#include <string>

namespace Acts {

namespace Logging {

#ifdef ACTS_ENABLE_LOG_FAILURE_THRESHOLD

namespace {

std::optional<Level> parseFailureThreshold(const char* value) {
  const std::string level = value;
  if (level == "VERBOSE") {
    return Level::VERBOSE;
  } else if (level == "DEBUG") {
    return Level::DEBUG;
  } else if (level == "INFO") {
    return Level::INFO;
  } else if (level == "WARNING") {
    return Level::WARNING;
  } else if (level == "ERROR") {
    return Level::ERROR;
  } else if (level == "FATAL") {
    return Level::FATAL;
  } else if (level == "MAX") {
    return Level::MAX;
  }
  std::cerr << "ACTS_LOG_FAILURE_THRESHOLD is set to unknown value: " << level
            << std::endl;
  return std::nullopt;
}

std::optional<Level>& defaultFailureThresholdMutable() {
  static std::optional<Level> s_level = []() -> std::optional<Level> {
#ifdef ACTS_LOG_FAILURE_THRESHOLD
    // Seeded from the deprecated CMake option. Unlike before, this is only a
    // default: it can still be overridden at runtime.
    std::optional<Level> level = Level::ACTS_LOG_FAILURE_THRESHOLD;
#else
    std::optional<Level> level = std::nullopt;
#endif
    if (const char* envvar = std::getenv("ACTS_LOG_FAILURE_THRESHOLD");
        envvar != nullptr) {
      if (std::optional<Level> parsed = parseFailureThreshold(envvar);
          parsed.has_value()) {
        level = parsed;
      }
    }
    return level;
  }();

  return s_level;
}

}  // namespace

namespace detail {

std::optional<Level> getDefaultFailureThreshold() {
  return defaultFailureThresholdMutable();
}

void setDefaultFailureThreshold(std::optional<Level> level) {
  defaultFailureThresholdMutable() = level;
}

}  // namespace detail

#else

namespace detail {

void setDefaultFailureThreshold(std::optional<Level> /*level*/) {
  // ACTS_ENABLE_LOG_FAILURE_THRESHOLD=OFF: this build can never be armed from
  // the outside, so setting the process-wide default has no effect.
}

}  // namespace detail

#endif

ACTS_PUSH_IGNORE_DEPRECATED()

Level getFailureThreshold() {
  return detail::getDefaultFailureThreshold().value_or(Level::MAX);
}

void setFailureThreshold(Level level) {
  detail::setDefaultFailureThreshold(level);
}

ScopedFailureThreshold::~ScopedFailureThreshold() noexcept {
  try {
    detail::setDefaultFailureThreshold(m_previousLevel);
  } catch (const std::bad_alloc&) {
    // bad alloc can be thrown when initializing the global static variable
    std::cerr << "Failed to reset log failure threshold (bad_alloc)"
              << std::endl;
    std::terminate();
  }
}

ACTS_POP_IGNORE_DEPRECATED()

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

class DummyPrintPolicy final : public OutputPrintPolicy {
 public:
  void flush(const Level& /*lvl*/, const std::string& /*input*/) override {}

  const std::string& name() const override {
    const static std::string s_name = "Dummy";
    return s_name;
  }

  std::unique_ptr<OutputPrintPolicy> clone(
      const std::string& /*name*/) const override {
    return std::make_unique<DummyPrintPolicy>();
  }
};

std::unique_ptr<const Logger> makeDummyLogger() {
  using namespace Logging;
  auto output = std::make_unique<DummyPrintPolicy>();
  auto print = std::make_unique<NeverFilterPolicy>();
  return std::make_unique<const Logger>(std::move(output), std::move(print),
                                        Level::MAX);
}

}  // namespace
}  // namespace Logging

std::unique_ptr<const Logger> getDefaultLogger(
    const std::string& name, const Logging::Level& lvl,
    std::ostream* log_stream, std::optional<Logging::Level> failureThreshold) {
  using namespace Logging;
  auto output = std::make_unique<LevelOutputDecorator>(
      std::make_unique<NamedOutputDecorator>(
          std::make_unique<TimedOutputDecorator>(
              std::make_unique<DefaultPrintPolicy>(log_stream)),
          name));
  auto print = std::make_unique<DefaultFilterPolicy>(lvl);
  return std::make_unique<const Logger>(std::move(output), std::move(print),
                                        failureThreshold);
}

const Logger& getDummyLogger() {
  static const std::unique_ptr<const Logger> logger =
      Logging::makeDummyLogger();

  return *logger;
}
}  // namespace Acts
