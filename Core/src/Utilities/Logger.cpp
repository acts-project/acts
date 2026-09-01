// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Logger.hpp"

#include "Acts/Utilities/Diagnostics.hpp"

#include <atomic>
#include <cstdlib>
#include <optional>
#include <string>

namespace Acts {

namespace Logging {

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

}  // namespace

namespace {

// Read when a logger is built, not when a message is logged, so this is off
// every hot path. constinit so that a read during static initialisation --
// before the seeder below has run -- reports MAX rather than an indeterminate
// value. Atomic because a write from one thread may race a read from another,
// even though the intended use is a single write at startup.
constinit std::atomic<Level> s_defaultFailureThreshold{Level::MAX};

}  // namespace

namespace detail {

Level getDefaultFailureThreshold() {
  return s_defaultFailureThreshold.load(std::memory_order_relaxed);
}

void setDefaultFailureThreshold(Level level) {
  s_defaultFailureThreshold.store(level, std::memory_order_relaxed);
}

}  // namespace detail

namespace {

// Seeds the default from the environment during dynamic initialisation. Given
// the earliest available init priority so that it runs ahead of ordinary
// dynamic initialisers that might log.
struct DefaultFailureThresholdSeeder {
  DefaultFailureThresholdSeeder() {
    const char* envvar = std::getenv("ACTS_LOG_FAILURE_THRESHOLD");
    if (envvar == nullptr) {
      return;
    }
    if (std::optional<Level> parsed = parseFailureThreshold(envvar);
        parsed.has_value()) {
      detail::setDefaultFailureThreshold(*parsed);
    }
  }
};

#if defined(__GNUC__) || defined(__clang__)
[[gnu::init_priority(101)]]
#endif
const DefaultFailureThresholdSeeder s_defaultFailureThresholdSeeder;

}  // namespace

ACTS_PUSH_IGNORE_DEPRECATED()

Level getFailureThreshold() {
  return detail::getDefaultFailureThreshold();
}

void setFailureThreshold(Level level) {
  detail::setDefaultFailureThreshold(level);
}

ScopedFailureThreshold::~ScopedFailureThreshold() noexcept {
  detail::setDefaultFailureThreshold(m_previousLevel);
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
  // The process-wide default is resolved here, once, and copied into the
  // logger. Nothing reads it again afterwards.
  return std::make_unique<const Logger>(
      std::move(output), std::move(print),
      failureThreshold.value_or(Logging::detail::getDefaultFailureThreshold()));
}

const Logger& getDummyLogger() {
  static const std::unique_ptr<const Logger> logger =
      Logging::makeDummyLogger();

  return *logger;
}
}  // namespace Acts
