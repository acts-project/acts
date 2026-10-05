// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ScopedTimer.hpp"

#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

namespace ActsTests {
namespace {

class ThrowingFilterPolicy final : public Acts::Logging::OutputFilterPolicy {
 public:
  bool doPrint(const Acts::Logging::Level& /*level*/) const override {
    // Use a non-standard exception to exercise the catch-all handler.
    throw 42;
  }

  Acts::Logging::Level level() const override { return Acts::Logging::INFO; }

  std::unique_ptr<Acts::Logging::OutputFilterPolicy> clone(
      Acts::Logging::Level /*level*/) const override {
    return std::make_unique<ThrowingFilterPolicy>();
  }
};

class ThrowingStreamBuffer final : public std::streambuf {
 public:
  int_type overflow(int_type /*character*/) override {
    throw std::runtime_error("timing output failed");
  }
};

// Cover the simple timer and both averaging timer output branches, including
// destruction while an unrelated exception is already being unwound.
void checkDestruction(const Acts::Logger& logger) {
  auto run = [&](bool unwind) {
    Acts::ScopedTimer timer("simple", logger);
    Acts::AveragingScopedTimer empty("empty", logger);
    Acts::AveragingScopedTimer populated("populated", logger);
    {
      auto sample = populated.sample();
    }
    if (unwind) {
      throw std::logic_error("original exception");
    }
  };

  BOOST_CHECK_NO_THROW(run(false));
  BOOST_CHECK_EXCEPTION(run(true), std::logic_error, [](const auto& error) {
    return std::string(error.what()) == "original exception";
  });
}

}  // namespace

BOOST_AUTO_TEST_SUITE(ScopedTimerTests)

BOOST_AUTO_TEST_CASE(LoggingOutput) {
  std::ostringstream output;
  auto logger = Acts::getDefaultLogger("timers", Acts::Logging::INFO, &output);
  {
    Acts::ScopedTimer timer("simple", *logger);
  }
  BOOST_CHECK(output.str().find("simple took ") != std::string::npos);
  BOOST_CHECK(output.str().find(" ms") != std::string::npos);

  output.str("");
  {
    Acts::AveragingScopedTimer timer("empty", *logger);
  }
  BOOST_CHECK(output.str().find("empty took 0 ms total (no samples)") !=
              std::string::npos);

  output.str("");
  {
    Acts::AveragingScopedTimer timer("populated", *logger);
    auto sample = timer.sample();
  }
  BOOST_CHECK(output.str().find("populated took ") != std::string::npos);
  BOOST_CHECK(output.str().find("us per sample (#1)") != std::string::npos);
}

BOOST_AUTO_TEST_CASE(ThrowingFilter) {
  std::ostringstream output;
  Acts::Logger logger(
      std::make_unique<Acts::Logging::DefaultPrintPolicy>(&output),
      std::make_unique<ThrowingFilterPolicy>());
  BOOST_CHECK_THROW(logger.doPrint(Acts::Logging::INFO), int);
  checkDestruction(logger);
  BOOST_CHECK(output.str().empty());
}

BOOST_AUTO_TEST_CASE(ThrowingOutput) {
  ThrowingStreamBuffer buffer;
  std::ostream output(&buffer);
  output.exceptions(std::ios::badbit | std::ios::failbit);
  auto logger = Acts::getDefaultLogger("timers", Acts::Logging::INFO, &output);
  BOOST_CHECK_THROW(logger->log(Acts::Logging::INFO, "test"), std::exception);
  output.clear();
  checkDestruction(*logger);
}

BOOST_AUTO_TEST_CASE(DisabledOutput) {
  checkDestruction(Acts::getDummyLogger());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
