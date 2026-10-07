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

#include <cerrno>
#include <cstdlib>
#include <exception>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/wait.h>
#include <unistd.h>
#endif

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

#if defined(__unix__) || defined(__APPLE__)
// Each timer/branch runs in its own child so termination cannot mask failures
// in another destructor or stop the test runner.
void checkTermination(const Acts::Logger& logger) {
  for (int kind : {0, 1, 2}) {
    for (bool unwind : {false, true}) {
      const auto child = fork();
      BOOST_REQUIRE(child >= 0);
      if (child == 0) {
        std::set_terminate([] { std::_Exit(42); });
        try {
          auto failIfRequested = [&] {
            if (unwind) {
              throw std::logic_error("original exception");
            }
          };
          if (kind == 0) {
            Acts::ScopedTimer timer("simple", logger);
            failIfRequested();
          } else {
            Acts::AveragingScopedTimer timer("averaging", logger);
            if (kind == 2) {
              auto sample = timer.sample();
            }
            failIfRequested();
          }
        } catch (...) {
          std::_Exit(43);  // A swallowed logging failure during unwinding.
        }
        std::_Exit(44);  // A swallowed logging failure on normal destruction.
      }
      int status = 0;
      pid_t waited;
      do {
        waited = waitpid(child, &status, 0);
      } while (waited == -1 && errno == EINTR);
      BOOST_REQUIRE_EQUAL(waited, child);
      BOOST_REQUIRE(WIFEXITED(status));
      BOOST_CHECK_EQUAL(WEXITSTATUS(status), 42);
    }
  }
}
#endif

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

#if defined(__unix__) || defined(__APPLE__)
BOOST_AUTO_TEST_CASE(ThrowingFilterTerminates) {
  std::ostringstream output;
  Acts::Logger logger(
      std::make_unique<Acts::Logging::DefaultPrintPolicy>(&output),
      std::make_unique<ThrowingFilterPolicy>());
  BOOST_CHECK_THROW(logger.doPrint(Acts::Logging::INFO), int);
  checkTermination(logger);
  BOOST_CHECK(output.str().empty());
}

BOOST_AUTO_TEST_CASE(ThrowingOutputTerminates) {
  ThrowingStreamBuffer buffer;
  std::ostream output(&buffer);
  output.exceptions(std::ios::badbit | std::ios::failbit);
  auto logger = Acts::getDefaultLogger("timers", Acts::Logging::INFO, &output);
  BOOST_CHECK_THROW(logger->log(Acts::Logging::INFO, "test"), std::exception);
  output.clear();
  checkTermination(*logger);
}

#ifdef ACTS_ENABLE_LOG_FAILURE_THRESHOLD
BOOST_AUTO_TEST_CASE(FailureThresholdTerminates) {
  Acts::Logging::ScopedFailureThreshold threshold(Acts::Logging::INFO);
  std::ostringstream output;
  auto logger = Acts::getDefaultLogger("timers", Acts::Logging::INFO, &output);
  checkTermination(*logger);
}
#endif
#endif

BOOST_AUTO_TEST_CASE(DisabledOutput) {
  checkDestruction(Acts::getDummyLogger());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
