// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace Acts;
using namespace Acts::Logging;

namespace ActsTests {

namespace {
/// Stand-in for an experiment-supplied print policy, which never carried the
/// failure-threshold check itself.
class CountingPrintPolicy final : public OutputPrintPolicy {
 public:
  void flush(const Level& /*lvl*/, const std::string& /*input*/) override {
    ++m_count;
  }
  const std::string& name() const override {
    static const std::string s_name = "Counting";
    return s_name;
  }
  std::unique_ptr<OutputPrintPolicy> clone(
      const std::string& /*name*/) const override {
    return std::make_unique<CountingPrintPolicy>();
  }
  std::size_t count() const { return m_count; }

 private:
  std::size_t m_count = 0;
};
}  // namespace

/// @cond
namespace detail {
std::unique_ptr<const Logger> create_logger(const std::string& logger_name,
                                            std::ostream* logfile,
                                            Logging::Level lvl) {
  auto output = std::make_unique<LevelOutputDecorator>(
      std::make_unique<NamedOutputDecorator>(
          std::make_unique<DefaultPrintPolicy>(logfile), logger_name, 30));
  auto print = std::make_unique<DefaultFilterPolicy>(lvl);
  // Disarmed: this test is about output formatting and must behave the same
  // whether or not the job as a whole is armed.
  return std::make_unique<const Logger>(std::move(output), std::move(print),
                                        Logging::Level::MAX);
}

}  // namespace detail
/// @endcond

/// @brief unit test for a certain debug level
///
/// This test checks for the expected output when using the
/// specified debug level as threshold. It also tests
/// - #ACTS_LOCAL_LOGGER
/// - getDefaultLogger
void debug_level_test(const char* output_file, Logging::Level lvl) {
  // Logs will go to this file
  std::ofstream logfile(output_file);

  auto test = [&](std::unique_ptr<const Logger> log, const std::string& name) {
    // Set up local logger
    ACTS_LOCAL_LOGGER(std::move(log));

    // Test logging at a certain debug level
    auto test_logging = [](auto&& test_operation, Logging::Level /*test_lvl*/) {
      test_operation();
    };

    // Test logging at all debug levels
    test_logging([&] { ACTS_FATAL("fatal level"); }, FATAL);
    test_logging([&] { ACTS_ERROR("error level"); }, ERROR);
    test_logging([&] { ACTS_WARNING("warning level"); }, WARNING);
    test_logging([&] { ACTS_INFO("info level"); }, INFO);
    test_logging([&] { ACTS_DEBUG("debug level"); }, DEBUG);
    test_logging([&] { ACTS_VERBOSE("verbose level"); }, VERBOSE);
    logfile.close();

    std::string padded_name = name;
    padded_name.resize(30, ' ');

    // Compute expected output for current debug levels
    std::vector<std::string> lines{padded_name + "FATAL     fatal level",
                                   padded_name + "ERROR     error level",
                                   padded_name + "WARNING   warning level",
                                   padded_name + "INFO      info level",
                                   padded_name + "DEBUG     debug level",
                                   padded_name + "VERBOSE   verbose level"};
    lines.resize(static_cast<int>(Logging::Level::MAX) - static_cast<int>(lvl));

    // Check output
    std::ifstream infile(output_file, std::ios::in);
    std::size_t i = 0;
    for (std::string line; std::getline(infile, line); ++i) {
      BOOST_CHECK_EQUAL(line, lines.at(i));
    }
  };

  auto log = detail::create_logger("TestLogger", &logfile, lvl);
  BOOST_CHECK_EQUAL(log->name(), "TestLogger");
  auto copy = log->clone("TestLoggerClone");
  test(std::move(copy), "TestLoggerClone");
  BOOST_CHECK_EQUAL(log->name(), "TestLogger");

  auto copy2 = log->clone("TestLoggerClone");
  BOOST_CHECK_EQUAL(copy2->level(), log->level());

  auto copy3 = log->cloneWithSuffix("Suffix");
  BOOST_CHECK_EQUAL(log->level(), copy3->level());

  logfile = std::ofstream{output_file};  // clear output

  test(std::move(log), "TestLogger");
}

BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

/// @brief unit test for FATAL debug level
BOOST_AUTO_TEST_CASE(FATAL_test) {
  debug_level_test("fatal_log.txt", FATAL);
}

/// @brief unit test for ERROR debug level
BOOST_AUTO_TEST_CASE(ERROR_test) {
  debug_level_test("error_log.txt", ERROR);
}

/// @brief unit test for WARNING debug level
BOOST_AUTO_TEST_CASE(WARNING_test) {
  debug_level_test("warning_log.txt", WARNING);
}

/// @brief unit test for INFO debug level
BOOST_AUTO_TEST_CASE(INFO_test) {
  debug_level_test("info_log.txt", INFO);
}

/// @brief unit test for DEBUG debug level
BOOST_AUTO_TEST_CASE(DEBUG_test) {
  debug_level_test("debug_log.txt", DEBUG);
}

/// @brief unit test for VERBOSE debug level
BOOST_AUTO_TEST_CASE(VERBOSE_test) {
  debug_level_test("verbose_log.txt", VERBOSE);
}

/// @brief the failure threshold is a property of the logger
BOOST_AUTO_TEST_CASE(FailureThreshold_test) {
  std::ostringstream out;

  auto armed = Acts::getDefaultLogger("Armed", Logging::VERBOSE, &out,
                                      Logging::Level::WARNING);
  BOOST_CHECK(armed->failureThreshold() == Logging::Level::WARNING);
  BOOST_CHECK_THROW(armed->log(Logging::ERROR, "boom"),
                    Logging::ThresholdFailure);
  BOOST_CHECK_NO_THROW(armed->log(Logging::INFO, "fine"));

  // the message that fails the job is emitted before the throw
  BOOST_CHECK(out.str().find("boom") != std::string::npos);

  // clones inherit the threshold, so child loggers stay armed
  auto child = armed->cloneWithSuffix("::child");
  BOOST_CHECK(child->failureThreshold() == Logging::Level::WARNING);
  BOOST_CHECK_THROW(child->log(Logging::ERROR, "boom"),
                    Logging::ThresholdFailure);

  // an explicitly disarmed logger never throws
  auto disarmed = armed->withoutFailureThreshold();
  BOOST_CHECK_NO_THROW(disarmed->log(Logging::FATAL, "expected"));

  // a coarse filter level cannot hide a message that should fail the job
  auto quiet = Acts::getDefaultLogger("Quiet", Logging::FATAL, &out,
                                      Logging::Level::WARNING);
  BOOST_CHECK_THROW(quiet->log(Logging::ERROR, "not hidden"),
                    Logging::ThresholdFailure);

  // the check does not depend on the print policy, so an experiment-supplied
  // one is covered too
  Logger custom{std::make_unique<CountingPrintPolicy>(),
                std::make_unique<DefaultFilterPolicy>(Logging::VERBOSE),
                Logging::Level::WARNING};
  BOOST_CHECK_THROW(custom.log(Logging::ERROR, "custom"),
                    Logging::ThresholdFailure);

  // the dummy logger discards everything and never throws
  BOOST_CHECK_NO_THROW(Acts::getDummyLogger().log(Logging::FATAL, "ignored"));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
