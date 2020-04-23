// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/// @file Logger_tests.cpp

#include <boost/test/unit_test.hpp>
#include <fstream>
#include <string>

#include "Acts/Utilities/Logger.hpp"

namespace Acts {
namespace Test {

using namespace Acts::Logging;

/// @cond
namespace detail {
std::unique_ptr<const Logger> create_logger(const std::string& logger_name,
                                            std::ostream* logfile,
                                            Logging::Level lvl) {
  auto output = std::make_unique<LevelOutputDecorator>(
      std::make_unique<NamedOutputDecorator>(
          std::make_unique<DefaultPrintPolicy>(logfile), logger_name));
  auto print = std::make_unique<DefaultFilterPolicy>(lvl);
  return std::make_unique<const Logger>(std::move(output), std::move(print));
}

std::string failure_msg(const std::string& expected, const std::string& found) {
  return std::string("'") + expected + "' != '" + found + "'";
}
}  // namespace detail
/// @endcond

/// @brief unit test for FATAL debug level
///
/// This test checks for the expected output when using the
/// Acts::Logging::FATAL debug level as threshold. It also tests
/// - #ACTS_LOCAL_LOGGER
/// - Acts::getDefaultLogger
BOOST_AUTO_TEST_CASE(FATAL_test) {
  std::ofstream logfile("fatal_log.txt");

  auto log = detail::create_logger("TestLogger", &logfile, FATAL);
  ACTS_LOCAL_LOGGER(std::move(log));
  ACTS_FATAL("fatal level");
  ACTS_ERROR("error level");
  ACTS_WARNING("warning level");
  ACTS_INFO("info level");
  ACTS_DEBUG("debug level");
  ACTS_VERBOSE("verbose level");
  logfile.close();

  std::vector<std::string> lines;
  lines.push_back("TestLogger     FATAL     fatal level");

  std::ifstream infile("fatal_log.txt", std::ios::in);
  size_t i = 0;
  for (std::string line; std::getline(infile, line); ++i) {
    BOOST_TEST(line == lines.at(i), detail::failure_msg(line, lines.at(i)));
  }
}

/// @brief unit test for ERROR debug level
///
/// This test checks for the expected output when using the
/// Acts::Logging::ERROR debug level as threshold. It also tests
/// - #ACTS_LOCAL_LOGGER
/// - Acts::getDefaultLogger
BOOST_AUTO_TEST_CASE(ERROR_test) {
  std::ofstream logfile("error_log.txt");

  auto log = detail::create_logger("TestLogger", &logfile, ERROR);
  ACTS_LOCAL_LOGGER(std::move(log));
  ACTS_FATAL("fatal level");
  ACTS_ERROR("error level");
  ACTS_WARNING("warning level");
  ACTS_INFO("info level");
  ACTS_DEBUG("debug level");
  ACTS_VERBOSE("verbose level");
  logfile.close();

  std::vector<std::string> lines;
  lines.push_back("TestLogger     FATAL     fatal level");
  lines.push_back("TestLogger     ERROR     error level");

  std::ifstream infile("error_log.txt", std::ios::in);
  size_t i = 0;
  for (std::string line; std::getline(infile, line); ++i) {
    BOOST_TEST(line == lines.at(i), detail::failure_msg(line, lines.at(i)));
  }
}

/// @brief unit test for WARNING debug level
///
/// This test checks for the expected output when using the
/// Acts::Logging::WARNING debug level as threshold. It also tests
/// - #ACTS_LOCAL_LOGGER
/// - Acts::getDefaultLogger
BOOST_AUTO_TEST_CASE(WARNING_test) {
  std::ofstream logfile("warning_log.txt");

  auto log = detail::create_logger("TestLogger", &logfile, WARNING);
  ACTS_LOCAL_LOGGER(std::move(log));
  ACTS_FATAL("fatal level");
  ACTS_ERROR("error level");
  ACTS_WARNING("warning level");
  ACTS_INFO("info level");
  ACTS_DEBUG("debug level");
  ACTS_VERBOSE("verbose level");
  logfile.close();

  std::vector<std::string> lines;
  lines.push_back("TestLogger     FATAL     fatal level");
  lines.push_back("TestLogger     ERROR     error level");
  lines.push_back("TestLogger     WARNING   warning level");

  std::ifstream infile("warning_log.txt", std::ios::in);
  size_t i = 0;
  for (std::string line; std::getline(infile, line); ++i) {
    BOOST_TEST(line == lines.at(i), detail::failure_msg(line, lines.at(i)));
  }
}

/// @brief unit test for INFO debug level
///
/// This test checks for the expected output when using the
/// Acts::Logging::INFO debug level as threshold. It also tests
/// - #ACTS_LOCAL_LOGGER
/// - Acts::getDefaultLogger
BOOST_AUTO_TEST_CASE(INFO_test) {
  std::ofstream logfile("info_log.txt");

  auto log = detail::create_logger("TestLogger", &logfile, INFO);
  ACTS_LOCAL_LOGGER(std::move(log));
  ACTS_FATAL("fatal level");
  ACTS_ERROR("error level");
  ACTS_WARNING("warning level");
  ACTS_INFO("info level");
  ACTS_DEBUG("debug level");
  ACTS_VERBOSE("verbose level");
  logfile.close();

  std::vector<std::string> lines;
  lines.push_back("TestLogger     FATAL     fatal level");
  lines.push_back("TestLogger     ERROR     error level");
  lines.push_back("TestLogger     WARNING   warning level");
  lines.push_back("TestLogger     INFO      info level");

  std::ifstream infile("info_log.txt", std::ios::in);
  size_t i = 0;
  for (std::string line; std::getline(infile, line); ++i) {
    BOOST_TEST(line == lines.at(i), detail::failure_msg(line, lines.at(i)));
  }
}

/// @brief unit test for DEBUG debug level
///
/// This test checks for the expected output when using the
/// Acts::Logging::DEBUG debug level as threshold. It also tests
/// - #ACTS_LOCAL_LOGGER
/// - Acts::getDefaultLogger
BOOST_AUTO_TEST_CASE(DEBUG_test) {
  std::ofstream logfile("debug_log.txt");

  auto log = detail::create_logger("TestLogger", &logfile, DEBUG);
  ACTS_LOCAL_LOGGER(std::move(log));
  ACTS_FATAL("fatal level");
  ACTS_ERROR("error level");
  ACTS_WARNING("warning level");
  ACTS_INFO("info level");
  ACTS_DEBUG("debug level");
  ACTS_VERBOSE("verbose level");
  logfile.close();

  std::vector<std::string> lines;
  lines.push_back("TestLogger     FATAL     fatal level");
  lines.push_back("TestLogger     ERROR     error level");
  lines.push_back("TestLogger     WARNING   warning level");
  lines.push_back("TestLogger     INFO      info level");
  lines.push_back("TestLogger     DEBUG     debug level");

  std::ifstream infile("debug_log.txt", std::ios::in);
  size_t i = 0;
  for (std::string line; std::getline(infile, line); ++i) {
    BOOST_TEST(line == lines.at(i), detail::failure_msg(line, lines.at(i)));
  }
}

/// @brief unit test for VERBOSE debug level
///
/// This test checks for the expected output when using the
/// Acts::Logging::VERBOSE debug level as threshold. It also tests
/// - #ACTS_LOCAL_LOGGER
/// - Acts::getDefaultLogger
BOOST_AUTO_TEST_CASE(VERBOSE_test) {
  std::ofstream logfile("verbose_log.txt");

  auto log = detail::create_logger("TestLogger", &logfile, VERBOSE);
  ACTS_LOCAL_LOGGER(std::move(log));
  ACTS_FATAL("fatal level");
  ACTS_ERROR("error level");
  ACTS_WARNING("warning level");
  ACTS_INFO("info level");
  ACTS_DEBUG("debug level");
  ACTS_VERBOSE("verbose level");
  logfile.close();

  std::vector<std::string> lines;
  lines.push_back("TestLogger     FATAL     fatal level");
  lines.push_back("TestLogger     ERROR     error level");
  lines.push_back("TestLogger     WARNING   warning level");
  lines.push_back("TestLogger     INFO      info level");
  lines.push_back("TestLogger     DEBUG     debug level");
  lines.push_back("TestLogger     VERBOSE   verbose level");

  std::ifstream infile("verbose_log.txt", std::ios::in);
  size_t i = 0;
  for (std::string line; std::getline(infile, line); ++i) {
    BOOST_TEST(line == lines.at(i), detail::failure_msg(line, lines.at(i)));
  }
}
}  // namespace Test
}  // namespace Acts
