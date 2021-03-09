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

LoggerWrapper::LoggerWrapper(const Logger& logger) : m_logger(&logger) {}

bool LoggerWrapper::doPrint(const Logging::Level& lvl) const {
  assert(m_logger != nullptr);
  return m_logger->doPrint(lvl);
}

void LoggerWrapper::log(const Logging::Level& lvl,
                        const std::string& input) const {
  assert(m_logger != nullptr);
  return m_logger->log(lvl, input);
}

const Logger& LoggerWrapper::operator()() const {
  assert(m_logger != nullptr);
  return *m_logger;
}

namespace Logging {

namespace {
class NeverFilterPolicy final : public OutputFilterPolicy {
 public:
  ~NeverFilterPolicy() override = default;

  bool doPrint(const Level& /*lvl*/) const override { return false; }
};

std::unique_ptr<const Logger> makeDummyLogger() {
  using namespace Logging;
  auto output = std::make_unique<DefaultPrintPolicy>(&std::cout);
  auto print = std::make_unique<NeverFilterPolicy>();
  return std::make_unique<const Logger>(std::move(output), std::move(print));
}

std::unique_ptr<const Logger> s_dummyLogger{makeDummyLogger()};
LoggerWrapper s_dummyLoggerWrapper{*s_dummyLogger};

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

LoggerWrapper getDummyLogger() {
  return Logging::s_dummyLoggerWrapper;
}
}  // namespace Acts
