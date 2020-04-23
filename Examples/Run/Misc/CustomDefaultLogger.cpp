// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <algorithm>

#include "Acts/Utilities/Logger.hpp"

namespace Acts {

namespace Logging {

/// @brief mirror debug message
///
/// This is just a fun example decorator which mirrors the debug message.
class MirrorOutputDecorator final : public OutputDecorator {
 public:
  /// @brief constructor
  ///
  /// @param [in] wrappee  output print policy object to be wrapped
  /// @param [in] maxWidth maximum width of field used for name
  MirrorOutputDecorator(std::unique_ptr<OutputPrintPolicy> wrappee,
                        unsigned int maxWidth = 180)
      : OutputDecorator(std::move(wrappee)), m_maxWidth(maxWidth) {}

  /// @brief flush the debug message to the destination stream
  ///
  /// @param [in] lvl   debug level of debug message
  /// @param [in] input text of debug message
  ///
  /// This function inverts the given string and flushes it to the right.
  void flush(const Level& lvl, const std::ostringstream& input) override {
    std::ostringstream os;
    std::string text = input.str();
    std::reverse(text.begin(), text.end());
    os << std::right << std::setw(m_maxWidth) << text;
    OutputDecorator::flush(lvl, os);
  }

 private:
  /// maximum width of field for printing the name
  unsigned int m_maxWidth;
};

}  // namespace Logging

/// @brief alternative implementation of default debug output logger
///
/// @param [in] name       name of the logger instance
/// @param [in] lvl        debug threshold level
/// @param [in] log_stream output stream used for printing debug messages
///
/// This function returns a pointer to an alternative Logger instance to
/// the default Acts logger. This instance prints the log output mirrored
/// from right to left.
///
/// @return pointer to logging instance
std::unique_ptr<const Logger> getDefaultLogger(const std::string& name,
                                               const Logging::Level& lvl,
                                               std::ostream* log_stream) {
  using namespace Logging;
  auto output = std::make_unique<LevelOutputDecorator>(
      std::make_unique<NamedOutputDecorator>(
          std::make_unique<TimedOutputDecorator>(
              std::make_unique<MirrorOutputDecorator>(
                  std::make_unique<DefaultPrintPolicy>(log_stream))),
          name));
  auto print = std::make_unique<DefaultFilterPolicy>(lvl);
  return std::make_unique<const Logger>(std::move(output), std::move(print));
}

}  // namespace Acts
