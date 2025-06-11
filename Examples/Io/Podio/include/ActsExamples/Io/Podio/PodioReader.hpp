// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/IReader.hpp"

#include <filesystem>
#include <memory>

namespace ActsExamples {

namespace detail {
class PodioReaderImpl;
}

/// This reader reads events from a PODIO file in the form of frame.
/// On it's own, it does not do anything with the frame, but is used to
/// trigger the reading of the frame from the file.
///
/// The frame is then written to the event store for use by other algorithms.
///
/// The reader supports parallel execution by opening a reader per thread.
///
/// @note This class may throw exceptions if the configuration is invalid,
///       such as when the input file doesn't exist or required parameters are
///       empty.
class PodioReader : public IReader {
 public:
  struct Config {
    /// The path to the PODIO file to read.
    std::filesystem::path inputPath;
    /// The name of the frame to write to the event store.
    std::string outputFrame = "events";
    /// The podio `category` name to read the frame from.
    std::string category = "events";
  };

  /// Construct the reader.
  ///
  /// @param config The configuration struct.
  /// @param level The logging level.
  /// @throws std::invalid_argument if the configuration is invalid
  PodioReader(const Config& config, Acts::Logging::Level level);

  /// Destruct the reader.
  ~PodioReader() override;

  /// Get the name of the reader.
  ///
  /// @return The name of the reader.
  std::string name() const final;

  /// Get the range of events available in the file.
  ///
  /// @return The range of events.
  std::pair<std::size_t, std::size_t> availableEvents() const final;

  /// Read an event from the file.
  ///
  /// @param context The algorithm context.
  /// @return The process code. Returns SUCCESS when the frame is successfully read and written.
  ProcessCode read(const ActsExamples::AlgorithmContext& context) final;

  /// Get the configuration of the reader.
  ///
  /// @return The configuration of the reader.
  const Config& config() const;

 private:
  std::unique_ptr<detail::PodioReaderImpl> m_impl;
  std::unique_ptr<const Acts::Logger> m_logger;

  /// Get the logger instance.
  ///
  /// @return The logger instance.
  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
