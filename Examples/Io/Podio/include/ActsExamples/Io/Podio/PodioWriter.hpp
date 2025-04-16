// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IWriter.hpp"

#include <memory>
#include <string>

namespace ActsExamples {

namespace detail {
class PodioWriterImpl;
}  // namespace detail

/// This writer writes events to a PODIO file in the form of frame.
///
/// The writer supports parallel execution by serializing the writes to the
/// file.
/// The writer is configured with a @c podio::Frame name to read from the white
/// board. If empty, it will create a new frame.
///
/// The writer is also configured with a list of collection names to write to
/// the file. The collections must be present in the event store and be of type
/// @c podio::CollectionBase.
///
/// @note The writer uses a mutex to ensure thread safety when writing to the file.
class PodioWriter final : public IWriter {
 public:
  struct Config {
    /// The path to the output file.
    std::string outputPath;

    /// Retrieve a @c podio::Frame from the event store using this name.
    /// @note If not set, a new frame will be created.
    std::optional<std::string> inputFrame = std::nullopt;

    /// The podio `category` name to write the frame to.
    /// This is used to organize frames in the output file.
    /// @note This field must not be empty.
    std::string category;

    /// The collection names to write to the output file.
    /// @note Collection names must not be empty and must be unique.
    /// The collections must be present in the event store.
    std::vector<std::string> collections;
  };

  /// Construct the writer.
  ///
  /// @param config The configuration struct.
  /// @param level The logging level.
  /// @throw std::invalid_argument if category is empty or if collection names are empty or duplicate
  PodioWriter(const Config& config, Acts::Logging::Level level);

  /// Destruct the writer.
  ~PodioWriter() override;

  /// Get the name of the writer.
  ///
  /// @return The name of the writer.
  std::string name() const override;

  /// Readonly access to the config
  ///
  /// @return The configuration of the writer.
  const Config& config() const;

  /// Finalize the writer.
  /// This closes the output file and releases any resources.
  ///
  /// @return The process code.
  ProcessCode finalize() override;

  /// Write an event to the output file.
  ///
  /// @param ctx The algorithm context.
  /// @return The process code.
  /// @return ProcessCode::ABORT if any collection is not initialized
  ProcessCode write(const AlgorithmContext& ctx) override;

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  std::unique_ptr<const Acts::Logger> m_logger;

  std::unique_ptr<detail::PodioWriterImpl> m_impl;
};

}  // namespace ActsExamples
