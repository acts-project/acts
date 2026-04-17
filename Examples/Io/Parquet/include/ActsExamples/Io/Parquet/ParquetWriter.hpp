// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/IWriter.hpp"

#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>

namespace ActsExamples {

/// Writer for a set of Parquet files, one per collection, using the nested
/// layout: one row per event, each column a @c list<T> of the per-object
/// values for that event.
///
/// For each configured collection, the writer reads a 1-row
/// @c std::shared_ptr<arrow::Table> from the whiteboard under the collection
/// key, prepends an @c event_id column, and buffers the row until
/// @c eventsPerRowGroup events have accumulated. The buffer is then
/// concatenated and flushed as a single Parquet row group to the configured
/// output path. The file is opened lazily on the first event so the schema
/// can be inferred from the data.
///
/// Writes are serialized with a mutex.
class ParquetWriter final : public IWriter {
 public:
  struct Config {
    /// Base output directory. Relative @c collections paths are resolved
    /// against this directory; absolute paths are used as-is.
    std::filesystem::path outputDir;

    /// Collections to write, keyed by whiteboard collection name. The value
    /// is the output Parquet file path. A relative path is interpreted
    /// relative to @c outputDir; an absolute path is used directly. No two
    /// collections may resolve to the same output path.
    std::unordered_map<std::string, std::filesystem::path> collections;

    /// Number of events to accumulate before flushing as one Parquet row
    /// group. Smaller values improve per-event seek granularity at the cost
    /// of larger file metadata and worse compression; larger values compress
    /// better but force more events through the scanner per read at input
    /// time.
    ///
    /// Any unflushed events are written in @c finalize().
    std::size_t eventsPerRowGroup = 1000;
  };

  /// Construct the writer.
  ///
  /// @param config The configuration.
  /// @param logger The logger to use.
  /// @throws std::invalid_argument if the configuration is invalid.
  ParquetWriter(const Config& config,
                std::unique_ptr<const Acts::Logger> logger);

  /// Convenience overload: build a default logger at @p level.
  ParquetWriter(const Config& config, Acts::Logging::Level level);

  ~ParquetWriter() override;

  std::string name() const override;

  const Config& config() const;

  ProcessCode write(const AlgorithmContext& ctx) override;

  ProcessCode finalize() override;

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  class Impl;

  std::unique_ptr<const Acts::Logger> m_logger;
  std::unique_ptr<Impl> m_impl;
};

}  // namespace ActsExamples
