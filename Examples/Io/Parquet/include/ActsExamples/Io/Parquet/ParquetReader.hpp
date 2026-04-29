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
#include "ActsPlugins/Arrow/Export.hpp"

#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>

namespace ActsExamples {

/// Reader for a set of Parquet files, one per collection, in the nested
/// layout: one row per event, each column a @c list<T> of the per-object
/// values for that event.
///
/// For each configured collection, the reader opens the configured Parquet
/// file path, pre-loads the full table, and on each @c read() call filters out
/// the row with the matching @c event_id and parks it on the whiteboard under
/// the collection name. The @c event_id column is stripped on the way out.
///
/// The event count is taken from the Parquet footer (@c num_rows). All input
/// files must agree on their row count.
class ACTS_ARROW_EXPORT ParquetReader : public IReader {
 public:
  struct Config {
    /// Base input directory. Relative @c collections paths are resolved
    /// against this directory; absolute paths are used as-is.
    std::filesystem::path inputDir;

    /// Collections to read, keyed by whiteboard collection name. The value is
    /// the input Parquet file path. A relative path is interpreted relative to
    /// @c inputDir; an absolute path is used directly. No two collections may
    /// resolve to the same input path.
    std::unordered_map<std::string, std::filesystem::path> collections;
  };

  /// Construct the reader.
  ///
  /// @param config The configuration.
  /// @param logger The logger to use.
  /// @throws std::invalid_argument if the configuration is invalid or a file
  ///         is missing.
  ParquetReader(const Config& config,
                std::unique_ptr<const Acts::Logger> logger);

  /// Convenience overload: build a default logger at @p level.
  ParquetReader(const Config& config, Acts::Logging::Level level);

  ~ParquetReader() override;

  std::string name() const final;

  std::pair<std::size_t, std::size_t> availableEvents() const final;

  ProcessCode read(const AlgorithmContext& context) final;

  const Config& config() const;

 private:
  class Impl;

  std::unique_ptr<Impl> m_impl;
  std::unique_ptr<const Acts::Logger> m_logger;

  const Acts::Logger& logger() const { return *m_logger; }
};

}  // namespace ActsExamples
