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
#include "ActsPlugins/Arrow/ArrowUtil.hpp"
#include "ActsPlugins/Arrow/Export.hpp"

#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>

namespace ActsExamples {

/// Reader for a sharded Parquet @c arrow::dataset, one directory per
/// collection. Mirror of @c ParquetWriter: each collection is read from
/// a directory of @c shard_*.parquet files, each holding a contiguous
/// disjoint range of event ids in an @c event_id column.
///
/// Per-event lookups push a filter @c event_id == N into the dataset
/// scanner, which prunes to a single fragment and a single row group via
/// Parquet footer/page-index statistics. The matching row's other
/// columns are late-materialized via the page offset index, so the read
/// cost per event is roughly the size of one event's payload.
///
/// The total event count is the sum of per-fragment row counts (footer
/// reads only). All collections must agree on their total row count.
///
/// Every configured collection must declare an expected schema. The
/// scanner projects each fragment to that schema: columns the fragment
/// lacks are materialized as nulls (added-column schema evolution), and
/// columns the fragment has but the schema doesn't are dropped. This
/// makes the reader's output deterministic — independent of which
/// optional columns happened to be written into the on-disk shards —
/// and gives downstream input converters a contract they can rely on.
class ACTS_ARROW_EXPORT ParquetReader : public IReader {
 public:
  struct Config {
    /// Base input directory. Relative @c collections paths are resolved
    /// against this directory; absolute paths are used as-is.
    std::filesystem::path inputDir;

    /// Collections to read, keyed by whiteboard collection name. The
    /// value is the input directory containing the collection's shard
    /// files. A relative path is interpreted relative to @c inputDir;
    /// an absolute path is used directly. No two collections may
    /// resolve to the same input directory.
    std::unordered_map<std::string, std::filesystem::path> collections;

    /// Expected schema per collection — same keys as @c collections,
    /// must cover every key. The scanner uses this as the dataset's
    /// target schema (sans @c event_id, which the reader manages
    /// internally), so missing columns become nulls and extra columns
    /// are dropped. The handles should wrap the schemas that the
    /// corresponding output converter / @c ArrowSchemas helper
    /// produces. Wrapped in @c ArrowSchemaHandle so the same map type
    /// is usable from Python (where exposing @c arrow::Schema directly
    /// would re-introduce the symbol-collision risk the arrow-plugin
    /// isolation is built to prevent).
    std::unordered_map<std::string, ActsPlugins::ArrowUtil::ArrowSchemaHandle>
        expectedSchemas;
  };

  /// Construct the reader.
  ///
  /// @param config The configuration.
  /// @param logger The logger to use.
  /// @throws std::invalid_argument if the configuration is invalid or a
  ///         directory is missing.
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
