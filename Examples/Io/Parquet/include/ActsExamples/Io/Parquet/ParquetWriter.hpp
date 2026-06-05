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
#include "ActsPlugins/Arrow/ArrowUtil.hpp"
#include "ActsPlugins/Arrow/Export.hpp"

#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>

namespace ActsExamples {

/// Writer for a sharded Parquet @c arrow::dataset, one directory per
/// collection. Events are routed to shard files by
/// @c event_id / eventsPerShard, so each shard owns a contiguous,
/// disjoint event-id range. Each shard is written as a single Parquet
/// file with one row group; file-level min/max statistics over
/// @c event_id are then tight by construction, which lets the reader
/// prune to a single fragment per per-event lookup.
///
/// Shard files inside a collection's directory are named
/// `<prefix>_<startEvent>-<endEvent>.parquet`, with the prefix taken
/// from the collection directory's stem (e.g. `particles.parquet/` →
/// `particles_000000-001000.parquet` for the first shard of an
/// `eventsPerShard=1000` job) and the bounds reflecting the *planned*
/// half-open event window the shard owns.
///
/// For each configured collection, the writer reads a 1-row
/// @c std::shared_ptr<arrow::Table> from the whiteboard under the
/// collection key, prepends an @c event_id column, buffers it under
/// the appropriate shard, and flushes the shard once it has received
/// @c eventsPerShard events.
///
/// At most @c maxOpenShards shards are kept open simultaneously per
/// collection. Under the Sequencer's sliding-window scheduling, the
/// in-flight event set spans at most ~thread-count event ids and so
/// touches at most one or two shards; this is a hard tripwire that
/// throws if it is ever exceeded.
class ACTS_ARROW_EXPORT ParquetWriter final : public IWriter {
 public:
  struct Config {
    /// Base output directory. Relative @c collections paths are resolved
    /// against this directory; absolute paths are used as-is.
    std::filesystem::path outputDir;

    /// Collections to write, keyed by whiteboard collection name. The value
    /// is the output directory (NOT a file path) into which shard files are
    /// written. A relative path is interpreted relative to @c outputDir; an
    /// absolute path is used directly. No two collections may resolve to
    /// the same output directory.
    std::unordered_map<std::string, std::filesystem::path> collections;

    /// Expected schema per collection — same keys as @c collections, must
    /// cover every key. Each per-event @c arrow::Table read off the
    /// whiteboard is checked for exact equality against the declared
    /// schema before stamping @c event_id and serialising. The schemas
    /// here MUST NOT include the @c event_id column — the writer
    /// prepends it. Symmetric with
    /// @c ParquetReader::Config::expectedSchemas: any pair of writer +
    /// reader for a collection should reference the same handle.
    std::unordered_map<std::string, ActsPlugins::ArrowUtil::ArrowSchemaHandle>
        expectedSchemas;

    /// Number of events per shard file. Drives file-level granularity
    /// (the read-side pruning unit).
    std::size_t eventsPerShard = 1000;

    /// Number of events per row group within a shard file. Bounds the
    /// in-memory write buffer: the writer holds up to this many 1-row
    /// tables per active shard before serialising them as a single row
    /// group and resetting the buffer. Lower this for collections with
    /// large per-event payloads (e.g. SimHits with millions of values
    /// per event) to keep peak memory tractable. Must satisfy
    /// @c 0 < eventsPerRowGroup <= eventsPerShard.
    ///
    /// @c 0 means "use the same value as @c eventsPerShard", which
    /// gives one row group per file. That's the simplest layout but
    /// pins peak memory to a full shard's worth of buffered data.
    std::size_t eventsPerRowGroup = 0;

    /// Hard cap on simultaneously open shards per collection. Defense in
    /// depth: with the Sequencer's sliding-window scheduling the working
    /// set is ~1-2 shards, so this should never be hit. If it is, the
    /// writer throws — silently producing files with overlapping
    /// event-id ranges would degrade the reader's pruning to a full
    /// scan.
    std::size_t maxOpenShards = 10;
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
