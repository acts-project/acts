// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Parquet/ParquetWriter.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsPlugins/Arrow/ArrowUtil.hpp"

#include <cstdint>
#include <format>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <arrow/api.h>

namespace ActsExamples {

class ParquetWriter::Impl {
 public:
  using TableHandle = ConsumeDataHandle<ActsPlugins::ArrowUtil::ArrowTable>;

  struct ShardState {
    std::filesystem::path path;
    std::mutex mutex;
    ActsPlugins::ArrowUtil::ParquetFileWriter writer;
    std::vector<std::shared_ptr<arrow::Table>> buffer;
    std::uint64_t eventsAccepted = 0;

    explicit ShardState(std::filesystem::path p)
        : path(std::move(p)), writer(path) {}
  };

  struct CollectionState {
    std::string name;
    std::filesystem::path directory;
    std::string filePrefix;
    std::shared_ptr<arrow::Schema> expectedSchema;
    std::unique_ptr<TableHandle> handle;
    std::mutex shardsMutex;
    std::map<std::uint64_t, std::unique_ptr<ShardState>> shards;
  };

  Impl(const ParquetWriter::Config& config, ParquetWriter& parent)
      : m_cfg(config) {
    if (m_cfg.collections.empty()) {
      throw std::invalid_argument("ParquetWriter: no collections configured");
    }
    if (m_cfg.eventsPerShard == 0) {
      throw std::invalid_argument("ParquetWriter: eventsPerShard must be > 0");
    }
    if (m_cfg.maxOpenShards == 0) {
      throw std::invalid_argument("ParquetWriter: maxOpenShards must be > 0");
    }
    // Sentinel: 0 means "one row group per shard".
    if (m_cfg.eventsPerRowGroup == 0) {
      m_cfg.eventsPerRowGroup = m_cfg.eventsPerShard;
    }
    if (m_cfg.eventsPerRowGroup > m_cfg.eventsPerShard) {
      throw std::invalid_argument(std::format(
          "ParquetWriter: eventsPerRowGroup ({}) must be <= eventsPerShard "
          "({})",
          m_cfg.eventsPerRowGroup, m_cfg.eventsPerShard));
    }

    std::unordered_set<std::string> seenDirs;
    for (const auto& [name, rawPath] : m_cfg.collections) {
      if (name.empty()) {
        throw std::invalid_argument("ParquetWriter: empty collection name");
      }
      if (rawPath.empty()) {
        throw std::invalid_argument(std::format(
            "ParquetWriter: empty output directory for collection '{}'", name));
      }

      std::filesystem::path resolved =
          rawPath.is_absolute() ? rawPath : m_cfg.outputDir / rawPath;
      if (!seenDirs.insert(resolved.lexically_normal().string()).second) {
        throw std::invalid_argument(
            std::format("ParquetWriter: duplicate output directory '{}'",
                        resolved.string()));
      }

      // Derive a file-name prefix from the configured per-collection
      // directory: e.g. `<outputDir>/particles.parquet/` → `particles`,
      // `<outputDir>/particles/` → `particles`. Walk past empty trailing
      // components so paths with trailing separators still decompose.
      auto trimmed = resolved;
      while (!trimmed.empty() && trimmed.filename().empty()) {
        trimmed = trimmed.parent_path();
      }
      std::string filePrefix = trimmed.filename().stem().string();
      if (filePrefix.empty()) {
        throw std::invalid_argument(std::format(
            "ParquetWriter: output directory '{}' for collection '{}' has "
            "no usable filename stem to derive a shard prefix from",
            resolved.string(), name));
      }

      auto schemaIt = m_cfg.expectedSchemas.find(name);
      if (schemaIt == m_cfg.expectedSchemas.end() || !schemaIt->second) {
        throw std::invalid_argument(std::format(
            "ParquetWriter: collection '{}' has no expected schema. Every "
            "configured collection must declare an expected schema; the "
            "writer validates each per-event table against it before "
            "stamping event_id and serialising.",
            name));
      }
      auto expected = schemaIt->second.schema();
      if (expected->GetFieldIndex(
              std::string{ActsPlugins::ArrowUtil::kEventIdColumn}) >= 0) {
        throw std::invalid_argument(std::format(
            "ParquetWriter: expected schema for '{}' must not contain "
            "event_id; the writer prepends it.",
            name));
      }

      std::filesystem::create_directories(resolved);

      auto state = std::make_unique<CollectionState>();
      state->name = name;
      state->directory = std::move(resolved);
      state->filePrefix = std::move(filePrefix);
      state->expectedSchema = std::move(expected);
      state->handle = std::make_unique<TableHandle>(&parent, name);
      state->handle->initialize(name);
      m_collectionStates.push_back(std::move(state));
    }
    for (const auto& [name, _] : m_cfg.expectedSchemas) {
      if (!m_cfg.collections.contains(name)) {
        throw std::invalid_argument(std::format(
            "ParquetWriter: expectedSchemas has entry for '{}' but no "
            "matching collection",
            name));
      }
    }
  }

  /// Concatenate the shard's buffered tables and write them as a single
  /// row group. Caller must hold @c shard->mutex.
  void flushBuffer(ShardState& shard) {
    if (shard.buffer.empty()) {
      return;
    }
    auto result = arrow::ConcatenateTables(shard.buffer);
    if (!result.ok()) {
      throw std::runtime_error(
          std::format("ParquetWriter concat for shard '{}': {}",
                      shard.path.string(), result.status().ToString()));
    }
    shard.writer.write(*result.ValueOrDie());
    shard.buffer.clear();
  }

  std::filesystem::path shardPath(const CollectionState& state,
                                  std::uint64_t shardId) const {
    // Encode the *planned* event window covered by this shard: shardId is
    // assigned by `eventNumber / eventsPerShard`, so this shard owns
    // events [startEvent, endEvent). The final shard of a job may contain
    // fewer events than the name implies; consumers should read the file
    // and trust its row count, not the filename.
    const std::uint64_t startEvent = shardId * m_cfg.eventsPerShard;
    const std::uint64_t endEvent = startEvent + m_cfg.eventsPerShard;
    return state.directory / std::format("{}_{:06}-{:06}.parquet",
                                         state.filePrefix, startEvent,
                                         endEvent);
  }

  ParquetWriter::Config m_cfg;
  std::vector<std::unique_ptr<CollectionState>> m_collectionStates;
};

ParquetWriter::ParquetWriter(const Config& config,
                             std::unique_ptr<const Acts::Logger> logger)
    : m_logger(std::move(logger)),
      m_impl(std::make_unique<Impl>(config, *this)) {}

ParquetWriter::ParquetWriter(const Config& config, Acts::Logging::Level level)
    : ParquetWriter(config, Acts::getDefaultLogger("ParquetWriter", level)) {}

ParquetWriter::~ParquetWriter() = default;

std::string ParquetWriter::name() const {
  return "ParquetWriter";
}

const ParquetWriter::Config& ParquetWriter::config() const {
  return m_impl->m_cfg;
}

ProcessCode ParquetWriter::write(const AlgorithmContext& ctx) {
  using ShardState = Impl::ShardState;

  const std::uint64_t shardId = ctx.eventNumber / m_impl->m_cfg.eventsPerShard;

  for (const auto& state : m_impl->m_collectionStates) {
    auto handle = (*state->handle)(ctx);
    if (!handle) {
      ACTS_ERROR("ParquetWriter: null table for collection " << state->name);
      return ProcessCode::ABORT;
    }
    const auto& tableSchema = *handle.table()->schema();
    if (!tableSchema.Equals(*state->expectedSchema, /*check_metadata=*/false)) {
      ACTS_ERROR("ParquetWriter: schema mismatch for collection '"
                 << state->name << "' at event " << ctx.eventNumber
                 << ".\n  expected: " << state->expectedSchema->ToString()
                 << "\n  actual:   " << tableSchema.ToString());
      return ProcessCode::ABORT;
    }
    auto stamped = ActsPlugins::ArrowUtil::withEventId(
        handle.table(), static_cast<std::uint64_t>(ctx.eventNumber));

    ShardState* shard = nullptr;
    {
      std::lock_guard<std::mutex> guard(state->shardsMutex);
      auto it = state->shards.find(shardId);
      if (it == state->shards.end()) {
        if (state->shards.size() >= m_impl->m_cfg.maxOpenShards) {
          std::string openIds;
          for (const auto& [id, _] : state->shards) {
            if (!openIds.empty()) {
              openIds += ", ";
            }
            openIds += std::to_string(id);
          }
          throw std::runtime_error(std::format(
              "ParquetWriter: collection '{}' would exceed maxOpenShards={} "
              "(currently open: [{}], requested shard: {}). This usually "
              "means a worker thread is far behind the event-id frontier.",
              state->name, m_impl->m_cfg.maxOpenShards, openIds, shardId));
        }
        auto created =
            std::make_unique<ShardState>(m_impl->shardPath(*state, shardId));
        shard = created.get();
        state->shards.emplace(shardId, std::move(created));
      } else {
        shard = it->second.get();
      }
    }

    bool full = false;
    {
      std::lock_guard<std::mutex> guard(shard->mutex);
      shard->buffer.push_back(std::move(stamped));
      shard->eventsAccepted += 1;
      full = (shard->eventsAccepted >= m_impl->m_cfg.eventsPerShard);
      if (full) {
        m_impl->flushBuffer(*shard);
        shard->writer.close();
      } else if (shard->buffer.size() >= m_impl->m_cfg.eventsPerRowGroup) {
        m_impl->flushBuffer(*shard);
      }
    }

    if (full) {
      std::lock_guard<std::mutex> guard(state->shardsMutex);
      state->shards.erase(shardId);
    }
  }

  return ProcessCode::SUCCESS;
}

ProcessCode ParquetWriter::finalize() {
  for (const auto& state : m_impl->m_collectionStates) {
    std::lock_guard<std::mutex> guard(state->shardsMutex);
    for (auto& [shardId, shard] : state->shards) {
      std::lock_guard<std::mutex> sguard(shard->mutex);
      m_impl->flushBuffer(*shard);
      shard->writer.close();
    }
    state->shards.clear();
  }
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
