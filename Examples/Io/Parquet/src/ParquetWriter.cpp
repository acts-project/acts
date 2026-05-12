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
  using TableHandle = ConsumeDataHandle<std::shared_ptr<arrow::Table>>;

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

      std::filesystem::create_directories(resolved);

      auto state = std::make_unique<CollectionState>();
      state->name = name;
      state->directory = std::move(resolved);
      state->handle = std::make_unique<TableHandle>(&parent, name);
      state->handle->initialize(name);
      m_states.push_back(std::move(state));
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

  std::filesystem::path shardPath(const std::filesystem::path& directory,
                                  std::uint64_t shardId) const {
    // Derive a file-name prefix from the configured per-collection
    // directory: e.g. `<outputDir>/particles.parquet/` → `particles`,
    // `<outputDir>/particles/` → `particles`. Falls back to `shard` if the
    // directory has no usable stem (root dir, trailing slash, etc.).
    std::string prefix = directory.filename().stem().string();
    if (prefix.empty()) {
      prefix = "shard";
    }
    // Encode the *planned* event window covered by this shard: shardId is
    // assigned by `eventNumber / eventsPerShard`, so this shard owns
    // events [startEvent, endEvent). The final shard of a job may contain
    // fewer events than the name implies; consumers should read the file
    // and trust its row count, not the filename.
    const std::uint64_t startEvent = shardId * m_cfg.eventsPerShard;
    const std::uint64_t endEvent = startEvent + m_cfg.eventsPerShard;
    return directory /
           std::format("{}_{:06}-{:06}.parquet", prefix, startEvent, endEvent);
  }

  ParquetWriter::Config m_cfg;
  std::vector<std::unique_ptr<CollectionState>> m_states;
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

  for (const auto& state : m_impl->m_states) {
    auto table = (*state->handle)(ctx);
    if (!table) {
      ACTS_ERROR("ParquetWriter: null table for collection " << state->name);
      return ProcessCode::ABORT;
    }
    auto stamped = ActsPlugins::ArrowUtil::withEventId(
        table, static_cast<std::uint64_t>(ctx.eventNumber));

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
        auto created = std::make_unique<ShardState>(
            m_impl->shardPath(state->directory, shardId));
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
  for (const auto& state : m_impl->m_states) {
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
