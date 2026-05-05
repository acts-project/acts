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

  struct PendingBuffer {
    std::vector<std::shared_ptr<arrow::Table>> tables;
    std::int64_t rows = 0;
  };

  Impl(const ParquetWriter::Config& config, ParquetWriter& parent)
      : m_cfg(config) {
    if (m_cfg.collections.empty()) {
      throw std::invalid_argument("ParquetWriter: no collections configured");
    }
    if (m_cfg.eventsPerRowGroup == 0) {
      throw std::invalid_argument(
          "ParquetWriter: eventsPerRowGroup must be > 0");
    }

    std::unordered_set<std::string> seenPaths;
    for (const auto& [name, rawPath] : m_cfg.collections) {
      if (name.empty()) {
        throw std::invalid_argument("ParquetWriter: empty collection name");
      }
      if (rawPath.empty()) {
        throw std::invalid_argument(std::format(
            "ParquetWriter: empty output path for collection '{}'", name));
      }

      std::filesystem::path resolved =
          rawPath.is_absolute() ? rawPath : m_cfg.outputDir / rawPath;
      if (!seenPaths.insert(resolved.lexically_normal().string()).second) {
        throw std::invalid_argument(std::format(
            "ParquetWriter: duplicate output path '{}'", resolved.string()));
      }

      if (resolved.has_parent_path()) {
        std::filesystem::create_directories(resolved.parent_path());
      }

      auto handle = std::make_unique<TableHandle>(&parent, name);
      handle->initialize(name);
      m_handles.push_back(std::move(handle));

      m_writers.emplace(
          name, std::make_unique<ActsPlugins::ArrowUtil::ParquetFileWriter>(
                    resolved));
      m_pending.emplace(name, PendingBuffer{});
    }
  }

  /// Concatenate and flush the pending buffer for @p name as one row group.
  /// Does nothing if the buffer is empty.
  void flush(const std::string& name) {
    auto& buf = m_pending.at(name);
    if (buf.tables.empty()) {
      return;
    }
    auto result = arrow::ConcatenateTables(buf.tables);
    if (!result.ok()) {
      throw std::runtime_error(std::format("ParquetWriter concat for '{}': {}",
                                           name, result.status().ToString()));
    }
    m_writers.at(name)->write(*result.ValueOrDie());
    buf.tables.clear();
    buf.rows = 0;
  }

  ParquetWriter::Config m_cfg;
  std::vector<std::unique_ptr<TableHandle>> m_handles;
  std::unordered_map<std::string,
                     std::unique_ptr<ActsPlugins::ArrowUtil::ParquetFileWriter>>
      m_writers;
  std::unordered_map<std::string, PendingBuffer> m_pending;
  std::mutex m_writeMutex;
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
  std::unique_lock<std::mutex> guard(m_impl->m_writeMutex);

  for (const auto& handle : m_impl->m_handles) {
    auto table = (*handle)(ctx);
    if (!table) {
      ACTS_ERROR("ParquetWriter: null table for collection " << handle->name());
      return ProcessCode::ABORT;
    }

    auto stamped = ActsPlugins::ArrowUtil::withEventId(
        table, static_cast<std::uint64_t>(ctx.eventNumber));

    auto& buf = m_impl->m_pending.at(handle->name());
    buf.rows += stamped->num_rows();
    buf.tables.push_back(std::move(stamped));

    if (static_cast<std::size_t>(buf.rows) >= m_impl->m_cfg.eventsPerRowGroup) {
      m_impl->flush(handle->name());
    }
  }

  return ProcessCode::SUCCESS;
}

ProcessCode ParquetWriter::finalize() {
  std::unique_lock<std::mutex> guard(m_impl->m_writeMutex);
  for (const auto& [name, _] : m_impl->m_cfg.collections) {
    m_impl->flush(name);
  }
  for (auto& [name, writer] : m_impl->m_writers) {
    writer->close();
  }
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
