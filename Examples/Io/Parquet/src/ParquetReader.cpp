// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Parquet/ParquetReader.hpp"

#include "Acts/Utilities/ScopedTimer.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsPlugins/Arrow/ArrowUtil.hpp"

#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include <arrow/api.h>

namespace ActsExamples {

class ParquetReader::Impl {
 public:
  explicit Impl(ParquetReader::Config cfg, ParquetReader& parent)
      : m_cfg(std::move(cfg)) {
    if (m_cfg.collections.empty()) {
      throw std::invalid_argument("ParquetReader: no collections configured");
    }
    std::unordered_set<std::string> seenPaths;
    for (const auto& [name, rawPath] : m_cfg.collections) {
      if (name.empty()) {
        throw std::invalid_argument("ParquetReader: empty collection name");
      }
      if (rawPath.empty()) {
        throw std::invalid_argument(std::format(
            "ParquetReader: empty input path for collection '{}'", name));
      }

      std::filesystem::path resolved =
          rawPath.is_absolute() ? rawPath : m_cfg.inputDir / rawPath;
      if (!seenPaths.insert(resolved.lexically_normal().string()).second) {
        throw std::invalid_argument(std::format(
            "ParquetReader: duplicate input path '{}'", resolved.string()));
      }
    }

    std::int64_t referenceEvents = -1;
    std::string referenceName;

    for (const auto& [name, rawPath] : m_cfg.collections) {
      std::filesystem::path path =
          rawPath.is_absolute() ? rawPath : m_cfg.inputDir / rawPath;
      if (!std::filesystem::exists(path)) {
        throw std::invalid_argument(
            std::format("ParquetReader: missing file '{}'", path.string()));
      }

      const auto events = ActsPlugins::ArrowUtil::numRowsInFile(path);
      if (referenceEvents < 0) {
        referenceEvents = events;
        referenceName = name;
      } else if (events != referenceEvents) {
        throw std::invalid_argument(std::format(
            "ParquetReader: event count mismatch across collections. '{}' has "
            "{} events, '{}' has {}. All collections must describe the same "
            "events.",
            referenceName, referenceEvents, name, events));
      }

      auto table = ActsPlugins::ArrowUtil::readTable(path);

      auto handle =
          std::make_unique<WriteDataHandle<std::shared_ptr<arrow::Table>>>(
              &parent, name);
      handle->initialize(name);

      m_tables.emplace(name, std::move(table));
      m_handles.push_back(std::move(handle));
    }

    m_eventsRange = {
        0,
        referenceEvents < 0 ? 0 : static_cast<std::size_t>(referenceEvents)};
  }

  std::shared_ptr<arrow::Table> sliceForEvent(const std::string& name,
                                              std::size_t eventNumber) const {
    return ActsPlugins::ArrowUtil::sliceByEventId(
        m_tables.at(name), static_cast<std::uint64_t>(eventNumber));
  }

  ParquetReader::Config m_cfg;
  std::unordered_map<std::string, std::shared_ptr<arrow::Table>> m_tables;
  std::vector<std::unique_ptr<WriteDataHandle<std::shared_ptr<arrow::Table>>>>
      m_handles;
  std::pair<std::size_t, std::size_t> m_eventsRange{0, 0};
};

ParquetReader::ParquetReader(const Config& config,
                             std::unique_ptr<const Acts::Logger> logger)
    : m_impl(std::make_unique<Impl>(config, *this)),
      m_logger(std::move(logger)) {}

ParquetReader::ParquetReader(const Config& config, Acts::Logging::Level level)
    : ParquetReader(config, Acts::getDefaultLogger("ParquetReader", level)) {}

ParquetReader::~ParquetReader() = default;

std::string ParquetReader::name() const {
  return "ParquetReader";
}

std::pair<std::size_t, std::size_t> ParquetReader::availableEvents() const {
  return m_impl->m_eventsRange;
}

ProcessCode ParquetReader::read(const AlgorithmContext& context) {
  Acts::ScopedTimer timer("Reading Parquet inputs", logger(),
                          Acts::Logging::DEBUG);

  for (const auto& handle : m_impl->m_handles) {
    auto slice = m_impl->sliceForEvent(handle->name(), context.eventNumber);
    (*handle)(context, std::move(slice));
  }

  return ProcessCode::SUCCESS;
}

const ParquetReader::Config& ParquetReader::config() const {
  return m_impl->m_cfg;
}

}  // namespace ActsExamples
