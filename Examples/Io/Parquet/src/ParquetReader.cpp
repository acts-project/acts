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

#include <format>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include <arrow/api.h>

namespace ActsExamples {

class ParquetReader::Impl {
 public:
  struct CollectionState {
    std::string name;
    std::unique_ptr<ActsPlugins::ArrowUtil::ParquetDatasetReader> reader;
    std::unique_ptr<WriteDataHandle<ActsPlugins::ArrowUtil::ArrowTable>> handle;
  };

  explicit Impl(ParquetReader::Config cfg, ParquetReader& parent)
      : m_cfg(std::move(cfg)) {
    if (m_cfg.collections.empty()) {
      throw std::invalid_argument("ParquetReader: no collections configured");
    }
    std::unordered_set<std::string> seenDirs;
    for (const auto& [name, rawPath] : m_cfg.collections) {
      if (name.empty()) {
        throw std::invalid_argument("ParquetReader: empty collection name");
      }
      if (rawPath.empty()) {
        throw std::invalid_argument(std::format(
            "ParquetReader: empty input directory for collection '{}'", name));
      }
      std::filesystem::path resolved =
          rawPath.is_absolute() ? rawPath : m_cfg.inputDir / rawPath;
      if (!seenDirs.insert(resolved.lexically_normal().string()).second) {
        throw std::invalid_argument(
            std::format("ParquetReader: duplicate input directory '{}'",
                        resolved.string()));
      }
      auto schemaIt = m_cfg.expectedSchemas.find(name);
      if (schemaIt == m_cfg.expectedSchemas.end() || !schemaIt->second) {
        throw std::invalid_argument(std::format(
            "ParquetReader: collection '{}' has no expected schema. Every "
            "configured collection must declare an expected schema; the "
            "scanner uses it as the dataset's target so missing columns "
            "become nulls and extras are dropped.",
            name));
      }
    }
    for (const auto& [name, _] : m_cfg.expectedSchemas) {
      if (!m_cfg.collections.contains(name)) {
        throw std::invalid_argument(std::format(
            "ParquetReader: expectedSchemas has entry for '{}' but no matching "
            "collection",
            name));
      }
    }

    std::int64_t referenceEvents = -1;
    std::string referenceName;

    for (const auto& [name, rawPath] : m_cfg.collections) {
      std::filesystem::path directory =
          rawPath.is_absolute() ? rawPath : m_cfg.inputDir / rawPath;
      if (!std::filesystem::exists(directory)) {
        throw std::invalid_argument(std::format(
            "ParquetReader: missing directory '{}'", directory.string()));
      }

      auto state = std::make_unique<CollectionState>();
      state->name = name;
      state->reader =
          std::make_unique<ActsPlugins::ArrowUtil::ParquetDatasetReader>(
              directory, m_cfg.expectedSchemas.at(name).schema());

      const auto events = state->reader->numEvents();
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

      state->handle =
          std::make_unique<WriteDataHandle<ActsPlugins::ArrowUtil::ArrowTable>>(
              &parent, name);
      state->handle->initialize(name);
      m_collectionStates.push_back(std::move(state));
    }

    m_eventsRange = {
        0, referenceEvents < 0 ? 0 : static_cast<std::size_t>(referenceEvents)};
  }

  ParquetReader::Config m_cfg;
  std::vector<std::unique_ptr<CollectionState>> m_collectionStates;
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

  for (const auto& state : m_impl->m_collectionStates) {
    auto table = state->reader->readEvent(
        static_cast<std::uint64_t>(context.eventNumber));
    if (table == nullptr || table->num_rows() == 0) {
      ACTS_ERROR("ParquetReader: no row matched event "
                 << context.eventNumber << " in collection '" << state->name
                 << "'");
      return ProcessCode::ABORT;
    }
    (*state->handle)(context,
                     ActsPlugins::ArrowUtil::ArrowTable{std::move(table)});
  }

  return ProcessCode::SUCCESS;
}

const ParquetReader::Config& ParquetReader::config() const {
  return m_impl->m_cfg;
}

}  // namespace ActsExamples
