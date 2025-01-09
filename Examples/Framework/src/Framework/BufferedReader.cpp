// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Framework/BufferedReader.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <random>
#include <utility>

namespace ActsExamples {

BufferedReader::BufferedReader(const Config &config, Acts::Logging::Level level)
    : m_cfg(config), m_logger(Acts::getDefaultLogger(name(), level)) {
  if (!m_cfg.downstreamReader) {
    throw std::invalid_argument("No downstream reader provided!");
  }

  // Register write and read handles of the downstream reader
  for (auto rh : m_cfg.downstreamReader->readHandles()) {
    registerReadHandle(*rh);
  }

  for (auto wh : m_cfg.downstreamReader->writeHandles()) {
    registerWriteHandle(*wh);
  }

  // Read the events
  auto [ebegin, eend] = m_cfg.downstreamReader->availableEvents();
  if (eend - ebegin < m_cfg.bufferSize) {
    throw std::runtime_error("Reader does not provide enough events");
  }

  ACTS_INFO("Start reading events into buffer...");

  m_buffer.reserve(eend - ebegin);
  for (auto i = ebegin; i < ebegin + m_cfg.bufferSize; ++i) {
    auto board = std::make_unique<ActsExamples::WhiteBoard>(m_logger->clone());
    ActsExamples::AlgorithmContext ctx(0, i, *board);

    ACTS_DEBUG("Read event " << i << " into buffer");
    m_cfg.downstreamReader->read(ctx);
    m_buffer.emplace_back(std::move(board));
  }

  ACTS_INFO("Filled " << m_buffer.size() << " events into the buffer");
}

ProcessCode BufferedReader::read(const AlgorithmContext &ctx) {
  // Set up a random event selection that is consistent if multiple
  // BufferedReader are used within a workflow The linear congruential engine is
  // chosen since it is cheap to instantiate. For each eventNumber, it is put in
  // a reproducible state.
  std::minstd_rand rng(m_cfg.selectionSeed);
  rng.discard(ctx.eventNumber);

  /// Sample from the buffer and transfer the content
  std::uniform_int_distribution<std::size_t> dist(0, m_cfg.bufferSize - 1);

  const auto entry = dist(rng);
  m_buffer.at(entry)->shareDataWith(ctx.eventStore);

  ACTS_DEBUG("Use buffer entry " << entry << " for event " << ctx.eventNumber);

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
