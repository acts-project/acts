// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Arrow/ArrowTableCheckAlg.hpp"

#include "Acts/Utilities/Logger.hpp"

#include <stdexcept>

namespace ActsExamples {

ArrowTableCheckAlg::ArrowTableCheckAlg(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm("ArrowTableCheckAlg", std::move(logger)), m_cfg(cfg) {
  if (m_cfg.inputTable.empty()) {
    throw std::invalid_argument(
        "ArrowTableCheckAlg: inputTable must not be empty");
  }
  m_inputTable.initialize(m_cfg.inputTable);
}

ProcessCode ArrowTableCheckAlg::execute(const AlgorithmContext& ctx) const {
  const auto& table = m_inputTable(ctx);
  if (table == nullptr) {
    ACTS_ERROR("null table on whiteboard for '" << m_cfg.inputTable << "'");
    return ProcessCode::ABORT;
  }

  for (const auto& name : m_cfg.requiredColumns) {
    if (table->schema()->GetFieldIndex(name) < 0) {
      ACTS_ERROR("event " << ctx.eventNumber << ": required column '" << name
                          << "' missing from table '" << m_cfg.inputTable
                          << "'. Schema: " << table->schema()->ToString());
      return ProcessCode::ABORT;
    }
  }

  for (const auto& name : m_cfg.allNullColumns) {
    auto column = table->GetColumnByName(name);
    if (column == nullptr) {
      ACTS_ERROR("event " << ctx.eventNumber << ": all-null column '" << name
                          << "' missing from table '" << m_cfg.inputTable
                          << "'. Schema: " << table->schema()->ToString());
      return ProcessCode::ABORT;
    }
    if (column->null_count() != column->length()) {
      ACTS_ERROR("event " << ctx.eventNumber << ": column '" << name
                          << "' on table '" << m_cfg.inputTable
                          << "' expected to be all-null but had "
                          << (column->length() - column->null_count())
                          << " non-null of " << column->length() << " values");
      return ProcessCode::ABORT;
    }
  }

  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
