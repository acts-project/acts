// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsPlugins/Arrow/Export.hpp"

#include <memory>
#include <string>
#include <vector>

#include <arrow/api.h>

namespace ActsExamples {

/// Test/debugging algorithm: pulls an @c arrow::Table off the whiteboard and
/// asserts shape properties against it. Aborts the sequencer on mismatch.
///
/// Useful for verifying that the @c ParquetReader's schema-projection
/// guarantees hold end-to-end (e.g. a column declared in the expected schema
/// but absent from the on-disk shards must arrive on the whiteboard with
/// nulls).
class ACTS_ARROW_EXPORT ArrowTableCheckAlg final : public IAlgorithm {
 public:
  struct Config {
    /// Whiteboard key of the @c std::shared_ptr<arrow::Table> to inspect.
    std::string inputTable;

    /// Columns that must be present on the table. A missing column aborts.
    std::vector<std::string> requiredColumns;

    /// Columns that must be present AND fully null (every chunked-array
    /// element is null). Useful for asserting that
    /// schema-projection-via-null-fill actually happened.
    std::vector<std::string> allNullColumns;
  };

  ArrowTableCheckAlg(const Config& cfg,
                     std::unique_ptr<const Acts::Logger> logger);

  ProcessCode execute(const AlgorithmContext& ctx) const override;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<std::shared_ptr<arrow::Table>> m_inputTable{this,
                                                             "InputTable"};
};

}  // namespace ActsExamples
