// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Parquet/ArrowInputConverter.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"

#include <stdexcept>

#include <arrow/api.h>

namespace ActsExamples {

class ArrowInputConverter::Impl {
 public:
  Impl(ArrowInputConverter& parent, const std::string& inputTable)
      : m_inputTable(&parent, "InputTable") {
    m_inputTable.initialize(inputTable);
  }

  ReadDataHandle<std::shared_ptr<arrow::Table>> m_inputTable;
};

ArrowInputConverter::ArrowInputConverter(
    const std::string& name, const std::string& inputTable,
    std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm(name, std::move(logger)),
      m_impl(std::make_unique<Impl>(*this, inputTable)) {}

ArrowInputConverter::~ArrowInputConverter() = default;

ProcessCode ArrowInputConverter::execute(const AlgorithmContext& ctx) const {
  const auto& table = m_impl->m_inputTable(ctx);
  if (table == nullptr) {
    ACTS_ERROR("ArrowInputConverter '" << name() << "' received null table");
    return ProcessCode::ABORT;
  }

  if (auto expected = expectedSchema(); expected != nullptr) {
    if (!table->schema()->Equals(*expected, /*check_metadata=*/false)) {
      ACTS_ERROR("ArrowInputConverter '"
                 << name() << "' schema mismatch.\n"
                 << "  expected: " << expected->ToString() << "\n"
                 << "  got:      " << table->schema()->ToString());
      return ProcessCode::ABORT;
    }
  }

  return convert(ctx, *table);
}

}  // namespace ActsExamples
