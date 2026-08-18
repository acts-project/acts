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
#include "ActsPlugins/Arrow/ArrowUtil.hpp"

#include <sstream>
#include <stdexcept>

#include <arrow/api.h>

namespace ActsExamples {

class ArrowInputConverter::Impl {
 public:
  Impl(ArrowInputConverter& parent, const std::string& inputTable)
      : m_inputTable(&parent, "InputTable") {
    m_inputTable.initialize(inputTable);
  }

  ReadDataHandle<ActsPlugins::ArrowUtil::ArrowTable> m_inputTable;
};

ArrowInputConverter::ArrowInputConverter(
    const std::string& name, const std::string& inputTable,
    std::unique_ptr<const Acts::Logger> logger)
    : IAlgorithm(name, std::move(logger)),
      m_impl(std::make_unique<Impl>(*this, inputTable)) {}

ArrowInputConverter::~ArrowInputConverter() = default;

ProcessCode ArrowInputConverter::execute(const AlgorithmContext& ctx) const {
  const auto& handle = m_impl->m_inputTable(ctx);
  if (!handle) {
    ACTS_ERROR("ArrowInputConverter '" << name() << "' received null table");
    return ProcessCode::ABORT;
  }
  const auto& table = handle.table();

  // Subset check, not equality: dataset reads may surface a unified schema
  // that includes fields from newer fragments which older converters don't
  // care about. We require every field the converter declares to be present
  // and to have a matching type, but allow extra fields.
  if (auto expected = expectedSchema(); expected != nullptr) {
    const auto& actual = *table->schema();
    std::ostringstream missing;
    bool ok = true;
    for (const auto& field : expected->fields()) {
      auto actualField = actual.GetFieldByName(field->name());
      if (actualField == nullptr) {
        if (!ok) {
          missing << ", ";
        }
        missing << field->name() << " (missing)";
        ok = false;
      } else if (!actualField->type()->Equals(*field->type())) {
        if (!ok) {
          missing << ", ";
        }
        missing << field->name() << " (type " << actualField->type()->ToString()
                << " != " << field->type()->ToString() << ")";
        ok = false;
      }
    }
    if (!ok) {
      ACTS_ERROR("ArrowInputConverter '"
                 << name() << "' schema mismatch: " << missing.str() << "\n"
                 << "  expected fields: " << expected->ToString() << "\n"
                 << "  actual schema:   " << actual.ToString());
      return ProcessCode::ABORT;
    }
  }

  return convert(ctx, *table);
}

}  // namespace ActsExamples
