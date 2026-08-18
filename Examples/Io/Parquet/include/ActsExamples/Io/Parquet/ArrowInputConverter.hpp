// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsPlugins/Arrow/Export.hpp"

#include <memory>
#include <string>

#include <arrow/api.h>

namespace ActsExamples {

/// Abstract base class for Arrow input converters.
///
/// Centralizes the retrieval of the input @c arrow::Table from the whiteboard
/// (where @c ParquetReader parked it, sliced to this event's rows) and
/// optionally validates its schema against subclass expectations before
/// delegating to @c convert(). Subclasses do not implement @c execute().
///
/// Unlike podio collections, @c arrow::Table has no intrinsic type identity —
/// "particles" and "tracks" have the same C++ type. Subclasses may override
/// @c expectedSchema() to opt into schema validation.
class ACTS_ARROW_EXPORT ArrowInputConverter : public IAlgorithm {
 public:
  /// Constructor for the Arrow input converter.
  ///
  /// @param name The name of the algorithm.
  /// @param inputTable The whiteboard key to read the input table from.
  /// @param logger The logger instance.
  ArrowInputConverter(const std::string& name, const std::string& inputTable,
                      std::unique_ptr<const Acts::Logger> logger = nullptr);

  ~ArrowInputConverter() override;

  /// Execute the algorithm. Subclasses do not implement this method; they
  /// implement @c convert() instead.
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Convert the input Arrow table to the internal EDM format.
  ///
  /// @param ctx The algorithm context.
  /// @param table The input Arrow table, already sliced to this event.
  /// @return The process code.
  virtual ProcessCode convert(const AlgorithmContext& ctx,
                              const arrow::Table& table) const = 0;

  /// Optional: declare the expected schema so the base class can validate the
  /// input table before calling @c convert(). Return @c nullptr to skip
  /// validation (the default).
  ///
  /// @return The expected schema, or @c nullptr for no validation.
  virtual std::shared_ptr<arrow::Schema> expectedSchema() const {
    return nullptr;
  }

 private:
  class Impl;

  std::unique_ptr<Impl> m_impl;
};

}  // namespace ActsExamples
