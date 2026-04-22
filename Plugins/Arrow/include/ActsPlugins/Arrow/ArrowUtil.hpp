// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsPlugins/Arrow/Export.hpp"

#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>

#include <arrow/api.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>

namespace ActsPlugins::ArrowUtil {

/// Canonical column name used by the Parquet reader/writer to stamp and slice
/// rows per event. Every @c arrow::Table handed to @c ParquetWriter must carry
/// this column; every table returned from @c ParquetReader is already filtered
/// to a single event's rows.
inline constexpr std::string_view kEventIdColumn = "event_id";

/// Field for the canonical event id column. Shared across all collections so
/// the reader/writer do not need to know per-collection schemas. @c uint32
/// matches the convention of existing ACTS Parquet outputs.
ACTS_ARROW_EXPORT std::shared_ptr<arrow::Field> eventIdField();

/// Read just the number of rows from a Parquet file's footer. Metadata-only:
/// no data columns are decoded.
///
/// @param path The path to the Parquet file.
/// @return The total row count across all row groups.
ACTS_ARROW_EXPORT std::int64_t numRowsInFile(const std::filesystem::path& path);

/// Prepend an @c event_id column to a 1-row (nested-layout) table.
///
/// @param table The input table. Must have exactly one row and must not
///              already contain an @c event_id column.
/// @param eventId The event id value to stamp.
/// @return A new table with the event id column prepended.
ACTS_ARROW_EXPORT std::shared_ptr<arrow::Table> withEventId(
    const std::shared_ptr<arrow::Table>& table, std::uint64_t eventId);

/// Return the subset of rows where @c event_id == @p eventId.
///
/// @param table The input table. Must contain an @c event_id column.
/// @param eventId The event id to select.
/// @return A new table with only the matching rows. The @c event_id column is
///         stripped on the way out.
ACTS_ARROW_EXPORT std::shared_ptr<arrow::Table> sliceByEventId(
    const std::shared_ptr<arrow::Table>& table, std::uint64_t eventId);

/// Thin RAII wrapper around @c parquet::arrow::FileWriter that opens lazily on
/// first write so the schema can be taken from the first event's table.
class ACTS_ARROW_EXPORT ParquetFileWriter {
 public:
  explicit ParquetFileWriter(std::filesystem::path path);
  ~ParquetFileWriter();

  ParquetFileWriter(const ParquetFileWriter&) = delete;
  ParquetFileWriter& operator=(const ParquetFileWriter&) = delete;

  /// Append a table as a single row group. Opens the underlying file on first
  /// call using the table's schema.
  ///
  /// @param table The table to write. Its schema must match across calls.
  void write(const arrow::Table& table);

  /// Close the underlying file. Idempotent.
  void close();

 private:
  class Impl;
  std::unique_ptr<Impl> m_impl;
};

/// Read a Parquet file into a single @c arrow::Table in memory. Convenience
/// wrapper used by the example reader; large files should use a streaming API.
///
/// @param path The path to the Parquet file.
/// @return The full table contents.
ACTS_ARROW_EXPORT std::shared_ptr<arrow::Table> readTable(
    const std::filesystem::path& path);

}  // namespace ActsPlugins::ArrowUtil
