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

// @TODO: Forward decl enough?
#include <arrow/api.h>

namespace ActsPlugins::ArrowUtil {

/// Canonical column name used by the Parquet reader/writer to stamp and
/// filter rows per event. Every @c arrow::Table handed to @c ParquetWriter
/// must carry this column; every table returned from @c ParquetReader has
/// already been filtered to a single event's rows and has the column
/// stripped on the way out.
inline constexpr std::string_view kEventIdColumn = "event_id";

/// Field for the canonical event id column. Shared across all collections so
/// the reader/writer do not need to know per-collection schemas. @c uint32
/// matches the convention of existing ACTS Parquet outputs.
ACTS_ARROW_EXPORT std::shared_ptr<arrow::Field> eventIdField();

/// Prepend an @c event_id column to a 1-row (nested-layout) table.
///
/// @param table The input table. Must have exactly one row and must not
///              already contain an @c event_id column.
/// @param eventId The event id value to stamp.
/// @return A new table with the event id column prepended.
ACTS_ARROW_EXPORT std::shared_ptr<arrow::Table> withEventId(
    const std::shared_ptr<arrow::Table>& table, std::uint64_t eventId);

/// Thin RAII wrapper around @c parquet::arrow::FileWriter that opens lazily
/// on first write so the schema can be taken from the first event's table.
/// Page index is enabled so the reader can locate the matching row's pages
/// in each column without decoding the entire row group on every read.
class ACTS_ARROW_EXPORT ParquetFileWriter {
 public:
  explicit ParquetFileWriter(std::filesystem::path path);
  ~ParquetFileWriter();

  ParquetFileWriter(const ParquetFileWriter&) = delete;
  ParquetFileWriter& operator=(const ParquetFileWriter&) = delete;

  /// Append a table as a single row group. Opens the underlying file on
  /// first call using the table's schema.
  ///
  /// @param table The table to write. Its schema must match across calls.
  void write(const arrow::Table& table);

  /// Close the underlying file. Idempotent.
  void close();

 private:
  class Impl;
  std::unique_ptr<Impl> m_impl;
};

/// Reader for an @c arrow::dataset::FileSystemDataset of Parquet shard
/// files. Discovers all @c *.parquet files under @p directory at
/// construction, validates that each carries the @c event_id column, and
/// answers per-event lookups via filter pushdown into the dataset scanner.
///
/// File-level pruning via Parquet footer min/max statistics reduces a
/// per-event lookup to opening exactly one fragment when shards own
/// disjoint event-id ranges (the layout produced by @c ParquetWriter).
class ACTS_ARROW_EXPORT ParquetDatasetReader {
 public:
  explicit ParquetDatasetReader(std::filesystem::path directory);
  ~ParquetDatasetReader();

  ParquetDatasetReader(const ParquetDatasetReader&) = delete;
  ParquetDatasetReader& operator=(const ParquetDatasetReader&) = delete;

  /// Sum of row counts across all parquet fragments in the directory.
  /// Footer-only reads, performed at construction time and cached.
  std::int64_t numEvents() const;

  /// Unified schema inspected across fragments. The @c event_id column is
  /// stripped (it's an internal routing column).
  std::shared_ptr<arrow::Schema> schema() const;

  /// Filtered scan returning the (1-row) table where @c event_id matches.
  /// The @c event_id column is stripped on the way out.
  std::shared_ptr<arrow::Table> readEvent(std::uint64_t eventId) const;

 private:
  class Impl;
  std::unique_ptr<Impl> m_impl;
};

}  // namespace ActsPlugins::ArrowUtil
