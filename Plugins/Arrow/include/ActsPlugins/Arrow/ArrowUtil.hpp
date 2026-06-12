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
#include <vector>

// @TODO: Forward decl enough?
#include <arrow/api.h>

// Forward declarations of the Arrow C Data Interface POD structs used by
// the export helpers below. The actual definitions live in
// `arrow/c/abi.h`; using a forward decl here keeps that header out of
// every public consumer's translation unit.
struct ArrowArray;
struct ArrowSchema;

namespace ActsPlugins::ArrowUtil {

/// Opaque, pybind-friendly wrapper around @c std::shared_ptr<arrow::Schema>.
///
/// We can't use @c py::class_<arrow::Schema, ...> directly: arrow is built
/// with hidden visibility inside @c libActsPluginArrow, so its typeinfo
/// isn't visible to other .so files (including the python module). Wrapping
/// it in this @c ACTS_ARROW_EXPORT type gives pybind a class it can resolve
/// across the .so boundary while keeping every arrow symbol private.
///
/// Also serves as the value type in @c ParquetReader::Config::expectedSchemas
/// so the same handle can be both produced by Python (from
/// @c acts.examples.arrow.particleSchema() etc.) and consumed by C++.
class ACTS_ARROW_EXPORT ArrowSchemaHandle {
 public:
  ArrowSchemaHandle() = default;
  explicit ArrowSchemaHandle(std::shared_ptr<arrow::Schema> schema)
      : m_schema(std::move(schema)) {}

  const std::shared_ptr<arrow::Schema>& schema() const { return m_schema; }
  explicit operator bool() const { return m_schema != nullptr; }

  std::string toString() const;
  std::vector<std::string> fieldNames() const;
  int numFields() const;

  /// Export the wrapped @c arrow::Schema through the Arrow C Data
  /// Interface. The caller-supplied @p out struct is populated; the
  /// arrow-side resources transfer to @p out and are released by its
  /// release callback per the C-Data spec.
  void exportToC(::ArrowSchema* out) const;

 private:
  std::shared_ptr<arrow::Schema> m_schema;
};

/// Opaque, pybind-friendly wrapper around @c std::shared_ptr<arrow::Table>.
///
/// Same rationale as @c ArrowSchemaHandle: arrow's typeinfo is hidden
/// inside @c libActsPluginArrow, so this wrapper is the visibility-exported
/// type that pybind can register and that downstream code (including
/// @c WhiteBoardRegistry-mediated Python algorithms) can resolve across the
/// @c .so boundary.
///
/// Used as the canonical whiteboard storage type for arrow tables: every
/// @c WriteDataHandle / @c ReadDataHandle / @c ConsumeDataHandle for a
/// table holds an @c ArrowTable rather than a bare
/// @c std::shared_ptr<arrow::Table>. This gives Python algorithms a typed
/// entry point via the registry, with the underlying buffers exposed
/// through the Arrow C Data Interface.
class ACTS_ARROW_EXPORT ArrowTable {
 public:
  ArrowTable() = default;
  explicit ArrowTable(std::shared_ptr<arrow::Table> table)
      : m_table(std::move(table)) {}

  const std::shared_ptr<arrow::Table>& table() const { return m_table; }
  explicit operator bool() const { return m_table != nullptr; }

  std::int64_t numRows() const;
  int numColumns() const;
  ArrowSchemaHandle schema() const;
  std::string toString() const;

  /// Export the wrapped @c arrow::Table through the Arrow C Data
  /// Interface. Multi-chunk tables are combined into a single record
  /// batch first; for the canonical 1-row-per-event layout this is a
  /// no-op. The caller-supplied @p out_schema and @p out_array structs
  /// are populated and ownership of the arrow-side resources transfers
  /// via each struct's release callback.
  void exportToC(::ArrowSchema* out_schema, ::ArrowArray* out_array) const;

  /// Build an @c ArrowTable from C Data Interface structs produced by
  /// another arrow implementation (typically pyarrow's
  /// @c __arrow_c_array__). Consumes @p in_schema and @p in_array per
  /// the C-Data spec — after this call their @c release callbacks are
  /// nulled out and the callers must not free them again.
  ///
  /// Throws if either struct is null or if the import fails.
  static ArrowTable importFromC(::ArrowSchema* in_schema,
                                ::ArrowArray* in_array);

 private:
  std::shared_ptr<arrow::Table> m_table;
};

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

/// Schema for the per-event particle table emitted by
/// @c ArrowParticleOutputConverter. Nested layout: one row per event, each
/// field a @c list<T> whose single list element holds all particles of that
/// event.
ACTS_ARROW_EXPORT std::shared_ptr<arrow::Schema> particleSchema();

/// Schema for the per-event track table emitted by
/// @c ArrowTrackOutputConverter. The @c t column is nullable at the outer
/// level so writers running with @c writeTime=false can emit a single null
/// per event row instead of an opened list of N inner nulls.
ACTS_ARROW_EXPORT std::shared_ptr<arrow::Schema> trackSchema();

/// Schema for the per-event TRUTH simulated-hit table emitted by
/// @c ArrowSimHitOutputConverter: one entry per sim-hit (ALL sim-hits, in
/// container order — the entry's position is the sim-hit id referenced by the
/// measurement table's @c simhit_ids).
ACTS_ARROW_EXPORT std::shared_ptr<arrow::Schema> simHitSchema();

/// Schema for the per-event RECO measurement table emitted by
/// @c ArrowMeasurementOutputConverter: one entry per measurement (the entry's
/// position is the measurement id referenced by track @c hit_ids), carrying
/// the measured local parameters, cluster-shape features, and truth links
/// (@c particle_ids and @c simhit_ids) into the particle / sim-hit tables.
ACTS_ARROW_EXPORT std::shared_ptr<arrow::Schema> measurementSchema();

/// Thin RAII wrapper around @c parquet::arrow::FileWriter that opens lazily
/// on first write so the schema can be taken from the first event's table.
/// Page index is enabled so the reader can locate the matching row's pages
/// in each column without decoding the entire row group on every read.
class ACTS_ARROW_EXPORT ParquetFileWriter {
 public:
  explicit ParquetFileWriter(std::filesystem::path path);
  ~ParquetFileWriter() noexcept;

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
///
/// Optionally accepts a target schema: when supplied, the dataset is
/// built with that schema instead of the inspected union, so fragments
/// missing one of its columns get the column materialized as nulls
/// (added-column schema evolution) and fragments carrying extra columns
/// have them dropped. The supplied schema must NOT include the
/// @c event_id column — the reader prepends it internally. When
/// omitted, the dataset schema is inferred from the fragments.
class ACTS_ARROW_EXPORT ParquetDatasetReader {
 public:
  explicit ParquetDatasetReader(
      std::filesystem::path directory,
      const std::shared_ptr<arrow::Schema>& targetSchema = nullptr);
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
