// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Arrow/ArrowUtil.hpp"

#include <algorithm>
#include <exception>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <vector>

#include <arrow/c/abi.h>
#include <arrow/c/bridge.h>
#include <arrow/compute/api.h>
#include <arrow/dataset/dataset.h>
#include <arrow/dataset/discovery.h>
#include <arrow/dataset/file_parquet.h>
#include <arrow/dataset/scanner.h>
#include <arrow/filesystem/localfs.h>
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>
#include <parquet/arrow/writer.h>

namespace ActsPlugins::ArrowUtil {

namespace {

[[noreturn]] void throwArrow(const std::string& what,
                             const arrow::Status& status) {
  throw std::runtime_error(what + ": " + status.ToString());
}

template <typename T>
T unwrap(arrow::Result<T> result, const std::string& what) {
  if (!result.ok()) {
    throwArrow(what, result.status());
  }
  return std::move(result).ValueOrDie();
}

// Arrow's compute kernels (e.g. `equal` used for the event_id filter
// pushdown below) are registered lazily. With the linker-isolated arrow
// island the registry is empty until we ask for it explicitly — without
// this call, scan-time filtering fails with "No function registered with
// name: equal".
void ensureComputeInitialized() {
  static std::once_flag flag;
  std::call_once(flag, [] {
    auto status = arrow::compute::Initialize();
    if (!status.ok()) {
      throwArrow("arrow compute init", status);
    }
  });
}

}  // namespace

std::string ArrowSchemaHandle::toString() const {
  return m_schema ? m_schema->ToString() : std::string{"<null schema>"};
}

std::vector<std::string> ArrowSchemaHandle::fieldNames() const {
  std::vector<std::string> names;
  if (!m_schema) {
    return names;
  }
  names.reserve(m_schema->num_fields());
  for (int i = 0; i < m_schema->num_fields(); ++i) {
    names.push_back(m_schema->field(i)->name());
  }
  return names;
}

int ArrowSchemaHandle::numFields() const {
  return m_schema ? m_schema->num_fields() : 0;
}

void ArrowSchemaHandle::exportToC(::ArrowSchema* out) const {
  if (out == nullptr) {
    throw std::invalid_argument("ArrowSchemaHandle::exportToC: null out");
  }
  if (!m_schema) {
    throw std::runtime_error(
        "ArrowSchemaHandle::exportToC: handle holds a null schema");
  }
  auto status = arrow::ExportSchema(*m_schema, out);
  if (!status.ok()) {
    throwArrow("ExportSchema", status);
  }
}

std::int64_t ArrowTable::numRows() const {
  return m_table ? m_table->num_rows() : 0;
}

int ArrowTable::numColumns() const {
  return m_table ? m_table->num_columns() : 0;
}

ArrowSchemaHandle ArrowTable::schema() const {
  return ArrowSchemaHandle{m_table ? m_table->schema() : nullptr};
}

std::string ArrowTable::toString() const {
  return m_table ? m_table->ToString() : std::string{"<null table>"};
}

void ArrowTable::exportToC(::ArrowSchema* out_schema,
                           ::ArrowArray* out_array) const {
  if (out_schema == nullptr || out_array == nullptr) {
    throw std::invalid_argument("ArrowTable::exportToC: null out");
  }
  if (!m_table) {
    throw std::runtime_error(
        "ArrowTable::exportToC: handle holds a null table");
  }
  // Combine to a single record batch — required by ExportRecordBatch and
  // matches what consumers receive through the `__arrow_c_array__`
  // protocol. For the canonical 1-row 1-chunk case this is essentially a
  // pointer copy.
  auto batch =
      unwrap(m_table->CombineChunksToBatch(), "table CombineChunksToBatch");
  auto status = arrow::ExportRecordBatch(*batch, out_array, out_schema);
  if (!status.ok()) {
    throwArrow("ExportRecordBatch", status);
  }
}

ArrowTable ArrowTable::importFromC(::ArrowSchema* in_schema,
                                   ::ArrowArray* in_array) {
  if (in_schema == nullptr || in_array == nullptr) {
    throw std::invalid_argument("ArrowTable::importFromC: null input");
  }
  // ImportRecordBatch consumes the C-Data structs (their release
  // callbacks are nulled out on success). Buffers are referenced via the
  // batch's internal release wiring, so the producer's memory stays
  // alive until our arrow::RecordBatch is destroyed.
  auto batch = unwrap(arrow::ImportRecordBatch(in_array, in_schema),
                      "ImportRecordBatch");
  auto table = unwrap(arrow::Table::FromRecordBatches({std::move(batch)}),
                      "Table::FromRecordBatches");
  return ArrowTable{std::move(table)};
}

std::shared_ptr<arrow::Field> eventIdField() {
  return arrow::field(std::string{kEventIdColumn}, arrow::uint32(),
                      /*nullable=*/false);
}

namespace {

// `list<T>` whose inner elements are nullable. Used for perigee parameters
// and track parameters so a per-element failure (no reference surface,
// failed local transform, failed propagation) emits a real null instead
// of a NaN sentinel.
std::shared_ptr<arrow::DataType> nullableFloatList() {
  return arrow::list(arrow::field("item", arrow::float32(), true));
}

}  // namespace

std::shared_ptr<arrow::Schema> particleSchema() {
  return arrow::schema({
      arrow::field("particle_id", arrow::list(arrow::uint64()), false),
      arrow::field("pdg_id", arrow::list(arrow::int64()), false),
      arrow::field("mass", arrow::list(arrow::float32()), false),
      arrow::field("energy", arrow::list(arrow::float32()), false),
      arrow::field("charge", arrow::list(arrow::float32()), false),
      arrow::field("vx", arrow::list(arrow::float32()), false),
      arrow::field("vy", arrow::list(arrow::float32()), false),
      arrow::field("vz", arrow::list(arrow::float32()), false),
      arrow::field("time", arrow::list(arrow::float32()), false),
      arrow::field("px", arrow::list(arrow::float32()), false),
      arrow::field("py", arrow::list(arrow::float32()), false),
      arrow::field("pz", arrow::list(arrow::float32()), false),
      arrow::field("perigee_d0", nullableFloatList(), false),
      arrow::field("perigee_z0", nullableFloatList(), false),
      arrow::field("vertex_primary", arrow::list(arrow::uint16()), false),
      arrow::field("parent_id", arrow::list(arrow::int64()), false),
      arrow::field("primary", arrow::list(arrow::boolean()), false),
  });
}

std::shared_ptr<arrow::Schema> trackSchema() {
  return arrow::schema({
      arrow::field("d0", nullableFloatList(), false),
      arrow::field("z0", nullableFloatList(), false),
      arrow::field("phi", nullableFloatList(), false),
      arrow::field("theta", nullableFloatList(), false),
      arrow::field("qop", nullableFloatList(), false),
      arrow::field("majority_particle_id", arrow::list(arrow::uint64()), false),
      arrow::field("hit_ids", arrow::list(arrow::list(arrow::uint32())), false),
      arrow::field("track_id", arrow::list(arrow::uint16()), false),
      arrow::field("t", nullableFloatList(), true),
  });
}

std::shared_ptr<arrow::Schema> simHitSchema() {
  return arrow::schema({
      arrow::field("x", arrow::list(arrow::float32()), false),
      arrow::field("y", arrow::list(arrow::float32()), false),
      arrow::field("z", arrow::list(arrow::float32()), false),
      arrow::field("true_x", arrow::list(arrow::float32()), false),
      arrow::field("true_y", arrow::list(arrow::float32()), false),
      arrow::field("true_z", arrow::list(arrow::float32()), false),
      arrow::field("time", arrow::list(arrow::float32()), false),
      arrow::field("particle_id", arrow::list(arrow::uint64()), false),
      arrow::field("detector", arrow::list(arrow::uint8()), false),
      arrow::field("volume_id", arrow::list(arrow::uint8()), false),
      arrow::field("layer_id", arrow::list(arrow::uint16()), false),
      arrow::field("surface_id", arrow::list(arrow::uint32()), false),
  });
}

std::shared_ptr<arrow::Schema> caloHitSchema() {
  return arrow::schema({
      arrow::field("detector", arrow::list(arrow::uint8()), false),
      arrow::field("total_energy", arrow::list(arrow::float32()), false),
      arrow::field("x", arrow::list(arrow::float32()), false),
      arrow::field("y", arrow::list(arrow::float32()), false),
      arrow::field("z", arrow::list(arrow::float32()), false),
      arrow::field("contrib_particle_ids",
                   arrow::list(arrow::list(arrow::uint64())), false),
      arrow::field("contrib_energies",
                   arrow::list(arrow::list(arrow::float32())), false),
      arrow::field("contrib_times", arrow::list(arrow::list(arrow::float32())),
                   false),
  });
}

std::shared_ptr<arrow::Table> withEventId(
    const std::shared_ptr<arrow::Table>& table, std::uint64_t eventId) {
  if (table == nullptr) {
    throw std::invalid_argument("withEventId: null table");
  }
  if (table->num_rows() != 1) {
    throw std::invalid_argument(
        "withEventId: expected a 1-row (nested-layout) table, got " +
        std::to_string(table->num_rows()) + " rows");
  }
  if (table->schema()->GetFieldIndex(kEventIdColumn) != -1) {
    throw std::invalid_argument("withEventId: table already has event_id");
  }

  arrow::UInt32Builder builder;
  if (auto status = builder.Append(static_cast<std::uint32_t>(eventId));
      !status.ok()) {
    throwArrow("event_id append", status);
  }
  auto array = unwrap(builder.Finish(), "event_id finish");
  auto chunked = std::make_shared<arrow::ChunkedArray>(std::move(array));

  return unwrap(table->AddColumn(0, eventIdField(), std::move(chunked)),
                "add event_id column");
}

class ParquetFileWriter::Impl {
 public:
  explicit Impl(std::filesystem::path path) : m_path(std::move(path)) {}

  void write(const arrow::Table& table) {
    if (!m_writer) {
      auto outfile = unwrap(arrow::io::FileOutputStream::Open(m_path.string()),
                            "open parquet");
      auto properties = parquet::WriterProperties::Builder()
                            .compression(parquet::Compression::ZSTD)
                            ->enable_write_page_index()
                            ->build();
      auto arrowProperties =
          parquet::ArrowWriterProperties::Builder().store_schema()->build();
      m_writer = unwrap(parquet::arrow::FileWriter::Open(
                            *table.schema(), arrow::default_memory_pool(),
                            outfile, properties, arrowProperties),
                        "open parquet writer");
    }
    auto status = m_writer->WriteTable(table, table.num_rows());
    if (!status.ok()) {
      throwArrow("parquet WriteTable", status);
    }
  }

  void close() {
    if (m_writer) {
      auto status = m_writer->Close();
      m_writer.reset();
      if (!status.ok()) {
        throwArrow("parquet close", status);
      }
    }
  }

 private:
  std::filesystem::path m_path;
  std::unique_ptr<parquet::arrow::FileWriter> m_writer;
};

ParquetFileWriter::ParquetFileWriter(std::filesystem::path path)
    : m_impl(std::make_unique<Impl>(std::move(path))) {}

ParquetFileWriter::~ParquetFileWriter() noexcept {
  if (m_impl) {
    try {
      m_impl->close();
    } catch (const std::exception& e) {
      std::cerr << "ParquetFileWriter::~ParquetFileWriter failed during close: "
                << e.what() << std::endl;
      std::terminate();
    } catch (...) {
      std::cerr << "ParquetFileWriter::~ParquetFileWriter failed during close "
                   "with an unknown exception"
                << std::endl;
      std::terminate();
    }
  }
}

void ParquetFileWriter::write(const arrow::Table& table) {
  m_impl->write(table);
}

void ParquetFileWriter::close() {
  m_impl->close();
}

class ParquetDatasetReader::Impl {
 public:
  Impl(std::filesystem::path directory,
       const std::shared_ptr<arrow::Schema>& targetSchema)
      : m_directory(std::move(directory)) {
    ensureComputeInitialized();

    if (!std::filesystem::exists(m_directory) ||
        !std::filesystem::is_directory(m_directory)) {
      throw std::invalid_argument("ParquetDatasetReader: not a directory: " +
                                  m_directory.string());
    }

    if (targetSchema != nullptr &&
        targetSchema->GetFieldIndex(kEventIdColumn) != -1) {
      throw std::invalid_argument(
          "ParquetDatasetReader: target schema must not contain event_id; "
          "the reader prepends it internally");
    }

    std::vector<std::string> files;
    for (const auto& e : std::filesystem::directory_iterator(m_directory)) {
      if (e.is_regular_file() && e.path().extension() == ".parquet") {
        files.push_back(e.path().string());
      }
    }
    if (files.empty()) {
      throw std::invalid_argument(
          "ParquetDatasetReader: no parquet files under " +
          m_directory.string());
    }
    std::ranges::sort(files);

    m_numEvents = 0;
    for (const auto& f : files) {
      auto infile =
          unwrap(arrow::io::ReadableFile::Open(f, arrow::default_memory_pool()),
                 "open parquet footer");
      auto reader =
          unwrap(parquet::arrow::OpenFile(infile, arrow::default_memory_pool()),
                 "open parquet reader");
      m_numEvents += reader->parquet_reader()->metadata()->num_rows();
    }

    auto fs = std::make_shared<arrow::fs::LocalFileSystem>();
    auto format = std::make_shared<arrow::dataset::ParquetFileFormat>();
    arrow::dataset::FileSystemFactoryOptions options;
    auto factory = unwrap(arrow::dataset::FileSystemDatasetFactory::Make(
                              std::move(fs), files, std::move(format), options),
                          "make dataset factory");

    arrow::dataset::FinishOptions finishOpts;
    if (targetSchema != nullptr) {
      // Prepend event_id; the dataset must carry it for the per-event
      // filter, but callers see it as an internal column and the
      // public schema below strips it again.
      auto withEventId =
          unwrap(targetSchema->AddField(0, eventIdField()), "prepend event_id");
      finishOpts.schema = std::move(withEventId);
    } else {
      arrow::dataset::InspectOptions inspectOpts;
      inspectOpts.fragments =
          arrow::dataset::InspectOptions::kInspectAllFragments;
      finishOpts.inspect_options = inspectOpts;
    }
    m_dataset = unwrap(factory->Finish(finishOpts), "finish dataset");

    auto fullSchema = m_dataset->schema();
    const int idx = fullSchema->GetFieldIndex(kEventIdColumn);
    if (idx < 0) {
      throw std::invalid_argument("ParquetDatasetReader: dataset under '" +
                                  m_directory.string() +
                                  "' lacks event_id column");
    }
    m_publicSchema =
        unwrap(fullSchema->RemoveField(idx), "strip event_id from schema");
  }

  std::int64_t numEvents() const { return m_numEvents; }
  std::shared_ptr<arrow::Schema> schema() const { return m_publicSchema; }

  std::shared_ptr<arrow::Table> readEvent(std::uint64_t eventId) const {
    auto builder = unwrap(m_dataset->NewScan(), "new scan");
    auto status = builder->Filter(arrow::compute::equal(
        arrow::compute::field_ref(kEventIdColumn.data()),
        arrow::compute::literal(static_cast<std::uint32_t>(eventId))));
    if (!status.ok()) {
      throwArrow("set scan filter", status);
    }
    status = builder->UseThreads(false);
    if (!status.ok()) {
      throwArrow("set use threads", status);
    }
    auto scanner = unwrap(builder->Finish(), "finish scanner");
    auto table = unwrap(scanner->ToTable(), "scan to table");

    const int idx = table->schema()->GetFieldIndex(std::string{kEventIdColumn});
    if (idx < 0) {
      throw std::runtime_error(
          "ParquetDatasetReader: scanned table lacks event_id column");
    }
    return unwrap(table->RemoveColumn(idx), "drop event_id column");
  }

 private:
  std::filesystem::path m_directory;
  std::int64_t m_numEvents = 0;
  std::shared_ptr<arrow::dataset::Dataset> m_dataset;
  std::shared_ptr<arrow::Schema> m_publicSchema;
};

ParquetDatasetReader::ParquetDatasetReader(
    std::filesystem::path directory,
    const std::shared_ptr<arrow::Schema>& targetSchema)
    : m_impl(std::make_unique<Impl>(std::move(directory), targetSchema)) {}

ParquetDatasetReader::~ParquetDatasetReader() = default;

std::int64_t ParquetDatasetReader::numEvents() const {
  return m_impl->numEvents();
}

std::shared_ptr<arrow::Schema> ParquetDatasetReader::schema() const {
  return m_impl->schema();
}

std::shared_ptr<arrow::Table> ParquetDatasetReader::readEvent(
    std::uint64_t eventId) const {
  return m_impl->readEvent(eventId);
}

}  // namespace ActsPlugins::ArrowUtil
