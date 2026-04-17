// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Arrow/ArrowUtil.hpp"

#include <stdexcept>

#include <arrow/compute/api.h>
#include <arrow/io/file.h>

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

}  // namespace

std::shared_ptr<arrow::Field> eventIdField() {
  return arrow::field(std::string{kEventIdColumn}, arrow::uint32(),
                      /*nullable=*/false);
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
  if (table->schema()->GetFieldIndex(std::string{kEventIdColumn}) != -1) {
    throw std::invalid_argument("withEventId: table already has event_id");
  }

  arrow::UInt32Builder builder;
  auto status = builder.Append(static_cast<std::uint32_t>(eventId));
  if (!status.ok()) {
    throwArrow("event_id append", status);
  }
  auto array = unwrap(builder.Finish(), "event_id finish");
  auto chunked = std::make_shared<arrow::ChunkedArray>(std::move(array));

  return unwrap(table->AddColumn(0, eventIdField(), std::move(chunked)),
                "add event_id column");
}

std::shared_ptr<arrow::Table> sliceByEventId(
    const std::shared_ptr<arrow::Table>& table, std::uint64_t eventId) {
  if (table == nullptr) {
    throw std::invalid_argument("sliceByEventId: null table");
  }
  const int idx = table->schema()->GetFieldIndex(std::string{kEventIdColumn});
  if (idx < 0) {
    throw std::invalid_argument("sliceByEventId: table lacks event_id column");
  }

  auto eventIdScalar = arrow::MakeScalar(static_cast<std::uint32_t>(eventId));
  auto mask = unwrap(arrow::compute::CallFunction(
                         "equal", {table->column(idx), eventIdScalar}),
                     "event_id equality");
  auto filtered =
      unwrap(arrow::compute::CallFunction("filter", {table, mask}),
             "event_id filter");

  auto filteredTable = filtered.table();
  return unwrap(filteredTable->RemoveColumn(idx), "drop event_id column");
}

class ParquetFileWriter::Impl {
 public:
  explicit Impl(std::filesystem::path path) : m_path(std::move(path)) {}

  void write(const arrow::Table& table) {
    if (!m_writer) {
      auto outfile = unwrap(
          arrow::io::FileOutputStream::Open(m_path.string()), "open parquet");
      auto properties = parquet::WriterProperties::Builder()
                            .compression(parquet::Compression::ZSTD)
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

ParquetFileWriter::~ParquetFileWriter() {
  if (m_impl) {
    m_impl->close();
  }
}

void ParquetFileWriter::write(const arrow::Table& table) {
  m_impl->write(table);
}

void ParquetFileWriter::close() {
  m_impl->close();
}

std::shared_ptr<arrow::Table> readTable(const std::filesystem::path& path) {
  auto infile = unwrap(
      arrow::io::ReadableFile::Open(path.string(), arrow::default_memory_pool()),
      "open parquet for read");

  auto reader = unwrap(
      parquet::arrow::OpenFile(infile, arrow::default_memory_pool()),
      "open parquet reader");

  std::shared_ptr<arrow::Table> table;
  auto status = reader->ReadTable(&table);
  if (!status.ok()) {
    throwArrow("read parquet table", status);
  }
  return table;
}

std::int64_t numRowsInFile(const std::filesystem::path& path) {
  auto infile = unwrap(
      arrow::io::ReadableFile::Open(path.string(), arrow::default_memory_pool()),
      "open parquet footer");

  auto reader = unwrap(
      parquet::arrow::OpenFile(infile, arrow::default_memory_pool()),
      "open parquet reader");
  return reader->parquet_reader()->metadata()->num_rows();
}

}  // namespace ActsPlugins::ArrowUtil
