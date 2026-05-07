// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Arrow/ArrowUtil.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

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

class ParquetFileWriter::Impl {
 public:
  explicit Impl(std::filesystem::path path) : m_path(std::move(path)) {}

  void write(const arrow::Table& table) {
    if (!m_writer) {
      auto outfile = unwrap(arrow::io::FileOutputStream::Open(m_path.string()),
                            "open parquet");
      auto properties = parquet::WriterProperties::Builder()
                            .compression(parquet::Compression::ZSTD)
                            // ->enable_write_page_index()
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

class ParquetDatasetReader::Impl {
 public:
  explicit Impl(std::filesystem::path directory)
      : m_directory(std::move(directory)) {
    if (!std::filesystem::exists(m_directory) ||
        !std::filesystem::is_directory(m_directory)) {
      throw std::invalid_argument("ParquetDatasetReader: not a directory: " +
                                  m_directory.string());
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

    arrow::dataset::InspectOptions inspectOpts;
    inspectOpts.fragments =
        arrow::dataset::InspectOptions::kInspectAllFragments;
    arrow::dataset::FinishOptions finishOpts;
    finishOpts.inspect_options = inspectOpts;
    m_dataset = unwrap(factory->Finish(finishOpts), "finish dataset");

    auto fullSchema = m_dataset->schema();
    const int idx = fullSchema->GetFieldIndex(std::string{kEventIdColumn});
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
        arrow::compute::field_ref(std::string{kEventIdColumn}),
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

ParquetDatasetReader::ParquetDatasetReader(std::filesystem::path directory)
    : m_impl(std::make_unique<Impl>(std::move(directory))) {}

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
