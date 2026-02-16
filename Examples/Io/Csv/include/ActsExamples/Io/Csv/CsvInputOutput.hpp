// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Concepts.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <limits>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include <boost/describe.hpp>
#include <boost/mp11.hpp>

namespace ActsExamples {

/// Write arbitrary data as comma-separated values into a text file.
class CsvWriter {
 public:
  static constexpr char Delimiter = ',';

  CsvWriter() = delete;
  CsvWriter(const CsvWriter&) = delete;
  CsvWriter(CsvWriter&&) noexcept = default;
  ~CsvWriter() = default;
  CsvWriter& operator=(const CsvWriter&) = delete;
  CsvWriter& operator=(CsvWriter&&) noexcept = default;

  /// Create a file at the given path. Overwrites existing data.
  ///
  /// \param columns    Column names, fixes the number of columns for the file
  /// \param path       Path to the output file
  /// \param precision  Output floating point precision
  CsvWriter(const std::vector<std::string>& columns, const std::string& path,
            int precision = std::numeric_limits<double>::max_digits10);

  /// Append arguments as a new row to the file.
  ///
  /// Each argument corresponds to one column. The writer ensures that the
  /// number of columns written match the number of columns that were specified
  /// during construction.
  ///
  /// \note `std::vector` arguments are automatically unpacked and each entry
  ///       is written as a separate column.
  template <typename Arg0, typename... Args>
  void append(Arg0&& arg0, Args&&... args);

 private:
  std::ofstream m_file;
  std::size_t m_num_columns;

  template <typename T>
  unsigned write(T&& x, std::ostream& os)
    requires(Acts::Concepts::arithmetic<std::decay_t<T>> ||
             std::convertible_to<T, std::string>);
  template <typename T, typename Allocator>
  static unsigned write(const std::vector<T, Allocator>& xs, std::ostream& os);
};

/// Read arbitrary data as comma-separated values from a text file.
class CsvReader {
 public:
  static constexpr char Delimiter = ',';

  CsvReader() = delete;
  CsvReader(const CsvReader&) = delete;
  CsvReader(CsvReader&&) noexcept = default;
  ~CsvReader() = default;
  CsvReader& operator=(const CsvReader&) = delete;
  CsvReader& operator=(CsvReader&&) noexcept = default;

  /// Open a file at the given path.
  ///
  /// \param path Path to the input file
  explicit CsvReader(const std::string& path);

  /// Read the next line from the file.
  ///
  /// \returns true   if the line was successfully read
  /// \returns false  if no more lines are available
  bool read(std::vector<std::string>& columns);

  /// Return the number of lines read so far.
  std::size_t num_lines() const { return m_num_lines; }

 private:
  std::ifstream m_file;
  std::string m_line;
  std::size_t m_num_lines = 0;
};

namespace detail {

template <typename T>
void parse(const std::string& str, T& value) {
  // TODO use something w/ lower overhead than stringstream e.g.
  // std::from_chars
  std::istringstream is(str);
  is >> value;
}

/// Return member names as a vector of strings for a Boost.Describe-annotated
/// struct.
template <typename T>
std::vector<std::string> member_names() {
  using members =
      boost::describe::describe_members<T, boost::describe::mod_public>;
  std::vector<std::string> names;
  boost::mp11::mp_for_each<members>(
      [&](auto D) { names.emplace_back(D.name); });
  return names;
}

/// Number of public members in a Boost.Describe-annotated struct.
template <typename T>
inline constexpr std::size_t member_count_v = boost::mp11::mp_size<
    boost::describe::describe_members<T, boost::describe::mod_public>>::value;

}  // namespace detail

/// Write Boost.Describe-annotated records as comma-separated values into a
/// text file.
template <typename T>
class BoostDescribeCsvWriter {
 public:
  BoostDescribeCsvWriter() = delete;
  BoostDescribeCsvWriter(const BoostDescribeCsvWriter&) = delete;
  BoostDescribeCsvWriter(BoostDescribeCsvWriter&&) noexcept = default;
  ~BoostDescribeCsvWriter() = default;
  BoostDescribeCsvWriter& operator=(const BoostDescribeCsvWriter&) = delete;
  BoostDescribeCsvWriter& operator=(BoostDescribeCsvWriter&&) noexcept =
      default;

  /// Create a file at the given path. Overwrites existing data.
  ///
  /// \param path       Path to the output file
  /// \param precision  Output floating point precision
  explicit BoostDescribeCsvWriter(
      const std::string& path,
      int precision = std::numeric_limits<double>::max_digits10)
      : m_writer(detail::member_names<T>(), path, precision) {}

  /// Append a record to the file.
  void append(const T& record) {
    using members =
        boost::describe::describe_members<T, boost::describe::mod_public>;
    using member_list = boost::mp11::mp_apply<boost::mp11::mp_list, members>;
    append_impl(record, member_list{});
  }

 private:
  CsvWriter m_writer;

  template <typename... Ds>
  void append_impl(const T& record, boost::mp11::mp_list<Ds...> /*unused*/) {
    m_writer.append(record.*Ds::pointer...);
  }
};

/// Read Boost.Describe-annotated records as comma-separated values from a
/// text file.
///
/// The reader is strict about its input format to avoid ambiguities. If
/// header verification is disabled, the first line will be skipped and each
/// line must contain exactly as many columns as there are struct members in
/// exactly the same order. If header verification is enabled, the first line
/// is interpreted as the header. Names in the header must match exactly to
/// the struct members but can be in arbitrary order. The file can contain
/// extra columns that are not struct members. Each following row must have
/// exactly the same number of columns as the header.
template <typename T>
class BoostDescribeCsvReader {
 public:
  BoostDescribeCsvReader() = delete;
  BoostDescribeCsvReader(const BoostDescribeCsvReader&) = delete;
  BoostDescribeCsvReader(BoostDescribeCsvReader&&) noexcept = default;
  ~BoostDescribeCsvReader() = default;
  BoostDescribeCsvReader& operator=(const BoostDescribeCsvReader&) = delete;
  BoostDescribeCsvReader& operator=(BoostDescribeCsvReader&&) noexcept =
      default;

  /// Open a file at the given path.
  ///
  /// \param path              Path to the input file
  /// \param optional_columns  Record columns that can be missing in the file
  /// \param verify_header     true to check header column names, false to skip
  ///
  /// The set of optional columns must match names in the record. When allowing
  /// optional columns, header verification must be set to true.
  explicit BoostDescribeCsvReader(
      const std::string& path,
      const std::vector<std::string>& optional_columns = {},
      bool verify_header = true);

  /// Read the next record from the file.
  ///
  /// Extra columns in the file will be ignored. Elements of the record that
  /// correspond to missing, optional columns will not be set and retain
  /// their value.
  ///
  /// \returns true   if a record was successfully read
  /// \returns false  if no more records are available
  bool read(T& record);

  /// Read the next record and any extra columns from the file.
  ///
  /// \returns true   if a record was successfully read
  /// \returns false  if no more records are available
  template <typename U>
  bool read(T& record, std::vector<U>& extra);

  /// Return the number of additional columns that are not part of the struct.
  std::size_t num_extra_columns() const { return m_extra_columns.size(); }
  /// Return the number of records read so far.
  std::size_t num_records() const { return m_reader.num_lines() - 1u; }

 private:
  static constexpr std::size_t NumMembers = detail::member_count_v<T>;

  CsvReader m_reader;
  std::vector<std::string> m_columns;
  // #columns is fixed to a reasonable value after reading the header
  std::size_t m_num_columns = SIZE_MAX;
  // map member index to column index in the file, SIZE_MAX for missing
  std::array<std::size_t, NumMembers> m_member_column_map;
  // column indices that do not map to a struct member
  std::vector<std::size_t> m_extra_columns;

  void use_default_columns();
  void parse_header(const std::vector<std::string>& optional_columns);
  void parse_record(T& record) const;
};

// implementation CsvWriter

inline CsvWriter::CsvWriter(const std::vector<std::string>& columns,
                             const std::string& path, int precision)
    : m_file(path,
             std::ios_base::binary | std::ios_base::out | std::ios_base::trunc),
      m_num_columns(columns.size()) {
  if (!m_file.is_open() || m_file.fail()) {
    throw std::runtime_error("Could not open file '" + path + "'");
  }
  m_file.precision(precision);
  if (m_num_columns == 0) {
    throw std::invalid_argument("No columns were specified");
  }
  // write column names as header row
  append(columns);
}

template <typename Arg0, typename... Args>
inline void CsvWriter::append(Arg0&& arg0, Args&&... args) {
  // we can only check how many columns were written after they have been
  // written. write to temporary first to prevent bad data on file.
  std::stringstream line;
  // ensure consistent formatting
  line.precision(m_file.precision());
  unsigned written_columns[] = {
      // write the first item without a delimiter and store columns written
      write(std::forward<Arg0>(arg0), line),
      // for all other items, write the delimiter followed by the item itself
      (line << Delimiter, write(std::forward<Args>(args), line))...,
  };
  line << '\n';
  // validate that the total number of written columns matches the specs.
  unsigned total_columns = 0;
  for (auto nc : written_columns) {
    total_columns += nc;
  }
  if (total_columns < m_num_columns) {
    throw std::invalid_argument("Not enough columns");
  }
  if (m_num_columns < total_columns) {
    throw std::invalid_argument("Too many columns");
  }
  // write the line to disk and check that it actually happened
  m_file << line.rdbuf();
  if (!m_file.good()) {
    throw std::runtime_error("Could not write data to file");
  }
}

template <typename T>
inline unsigned CsvWriter::write(T&& x, std::ostream& os)
  requires(Acts::Concepts::arithmetic<std::decay_t<T>> ||
           std::convertible_to<T, std::string>)
{
  os << x;
  return 1u;
}

template <typename T, typename Allocator>
inline unsigned CsvWriter::write(const std::vector<T, Allocator>& xs,
                                 std::ostream& os) {
  unsigned n = 0;
  for (const auto& x : xs) {
    if (0 < n) {
      os << Delimiter;
    }
    os << x;
    n += 1;
  }
  return n;
}

// implementation CsvReader

inline CsvReader::CsvReader(const std::string& path)
    : m_file(path, std::ios_base::binary | std::ios_base::in) {
  if (!m_file.is_open() || m_file.fail()) {
    throw std::runtime_error("Could not open file '" + path + "'");
  }
}

inline bool CsvReader::read(std::vector<std::string>& columns) {
  // read the next line and check for both end-of-file and errors
  std::getline(m_file, m_line);
  if (m_file.eof()) {
    return false;
  }
  if (m_file.fail()) {
    throw std::runtime_error(std::string("Could not read line ") +
                             std::to_string(m_num_lines));
  }
  m_num_lines += 1;

  // split the line into columns
  columns.clear();
  for (std::string::size_type pos = 0; pos < m_line.size();) {
    auto del = m_line.find_first_of(Delimiter, pos);
    if (del == std::string::npos) {
      // reached the end of the line; also determines the last column
      columns.emplace_back(m_line, pos);
      break;
    } else {
      columns.emplace_back(m_line, pos, del - pos);
      // start next column search after the delimiter
      pos = del + 1;
    }
  }
  return true;
}

// implementation BoostDescribeCsvReader

template <typename T>
inline BoostDescribeCsvReader<T>::BoostDescribeCsvReader(
    const std::string& path, const std::vector<std::string>& optional_columns,
    bool verify_header)
    : m_reader(path) {
  // optional columns only work if we verify the header
  if ((!optional_columns.empty()) && (!verify_header)) {
    throw std::runtime_error(
        "Optional columns can not be used without header verification");
  }
  // first line is always the header
  if (!m_reader.read(m_columns)) {
    throw std::runtime_error("Could not read header from '" + path + "'");
  }
  if (verify_header) {
    parse_header(optional_columns);
  } else {
    use_default_columns();
  }
}

template <typename T>
inline bool BoostDescribeCsvReader<T>::read(T& record) {
  if (!m_reader.read(m_columns)) {
    return false;
  }
  // check for consistent entries per-line
  if (m_columns.size() < m_num_columns) {
    throw std::runtime_error("Too few columns in line " +
                             std::to_string(m_reader.num_lines()));
  }
  if (m_num_columns < m_columns.size()) {
    throw std::runtime_error("Too many columns in line " +
                             std::to_string(m_reader.num_lines()));
  }
  parse_record(record);
  return true;
}

template <typename T>
template <typename U>
inline bool BoostDescribeCsvReader<T>::read(T& record, std::vector<U>& extra) {
  // parse columns belonging to the regular record
  if (!read(record)) {
    return false;
  }
  // parse extra columns
  extra.resize(m_extra_columns.size());
  for (std::size_t i = 0; i < m_extra_columns.size(); ++i) {
    detail::parse(m_columns[m_extra_columns[i]], extra[i]);
  }
  return true;
}

template <typename T>
inline void BoostDescribeCsvReader<T>::use_default_columns() {
  // assume row content is identical in content and order to the struct
  m_num_columns = NumMembers;
  for (std::size_t i = 0; i < m_member_column_map.size(); ++i) {
    m_member_column_map[i] = i;
  }
  // no extra columns by construction
  m_extra_columns.clear();
}

template <typename T>
inline void BoostDescribeCsvReader<T>::parse_header(
    const std::vector<std::string>& optional_columns) {
  const auto names = detail::member_names<T>();

  // the number of header columns fixes the expected number of data columns
  m_num_columns = m_columns.size();

  // check that all non-optional columns are available
  for (const auto& name : names) {
    // no need to check availability if the column is optional
    if (Acts::rangeContainsValue(optional_columns, name)) {
      continue;
    }
    // missing, non-optional column means we can not continue
    if (!Acts::rangeContainsValue(m_columns, name)) {
      throw std::runtime_error("Missing header column '" + name + "'");
    }
  }

  // ensure missing columns are correctly marked as such
  m_member_column_map.fill(SIZE_MAX);

  // determine column-member mapping and extra column indices
  m_extra_columns.clear();
  for (std::size_t i = 0; i < m_columns.size(); ++i) {
    // find the position of the column in the member names.
    auto it = std::ranges::find(names, m_columns[i]);
    if (it != names.end()) {
      // establish mapping between column and member position
      m_member_column_map[std::distance(names.begin(), it)] = i;
    } else {
      // record non-member columns
      m_extra_columns.push_back(i);
    }
  }
}

template <typename T>
inline void BoostDescribeCsvReader<T>::parse_record(T& record) const {
  using members =
      boost::describe::describe_members<T, boost::describe::mod_public>;
  std::size_t i = 0;
  boost::mp11::mp_for_each<members>([&](auto D) {
    if (m_member_column_map[i] != SIZE_MAX) {
      detail::parse(m_columns[m_member_column_map[i]], record.*D.pointer);
    }
    ++i;
  });
}

}  // namespace ActsExamples
