// SPDX-License-Identifier: MIT
// Copyright 2015,2018-2019 Moritz Kiehn
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

/// \file
/// \brief   Write numpy-compatible .npy binary files
/// \author  Moritz Kiehn <msmk@cern.ch>
/// \date    2019-09-08, Split numpy i/o from the namedtuple library

#pragma once

#include <array>
#include <cstdint>
#include <fstream>
#include <string>
#include <tuple>
#include <utility>

namespace dfe {

/// Write records into a binary NumPy-compatible `.npy` file.
///
/// See
/// https://docs.scipy.org/doc/numpy/reference/generated/numpy.lib.format.html
/// for an explanation of the file format.
template<typename NamedTuple>
class NamedTupleNumpyWriter {
public:
  NamedTupleNumpyWriter() = delete;
  NamedTupleNumpyWriter(const NamedTupleNumpyWriter&) = delete;
  NamedTupleNumpyWriter(NamedTupleNumpyWriter&&) = default;
  ~NamedTupleNumpyWriter();
  NamedTupleNumpyWriter& operator=(const NamedTupleNumpyWriter&) = delete;
  NamedTupleNumpyWriter& operator=(NamedTupleNumpyWriter&&) = default;

  /// Create a npy file at the given path. Overwrites existing data.
  NamedTupleNumpyWriter(const std::string& path);

  /// Append a record to the end of the file.
  void append(const NamedTuple& record);

private:
  // the equivalent std::tuple-like type
  using Tuple = typename NamedTuple::Tuple;

  std::ofstream m_file;
  std::size_t m_fixed_header_length;
  std::size_t m_num_tuples;

  void write_header(std::size_t num_tuples);
  template<std::size_t... I>
  void write_record(const NamedTuple& record, std::index_sequence<I...>);
  template<typename T>
  void write_bytes(const T* ptr);
};

// implementation helpers
namespace io_npy_impl {

template<typename T>
constexpr std::enable_if_t<false, T> kNumpyDtypeCode;
template<>
constexpr const char* kNumpyDtypeCode<uint8_t> = "u1";
template<>
constexpr const char* kNumpyDtypeCode<uint16_t> = "u2";
template<>
constexpr const char* kNumpyDtypeCode<uint32_t> = "u4";
template<>
constexpr const char* kNumpyDtypeCode<uint64_t> = "u8";
template<>
constexpr const char* kNumpyDtypeCode<int8_t> = "i1";
template<>
constexpr const char* kNumpyDtypeCode<int16_t> = "i2";
template<>
constexpr const char* kNumpyDtypeCode<int32_t> = "i4";
template<>
constexpr const char* kNumpyDtypeCode<int64_t> = "i8";
template<>
constexpr const char* kNumpyDtypeCode<float> = "f4";
template<>
constexpr const char* kNumpyDtypeCode<double> = "f8";
template<>
constexpr const char* kNumpyDtypeCode<bool> = "?";

template<typename... Types>
constexpr std::array<const char*, sizeof...(Types)>
dtypes_codes(const std::tuple<Types...>&) {
  return {kNumpyDtypeCode<typename std::decay<Types>::type>...};
}

// Determines endianness and return the corresponding dtype code modifier.
//
// Derived from:
// https://stackoverflow.com/questions/1001307/detecting-endianness-programmatically-in-a-c-program
inline char
dtype_endianness_modifier() {
  union {
    uint32_t i;
    char c[4];
  } x = {0x0A0B0C0D};
  bool is_little_endian =
    (x.c[0] == 0xD) and (x.c[1] == 0xC) and (x.c[2] == 0xB) and (x.c[3] == 0xA);
  // TODO this assumes that only little and big endian exists and checks only
  // for little. maybe verify that it always is one or the other?
  return is_little_endian ? '<' : '>';
}

template<typename NamedTuple>
inline std::string
dtypes_description(const NamedTuple& nt) {
  std::string descr;
  std::size_t n = std::tuple_size<typename NamedTuple::Tuple>::value;
  auto names = nt.names();
  auto codes = dtypes_codes(nt.tuple());
  auto endianness_modifier = dtype_endianness_modifier();
  descr += '[';
  for (decltype(n) i = 0; i < n; ++i) {
    descr += "('";
    descr += names[i];
    descr += "', '";
    descr += endianness_modifier;
    descr += codes[i];
    descr += "')";
    if ((i + 1) < n) {
      descr += ", ";
    }
  }
  descr += ']';
  return descr;
}

} // namespace io_npy_impl

// implementation

template<typename NamedTuple>
inline NamedTupleNumpyWriter<NamedTuple>::NamedTupleNumpyWriter(
  const std::string& path)
  : m_fixed_header_length(0), m_num_tuples(0) {
  // make our life easier. always throw on error
  m_file.exceptions(std::ofstream::badbit | std::ofstream::failbit);
  m_file.open(
    path, std::ios_base::binary | std::ios_base::out | std::ios_base::trunc);
  // write a header that uses the maximum amount of space, i.e. biggest
  // possible number of ntuples, so that we have enough space when we
  // overwrite it w/ the actual number of tuples at closing time.
  write_header(SIZE_MAX);
  write_header(0);
}

template<typename NamedTuple>
inline NamedTupleNumpyWriter<NamedTuple>::~NamedTupleNumpyWriter() {
  if (!m_file.is_open()) {
    return;
  }
  write_header(m_num_tuples);
  m_file.close();
}

template<typename NamedTuple>
inline void
NamedTupleNumpyWriter<NamedTuple>::append(const NamedTuple& record) {
  write_record(
    record, std::make_index_sequence<std::tuple_size<Tuple>::value>{});
  m_num_tuples += 1;
}

template<typename NamedTuple>
inline void
NamedTupleNumpyWriter<NamedTuple>::write_header(std::size_t num_tuples) {
  std::string header;
  // magic
  header += "\x93NUMPY";
  // fixed version number (major, minor), 1byte unsigned each
  header += static_cast<char>(0x1);
  header += static_cast<char>(0x0);
  // placeholder value for the header length, 2byte little endian unsigned
  header += static_cast<char>(0xAF);
  header += static_cast<char>(0xFE);
  // python dict w/ data type and size information
  header += "{'descr': ";
  header += io_npy_impl::dtypes_description(NamedTuple());
  header += ", 'fortran_order': False";
  header += ", 'shape': (";
  header += std::to_string(num_tuples);
  header += ",)}";
  // padd w/ spaces for 16 byte alignment of the whole header
  while (((header.size() + 1) % 16) != 0) {
    header += ' ';
  }
  // the initial header fixes the available header size. updated headers
  // must always occupy the same space and might require additional
  // padding spaces
  if (m_fixed_header_length == 0) {
    m_fixed_header_length = header.size();
  } else {
    while (header.size() < m_fixed_header_length) {
      header += ' ';
    }
  }
  header += '\n';
  // replace the header length place holder
  std::size_t header_length = header.size() - 10;
  header[8] = static_cast<char>(header_length >> 0);
  header[9] = static_cast<char>(header_length >> 8);
  m_file.seekp(0);
  m_file.write(header.data(), header.size());
}

template<typename NamedTuple>
template<std::size_t... I>
inline void
NamedTupleNumpyWriter<NamedTuple>::write_record(
  const NamedTuple& record, std::index_sequence<I...>) {
  // see namedtuple_impl::print_tuple for explanation
  using std::get;
  using Vacuum = int[];
  (void)Vacuum{(write_bytes(&get<I>(record)), 0)...};
}

template<typename NamedTuple>
template<typename T>
inline void
NamedTupleNumpyWriter<NamedTuple>::write_bytes(const T* ptr) {
  m_file.write(reinterpret_cast<const char*>(ptr), sizeof(T));
}

} // namespace dfe
