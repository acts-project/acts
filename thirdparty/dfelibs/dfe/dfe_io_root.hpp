// SPDX-License-Identifier: MIT
// Copyright 2019 Moritz Kiehn
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
/// \brief   Wrappers to simplify reading/writing ROOT TTree-based files
/// \author  Moritz Kiehn <msmk@cern.ch>
/// \date    2019-04-00, Initial version

#pragma once

#include <array>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>

#include <TFile.h>
#include <TTree.h>

namespace dfe {

/// Write records into a ROOT TTree.
template<typename NamedTuple>
class NamedTupleRootWriter {
public:
  NamedTupleRootWriter() = delete;
  NamedTupleRootWriter(const NamedTupleRootWriter&) = delete;
  NamedTupleRootWriter(NamedTupleRootWriter&&) = delete;
  NamedTupleRootWriter& operator=(const NamedTupleRootWriter&) = delete;
  NamedTupleRootWriter& operator=(NamedTupleRootWriter&&) = delete;

  /// Create a file at the given path. Overwrites existing data.
  ///
  /// \param path       Path to the output file
  /// \param tree_name  Name of the output tree within the file
  NamedTupleRootWriter(const std::string& path, const std::string& tree_name);
  /// Create a tree in a ROOT directory. Overwrites existing data.
  ///
  /// \param dir        Output directory for the tree
  /// \param tree_name  Name of the output tree relative to the directory
  ///
  /// When the writer is created with an existing ROOT directory, the user
  /// is responsible for ensuring the underlying file is closed.
  NamedTupleRootWriter(TDirectory* dir, const std::string& tree_name);
  /// Write the tree and close the owned file.
  ~NamedTupleRootWriter();

  /// Append a record to the file.
  void append(const NamedTuple& record);

private:
  // the equivalent std::tuple-like type
  using Tuple = typename NamedTuple::Tuple;

  TFile* m_file;
  TTree* m_tree;
  Tuple m_data;

  template<std::size_t... I>
  void setup_branches(std::index_sequence<I...>);
};

/// Read records from a ROOT TTree.
template<typename NamedTuple>
class NamedTupleRootReader {
public:
  NamedTupleRootReader() = delete;
  NamedTupleRootReader(const NamedTupleRootReader&) = delete;
  NamedTupleRootReader(NamedTupleRootReader&&) = delete;
  NamedTupleRootReader& operator=(const NamedTupleRootReader&) = delete;
  NamedTupleRootReader& operator=(NamedTupleRootReader&&) = delete;

  /// Open a file at the given path.
  ///
  /// \param path       Path to the input file
  /// \param tree_name  Name of the input tree within the file
  NamedTupleRootReader(const std::string& path, const std::string& tree_name);
  /// Open a tree from a ROOT directory.
  ///
  /// \param dir        Input directory for the tree
  /// \param tree_name  Name of the input tree relative to the directory
  ///
  /// When the reader is created with an existing ROOT directory, the user
  /// is responsible for ensuring the underlying file is closed.
  NamedTupleRootReader(TDirectory* dir, const std::string& tree_name);
  /// Write the tree and close the owned file.
  ~NamedTupleRootReader();

  /// Read the next record from the file.
  ///
  /// \returns true   if a record was successfully read
  /// \returns false  if no more records are available
  bool read(NamedTuple& record);

private:
  // the equivalent std::tuple-like type
  using Tuple = typename NamedTuple::Tuple;

  TFile* m_file;
  TTree* m_tree;
  int64_t m_next;
  Tuple m_data;

  template<std::size_t... I>
  void setup_branches(std::index_sequence<I...>);
};

// implementation writer

template<typename NamedTuple>
inline NamedTupleRootWriter<NamedTuple>::NamedTupleRootWriter(
  const std::string& path, const std::string& tree_name)
  : m_file(new TFile(path.c_str(), "RECREATE"))
  , m_tree(new TTree(tree_name.c_str(), "", 99, m_file)) {
  if (not m_file) {
    throw std::runtime_error("Could not create file");
  }
  if (not m_file->IsOpen()) {
    throw std::runtime_error("Could not open file");
  }
  if (not m_tree) {
    throw std::runtime_error("Could not create tree");
  }
  setup_branches(std::make_index_sequence<std::tuple_size<Tuple>::value>());
}

template<typename NamedTuple>
inline NamedTupleRootWriter<NamedTuple>::NamedTupleRootWriter(
  TDirectory* dir, const std::string& tree_name)
  : m_file(nullptr) // no file since it is not owned by the writer
  , m_tree(new TTree(tree_name.c_str(), "", 99, dir)) {
  if (not dir) {
    throw std::runtime_error("Invalid output directory given");
  }
  if (not m_tree) {
    throw std::runtime_error("Could not create tree");
  }
  setup_branches(std::make_index_sequence<std::tuple_size<Tuple>::value>());
}

namespace namedtuple_root_impl {

// see: https://cpppatterns.com/patterns/class-template-sfinae.html
template<typename T, typename Enable = void>
struct TypeCode;
// non-integer types
template<char Code>
struct TypeCodePlainImpl {
  static constexpr char value = Code;
};
template<>
struct TypeCode<bool> : TypeCodePlainImpl<'O'> {};
template<>
struct TypeCode<float> : TypeCodePlainImpl<'F'> {};
template<>
struct TypeCode<double> : TypeCodePlainImpl<'D'> {};
// integer types
// you might think that you could just define this for all the stdint types;
// but no, this breaks because ROOT [U]Long64_t might not be the same type as
// [u]int64_t depending on the machine, os, moon phase, ... .
// Why you ask? Because the universe hates you.
template<typename T, char Unsigned, char Signed>
struct TypeCodeIntImpl {
  static constexpr char value = std::is_unsigned<T>::value ? Unsigned : Signed;
};
template<typename T, std::size_t S>
constexpr bool is_integer_with_size_v = std::is_integral<T>::value
                                        and (sizeof(T) == S);
template<typename T>
struct TypeCode<T, typename std::enable_if_t<is_integer_with_size_v<T, 1>>>
  : TypeCodeIntImpl<T, 'b', 'B'> {};
template<typename T>
struct TypeCode<T, typename std::enable_if_t<is_integer_with_size_v<T, 2>>>
  : TypeCodeIntImpl<T, 's', 'S'> {};
template<typename T>
struct TypeCode<T, typename std::enable_if_t<is_integer_with_size_v<T, 4>>>
  : TypeCodeIntImpl<T, 'i', 'I'> {};
template<typename T>
struct TypeCode<T, typename std::enable_if_t<is_integer_with_size_v<T, 8>>>
  : TypeCodeIntImpl<T, 'l', 'L'> {};

} // namespace namedtuple_root_impl

template<typename NamedTuple>
template<std::size_t... I>
inline void
NamedTupleRootWriter<NamedTuple>::setup_branches(std::index_sequence<I...>) {
  static_assert(
    sizeof...(I) == std::tuple_size<Tuple>::value, "Something is very wrong");

  // construct leaf names w/ type info
  std::array<std::string, sizeof...(I)> names = NamedTuple::names();
  std::array<std::string, sizeof...(I)> leafs = {
    (names[I] + '/'
     + namedtuple_root_impl::TypeCode<
       std::tuple_element_t<I, Tuple>>::value)...};
  // construct branches
  // NOTE 2019-05-13 msmk:
  // the documentation suggests that ROOT can figure out the branch types on
  // its own, but doing so seems to break for {u}int64_t. do it manually for
  // now.
  (void)std::array<TBranch*, sizeof...(I)>{m_tree->Branch(
    names[I].c_str(), &std::get<I>(m_data), leafs[I].c_str())...};
}

template<typename NamedTuple>
inline NamedTupleRootWriter<NamedTuple>::~NamedTupleRootWriter() {
  // alway overwrite old data
  if (m_tree) {
    m_tree->Write(nullptr, TObject::kOverwrite);
  }
  // writer owns the file
  if (m_file) {
    m_file->Close();
    delete m_file;
  }
}

template<typename NamedTuple>
inline void
NamedTupleRootWriter<NamedTuple>::append(const NamedTuple& record) {
  m_data = record;
  if (m_tree->Fill() == -1) {
    throw std::runtime_error("Could not fill an entry");
  }
}

// implementation reader

template<typename NamedTuple>
inline NamedTupleRootReader<NamedTuple>::NamedTupleRootReader(
  const std::string& path, const std::string& tree_name)
  : m_file(new TFile(path.c_str(), "READ")), m_tree(nullptr), m_next(0) {
  if (not m_file) {
    throw std::runtime_error("Could not open file");
  }
  if (not m_file->IsOpen()) {
    throw std::runtime_error("Could not open file");
  }
  m_tree = static_cast<TTree*>(m_file->Get(tree_name.c_str()));
  if (not m_tree) {
    throw std::runtime_error("Could not read tree");
  }
  setup_branches(std::make_index_sequence<std::tuple_size<Tuple>::value>());
}

template<typename NamedTuple>
inline NamedTupleRootReader<NamedTuple>::NamedTupleRootReader(
  TDirectory* dir, const std::string& tree_name)
  : m_file(nullptr) // no file since it is not owned by the writer
  , m_tree(nullptr)
  , m_next(0) {
  if (not dir) {
    throw std::runtime_error("Invalid input directory given");
  }
  m_tree = static_cast<TTree*>(dir->Get(tree_name.c_str()));
  if (not m_tree) {
    throw std::runtime_error("Could not read tree");
  }
  setup_branches(std::make_index_sequence<std::tuple_size<Tuple>::value>());
}

namespace io_root_impl {

// WARNING this is a hack to get around inconsistent ROOT types for 8bit chars
// and 64bit intengers compared to the stdint types.
__attribute__((unused)) inline ULong64_t*
get_address(uint64_t& x) {
  static_assert(
    sizeof(ULong64_t) == sizeof(uint64_t), "Inconsistent type sizes");
  return reinterpret_cast<ULong64_t*>(&x);
}
__attribute__((unused)) inline char*
get_address(int8_t& x) {
  static_assert(sizeof(char) == sizeof(int8_t), "Inconsistent type sizes");
  return reinterpret_cast<char*>(&x);
}
__attribute__((unused)) inline Long64_t*
get_address(int64_t& x) {
  static_assert(sizeof(Long64_t) == sizeof(int64_t), "Inconsistent type sizes");
  return reinterpret_cast<Long64_t*>(&x);
}
template<typename T>
inline T*
get_address(T& x) {
  return &x;
}

} // namespace io_root_impl

template<typename NamedTuple>
template<std::size_t... I>
inline void
NamedTupleRootReader<NamedTuple>::setup_branches(std::index_sequence<I...>) {
  static_assert(
    sizeof...(I) == std::tuple_size<Tuple>::value, "Something is very wrong");

  using std::get;

  // construct leaf names w/ type info
  std::array<std::string, sizeof...(I)> names = NamedTuple::names();
  // construct branches
  (void)std::array<Int_t, sizeof...(I)>{m_tree->SetBranchAddress(
    names[I].c_str(), io_root_impl::get_address(get<I>(m_data)))...};
}

template<typename NamedTuple>
inline NamedTupleRootReader<NamedTuple>::~NamedTupleRootReader() {
  // reader owns the file
  if (m_file) {
    m_file->Close();
    delete m_file;
  }
}

template<typename NamedTuple>
inline bool
NamedTupleRootReader<NamedTuple>::read(NamedTuple& record) {
  auto ret = m_tree->GetEntry(m_next);
  // i/o error occured
  if (ret < 0) {
    throw std::runtime_error("Could not read entry");
  }
  // the entry does not exist, probably end-of-file reached
  if (ret == 0) {
    return false;
  }
  // GetEntry(...) has already filled the local buffer
  record = m_data;
  m_next += 1;
  return true;
}

} // namespace dfe
