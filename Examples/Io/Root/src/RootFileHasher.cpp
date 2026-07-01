// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootFileHasher.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <format>
#include <iostream>
#include <memory>
#include <set>
#include <span>
#include <stdexcept>
#include <vector>

#include <boost/algorithm/string/trim.hpp>

#include <TBranch.h>
#include <TFile.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#if defined(ACTS_ROOT_FILE_HASHER_USE_BOOST_HASH2)
#include <boost/hash2/sha2.hpp>
#else
#include <TMD5.h>
#endif

namespace ActsExamples {

namespace {

#if defined(ACTS_ROOT_FILE_HASHER_USE_BOOST_HASH2)
constexpr std::size_t kDigestSize = 32;
#else
constexpr std::size_t kDigestSize = 16;
#endif

using Digest = std::array<unsigned char, kDigestSize>;

/// Incremental hasher backed by Boost.Hash2 SHA-256 (preferred) or ROOT's
/// TMD5 (fallback for Boost < 1.86). Only the backend member and the two
/// method bodies differ between the two configurations.
class Hasher {
 public:
  void update(std::span<const std::byte> data) {
#if defined(ACTS_ROOT_FILE_HASHER_USE_BOOST_HASH2)
    m_hasher.update(reinterpret_cast<const unsigned char*>(data.data()),
                    data.size());
#else
    m_hasher.Update(reinterpret_cast<const UChar_t*>(data.data()), data.size());
#endif
  }

  Digest finalize() {
    Digest digest{};
#if defined(ACTS_ROOT_FILE_HASHER_USE_BOOST_HASH2)
    boost::hash2::digest<kDigestSize> result = m_hasher.result();
    std::copy(result.data(), result.data() + kDigestSize, digest.begin());
#else
    m_hasher.Final(digest.data());
#endif
    return digest;
  }

 private:
#if defined(ACTS_ROOT_FILE_HASHER_USE_BOOST_HASH2)
  boost::hash2::sha2_256 m_hasher;
#else
  TMD5 m_hasher;
#endif
};

/// Compute the digest of a contiguous byte buffer.
Digest digestOf(std::span<const std::byte> data) {
  Hasher hasher;
  hasher.update(data);
  return hasher.finalize();
}

/// Render a digest as a lowercase hex string.
std::string toHex(const Digest& digest) {
  static const char* const hexDigits = "0123456789abcdef";
  std::string out;
  out.reserve(2 * digest.size());
  for (unsigned char byte : digest) {
    out.push_back(hexDigits[byte >> 4]);
    out.push_back(hexDigits[byte & 0x0f]);
  }
  return out;
}

/// Append the raw object representation of `x` to a byte buffer.
template <typename T>
void appendBytes(std::vector<std::byte>& out, const T& x) {
  const auto* p = reinterpret_cast<const std::byte*>(&x);
  out.insert(out.end(), p, p + sizeof(T));
}

/// Type-erased reader appending the raw bytes of one branch's current-entry
/// value(s) to a buffer.
struct IBranchReader {
  virtual ~IBranchReader() = default;
  /// Append the bytes of the current entry to `out`.
  virtual void append(std::vector<std::byte>& out) = 0;
  /// Whether the underlying branch was matched successfully. Only meaningful
  /// after the first entry has been loaded.
  virtual bool valid() const = 0;
};

template <typename T>
bool setupOk(const T& proxy) {
  return proxy.GetSetupStatus() ==
         ROOT::Internal::TTreeReaderValueBase::kSetupMatch;
}

/// Reader for a fundamental scalar branch.
template <typename T>
struct ScalarReader final : IBranchReader {
  TTreeReaderValue<T> value;
  ScalarReader(TTreeReader& reader, const char* name) : value(reader, name) {}
  bool valid() const override { return setupOk(value); }
  void append(std::vector<std::byte>& out) override { appendBytes(out, *value); }
};

/// Reader for a 1D array branch (`std::vector<T>` or a C-style array leaf).
template <typename T>
struct ArrayReader final : IBranchReader {
  TTreeReaderArray<T> array;
  ArrayReader(TTreeReader& reader, const char* name) : array(reader, name) {}
  bool valid() const override { return setupOk(array); }
  void append(std::vector<std::byte>& out) override {
    for (T x : array) {
      appendBytes(out, x);
    }
  }
};

/// Reader for a nested `std::vector<std::vector<T>>` branch.
template <typename T>
struct NestedReader final : IBranchReader {
  TTreeReaderValue<std::vector<std::vector<T>>> value;
  NestedReader(TTreeReader& reader, const char* name) : value(reader, name) {}
  bool valid() const override { return setupOk(value); }
  void append(std::vector<std::byte>& out) override {
    for (const auto& inner : *value) {
      for (T x : inner) {
        appendBytes(out, x);
      }
    }
  }
};

enum class Kind { Scalar, Array, Nested };

template <typename T>
std::unique_ptr<IBranchReader> makeReader(Kind kind, TTreeReader& reader,
                                          const char* name) {
  switch (kind) {
    case Kind::Scalar:
      return std::make_unique<ScalarReader<T>>(reader, name);
    case Kind::Array:
      return std::make_unique<ArrayReader<T>>(reader, name);
    case Kind::Nested:
      return std::make_unique<NestedReader<T>>(reader, name);
  }
  return nullptr;
}

/// Map a (ROOT or canonical C++) element type name to a typed reader. The
/// instantiated C++ type must match the branch's stored type for ROOT to bind
/// the reader, so we keep the names aligned. Returns nullptr for unsupported
/// types.
std::unique_ptr<IBranchReader> dispatch(const std::string& type, Kind kind,
                                        TTreeReader& reader, const char* name) {
  if (type == "Float_t" || type == "float") {
    return makeReader<float>(kind, reader, name);
  }
  if (type == "Double_t" || type == "double") {
    return makeReader<double>(kind, reader, name);
  }
  if (type == "Int_t" || type == "int") {
    return makeReader<int>(kind, reader, name);
  }
  if (type == "UInt_t" || type == "unsigned int" || type == "unsigned") {
    return makeReader<unsigned int>(kind, reader, name);
  }
  if (type == "Long64_t" || type == "long long") {
    return makeReader<long long>(kind, reader, name);
  }
  if (type == "ULong64_t" || type == "unsigned long long") {
    return makeReader<unsigned long long>(kind, reader, name);
  }
  if (type == "Long_t" || type == "long") {
    return makeReader<long>(kind, reader, name);
  }
  if (type == "ULong_t" || type == "unsigned long") {
    return makeReader<unsigned long>(kind, reader, name);
  }
  if (type == "Short_t" || type == "short") {
    return makeReader<short>(kind, reader, name);
  }
  if (type == "UShort_t" || type == "unsigned short") {
    return makeReader<unsigned short>(kind, reader, name);
  }
  if (type == "Char_t" || type == "char") {
    return makeReader<char>(kind, reader, name);
  }
  if (type == "UChar_t" || type == "unsigned char") {
    return makeReader<unsigned char>(kind, reader, name);
  }
  if (type == "Bool_t" || type == "bool") {
    return makeReader<bool>(kind, reader, name);
  }
  return nullptr;
}

/// Return the first top-level template argument of a type name, e.g.
/// "vector<float>" -> "float" and "vector<vector<int> >" -> "vector<int>".
std::string firstTemplateArg(const std::string& s) {
  auto lt = s.find('<');
  if (lt == std::string::npos) {
    return "";
  }
  int depth = 0;
  std::size_t start = lt + 1;
  for (std::size_t i = lt; i < s.size(); ++i) {
    char c = s[i];
    if (c == '<') {
      ++depth;
      if (depth == 1) {
        start = i + 1;
      }
    } else if (c == '>') {
      --depth;
      if (depth == 0) {
        return boost::algorithm::trim_copy(s.substr(start, i - start));
      }
    } else if (c == ',' && depth == 1) {
      return boost::algorithm::trim_copy(s.substr(start, i - start));
    }
  }
  return "";
}

struct BranchInfo {
  std::string name;
  Kind kind = Kind::Scalar;
  std::string elementType;
};

/// Determine the kind and element type of every branch in a tree.
std::vector<BranchInfo> describeBranches(TTree& tree) {
  std::vector<BranchInfo> infos;
  TObjArray* branches = tree.GetListOfBranches();
  if (branches == nullptr) {
    return infos;
  }
  for (int i = 0; i < branches->GetEntriesFast(); ++i) {
    auto* branch = dynamic_cast<TBranch*>(branches->At(i));
    if (branch == nullptr) {
      continue;
    }
    BranchInfo info;
    info.name = branch->GetName();

    std::string className = branch->GetClassName();
    if (!className.empty()) {
      // STL collection branch.
      std::string arg = firstTemplateArg(className);
      if (arg.rfind("vector<", 0) == 0) {
        info.kind = Kind::Nested;
        info.elementType = firstTemplateArg(arg);
      } else {
        info.kind = Kind::Array;
        info.elementType = arg;
      }
    } else {
      // Fundamental leaf branch.
      auto* leaf = dynamic_cast<TLeaf*>(branch->GetListOfLeaves()->At(0));
      if (leaf == nullptr) {
        continue;
      }
      info.elementType = leaf->GetTypeName();
      bool isArray =
          (leaf->GetLeafCount() != nullptr) || (leaf->GetLenStatic() > 1);
      info.kind = isArray ? Kind::Array : Kind::Scalar;
    }
    infos.push_back(std::move(info));
  }
  std::sort(
      infos.begin(), infos.end(),
      [](const BranchInfo& a, const BranchInfo& b) { return a.name < b.name; });
  return infos;
}

struct Reader {
  std::string name;
  std::unique_ptr<IBranchReader> reader;
  std::unique_ptr<Hasher> hasher;  // only used in non-order-invariant mode
};

/// Compute the digest of a single tree.
Digest hashTree(TTree& tree, bool orderInvariant) {
  std::vector<BranchInfo> infos = describeBranches(tree);

  TTreeReader treeReader(&tree);
  std::vector<Reader> readers;
  for (const BranchInfo& info : infos) {
    auto branchReader =
        dispatch(info.elementType, info.kind, treeReader, info.name.c_str());
    if (branchReader == nullptr) {
      std::cerr << "RootFileHasher: unsupported type '" << info.elementType
                << "' for branch '" << info.name << "' (skipped)" << std::endl;
      continue;
    }
    readers.push_back({info.name, std::move(branchReader),
                       orderInvariant ? nullptr : std::make_unique<Hasher>()});
  }

  std::vector<Digest> rowDigests;
  std::vector<std::byte> buffer;
  bool pruned = false;
  while (treeReader.Next()) {
    if (!pruned) {
      // Drop branches that could not be bound (e.g. unmatched type). This can
      // only be checked once an entry has been loaded.
      std::vector<Reader> kept;
      for (auto& r : readers) {
        if (r.reader->valid()) {
          kept.push_back(std::move(r));
        } else {
          std::cerr << "RootFileHasher: branch '" << r.name
                    << "' could not be read (skipped)" << std::endl;
        }
      }
      readers = std::move(kept);
      pruned = true;
    }

    if (orderInvariant) {
      buffer.clear();
      for (auto& r : readers) {
        r.reader->append(buffer);
      }
      rowDigests.push_back(digestOf(buffer));
    } else {
      for (auto& r : readers) {
        buffer.clear();
        r.reader->append(buffer);
        r.hasher->update(buffer);
      }
    }
  }

  Hasher treeHash;
  std::string names;
  for (const auto& r : readers) {
    names += r.name;
  }
  treeHash.update(std::as_bytes(std::span(names)));

  if (orderInvariant) {
    std::sort(rowDigests.begin(), rowDigests.end());
    for (const Digest& d : rowDigests) {
      treeHash.update(std::as_bytes(std::span(d)));
    }
  } else {
    for (auto& r : readers) {
      Digest d = r.hasher->finalize();
      treeHash.update(std::as_bytes(std::span(d)));
    }
  }

  return treeHash.finalize();
}

}  // namespace

std::string hashRootFile(const std::filesystem::path& path,
                         bool orderInvariant) {
  TFile file(path.c_str(), "READ");
  if (file.IsZombie()) {
    throw std::runtime_error(
        std::format("RootFileHasher: could not open '{}'", path.string()));
  }

  // Collect unique top-level key names in sorted order. Every key name is
  // folded into the hash (so structural changes are detected even for files
  // that only contain histograms), and the numeric content of TTrees is hashed
  // on top.
  std::set<std::string> keyNames;
  TList* keys = file.GetListOfKeys();
  if (keys != nullptr) {
    for (int i = 0; i < keys->GetEntries(); ++i) {
      auto* key = dynamic_cast<TKey*>(keys->At(i));
      if (key != nullptr) {
        keyNames.insert(key->GetName());
      }
    }
  }

  Hasher global;
  for (const std::string& name : keyNames) {
    global.update(std::as_bytes(std::span(name)));

    auto* tree = dynamic_cast<TTree*>(file.Get(name.c_str()));
    if (tree == nullptr) {
      continue;
    }
    Digest treeDigest = hashTree(*tree, orderInvariant);
    global.update(std::as_bytes(std::span(treeDigest)));
  }

  return toHex(global.finalize());
}

}  // namespace ActsExamples
