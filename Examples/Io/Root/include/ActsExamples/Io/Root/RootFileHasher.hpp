// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <filesystem>
#include <string>

namespace ActsExamples {

/// Compute a hash of the numeric content of a ROOT file.
///
/// This is a C++ counterpart to Examples/Python/tests/helpers/hash_root.py and
/// is meant for regression testing: it produces a stable fingerprint of all
/// TTree branch values in a file. The hash is NOT byte-compatible with the
/// Python implementation; only the guarantees are the same:
///   - determinism for identical content,
///   - sensitivity to any change in numeric content,
///   - (in order-invariant mode) invariance under reordering of tree entries.
///
/// The hash is computed from the in-memory representation of the numeric types
/// and is therefore tied to the platform (endianness, type widths); it is not
/// portable across differing architectures.
///
/// The underlying digest algorithm is selected at build time: when Boost
/// >= 1.86 is available, Boost.Hash2 SHA-256 is used (64 hex characters);
/// otherwise it falls back to ROOT's TMD5 / MD5 (32 hex characters). The two
/// backends produce different digests, so hashes are only comparable within a
/// single build configuration.
///
/// @param path Path to the input ROOT file.
/// @param orderInvariant If true (default), the per-entry digests of each tree
///        are sorted before being combined, so the result does not depend on
///        the ordering of entries. If false, contents are hashed in storage
///        order (faster, but sensitive to reordering).
/// @return The hash as a lowercase hex string.
/// @throws std::runtime_error if the file cannot be opened.
std::string hashRootFile(const std::filesystem::path& path,
                         bool orderInvariant = true);

}  // namespace ActsExamples
