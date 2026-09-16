// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <filesystem>
#include <span>
#include <vector>

#include <nlohmann/json.hpp>

namespace Acts::detail {

/// How a JSON document is encoded on disk.
enum class JsonEncoding {
  /// Human readable JSON text
  Text,
  /// Binary CBOR, as produced by `nlohmann::json::to_cbor`
  Cbor,
};

/// Compression applied on top of the encoded payload.
enum class JsonCompression {
  /// The encoded payload is stored verbatim
  None,
  /// The encoded payload is wrapped in a single zstd frame
  Zstd,
};

/// The on-disk format of a JSON document: how it is encoded, and whether it
/// is compressed. The two are independent, so all four combinations are
/// valid.
struct JsonFileFormat {
  /// Payload encoding
  JsonEncoding encoding{JsonEncoding::Text};
  /// Payload compression
  JsonCompression compression{JsonCompression::None};
};

/// Whether this build can read and write zstd compressed payloads.
///
/// zstd is an optional dependency: when it is not found at configure time,
/// writing compressed output and reading a compressed file both throw.
///
/// @return true if zstd support was compiled in
bool zstdSupported();

/// Deduce the on-disk format from a file name.
///
/// The file name is the only specification of the output format, which keeps
/// it honest about its own contents:
///   - `.json`     text, uncompressed
///   - `.cbor`     CBOR, uncompressed
///   - `.json.zst` text, zstd compressed
///   - `.cbor.zst` CBOR, zstd compressed
///
/// @param path the file name to inspect
/// @throw std::invalid_argument if the extension is not one of the above
/// @return the deduced format
JsonFileFormat jsonFormatFromPath(const std::filesystem::path& path);

/// Encode a JSON payload into its on-disk representation.
///
/// @param payload the JSON document to encode
/// @param format the target encoding and compression
/// @param indentation indentation for text output, ignored for CBOR
/// @param compressionLevel zstd level, ignored when not compressing
/// @return the encoded bytes
std::vector<std::byte> encodeJson(const nlohmann::json& payload,
                                  const JsonFileFormat& format,
                                  unsigned indentation, int compressionLevel);

/// Decode a JSON payload, detecting encoding and compression from content.
///
/// Compression is detected from the zstd frame magic, and the encoding from
/// the leading byte of the (decompressed) payload, so the caller does not
/// need to know or declare the format. The file name is used only to make
/// error messages point at the offending file.
///
/// @param data the encoded bytes
/// @param origin file the data came from, used for diagnostics only
/// @return the decoded JSON document
nlohmann::json decodeJson(std::span<const std::byte> data,
                          const std::filesystem::path& origin = {});

/// Write a JSON payload to a file, in the format its name implies.
///
/// @param path destination, whose extension selects the format
/// @param payload the JSON document to write
/// @param indentation indentation for text output, ignored for CBOR
/// @param compressionLevel zstd level, ignored when not compressing
void writeJsonFile(const std::filesystem::path& path,
                   const nlohmann::json& payload, unsigned indentation,
                   int compressionLevel);

/// Read a JSON payload from a file, whatever format it is in.
///
/// @param path the file to read
/// @return the decoded JSON document
nlohmann::json readJsonFile(const std::filesystem::path& path);

}  // namespace Acts::detail
