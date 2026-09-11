// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Json/detail/JsonIo.hpp"

#include <algorithm>
#include <array>
#include <format>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>

#ifdef ACTS_JSON_ZSTD_SUPPORT
#include <zstd.h>
#endif

namespace {

using Acts::detail::JsonCompression;
using Acts::detail::JsonEncoding;
using Acts::detail::JsonFileFormat;

/// zstd frame magic number 0xFD2FB528, little endian as it appears on disk.
constexpr std::array<std::byte, 4> kZstdMagic{std::byte{0x28}, std::byte{0xB5},
                                              std::byte{0x2F}, std::byte{0xFD}};

bool hasZstdMagic(std::span<const std::byte> data) {
  return data.size() >= kZstdMagic.size() &&
         std::equal(kZstdMagic.begin(), kZstdMagic.end(), data.begin());
}

/// Copy any contiguous byte-like range into a byte buffer.
template <typename container_t>
std::vector<std::byte> toBytes(const container_t& container) {
  auto bytes = std::as_bytes(std::span{container});
  return {bytes.begin(), bytes.end()};
}

std::string_view asStringView(std::span<const std::byte> data) {
  return {reinterpret_cast<const char*>(data.data()), data.size()};
}

std::vector<std::byte> compressZstd(std::span<const std::byte> data,
                                    int level) {
#ifdef ACTS_JSON_ZSTD_SUPPORT
  std::vector<std::byte> out(ZSTD_compressBound(data.size()));
  std::size_t written =
      ZSTD_compress(out.data(), out.size(), data.data(), data.size(), level);
  if (ZSTD_isError(written) != 0u) {
    throw std::runtime_error(
        std::format("zstd compression failed: {}", ZSTD_getErrorName(written)));
  }
  out.resize(written);
  return out;
#else
  static_cast<void>(data);
  static_cast<void>(level);
  throw std::runtime_error(
      "zstd compressed output was requested, but ACTS was built without zstd "
      "support. Reconfigure with zstd available to write compressed files.");
#endif
}

std::vector<std::byte> decompressZstd(std::span<const std::byte> data,
                                      const std::filesystem::path& origin) {
#ifdef ACTS_JSON_ZSTD_SUPPORT
  unsigned long long size = ZSTD_getFrameContentSize(data.data(), data.size());
  if (size == ZSTD_CONTENTSIZE_ERROR) {
    throw std::runtime_error(std::format(
        "'{}' starts with a zstd magic number but is not a valid zstd frame",
        origin.string()));
  }
  if (size == ZSTD_CONTENTSIZE_UNKNOWN) {
    throw std::runtime_error(std::format(
        "'{}' is a zstd stream without a declared content size, which is not "
        "supported here; recompress it as a single frame",
        origin.string()));
  }
  std::vector<std::byte> out(size);
  std::size_t written =
      ZSTD_decompress(out.data(), out.size(), data.data(), data.size());
  if (ZSTD_isError(written) != 0u) {
    throw std::runtime_error(
        std::format("'{}' could not be decompressed, the zstd frame is likely "
                    "truncated or corrupt: {}",
                    origin.string(), ZSTD_getErrorName(written)));
  }
  out.resize(written);
  return out;
#else
  static_cast<void>(data);
  throw std::runtime_error(std::format(
      "'{}' is zstd compressed, but ACTS was built without zstd support. "
      "Reconfigure with zstd available to read this file.",
      origin.string()));
#endif
}

}  // namespace

bool Acts::detail::zstdSupported() {
#ifdef ACTS_JSON_ZSTD_SUPPORT
  return true;
#else
  return false;
#endif
}

Acts::detail::JsonFileFormat Acts::detail::jsonFormatFromPath(
    const std::filesystem::path& path) {
  std::string name = path.filename().string();
  JsonFileFormat format{};

  if (name.ends_with(".zst")) {
    format.compression = JsonCompression::Zstd;
    name.resize(name.size() - std::string_view{".zst"}.size());
  }

  if (name.ends_with(".json")) {
    format.encoding = JsonEncoding::Text;
  } else if (name.ends_with(".cbor")) {
    format.encoding = JsonEncoding::Cbor;
  } else {
    throw std::invalid_argument(std::format(
        "Cannot deduce the output format of '{}': expected one of '.json', "
        "'.cbor', '.json.zst' or '.cbor.zst'",
        path.string()));
  }

  return format;
}

std::vector<std::byte> Acts::detail::encodeJson(const nlohmann::json& payload,
                                                const JsonFileFormat& format,
                                                unsigned indentation,
                                                int compressionLevel) {
  std::vector<std::byte> encoded =
      format.encoding == JsonEncoding::Cbor
          ? toBytes(nlohmann::json::to_cbor(payload))
          : toBytes(payload.dump(static_cast<int>(indentation)));

  if (format.compression == JsonCompression::Zstd) {
    return compressZstd(encoded, compressionLevel);
  }
  return encoded;
}

nlohmann::json Acts::detail::decodeJson(std::span<const std::byte> data,
                                        const std::filesystem::path& origin) {
  std::vector<std::byte> decompressed;
  bool wasCompressed = hasZstdMagic(data);
  if (wasCompressed) {
    decompressed = decompressZstd(data, origin);
    data = decompressed;
  }

  auto isSpace = [](std::byte b) {
    char c = static_cast<char>(b);
    return c == ' ' || c == '\t' || c == '\n' || c == '\r';
  };
  auto it = std::ranges::find_if_not(data, isSpace);
  if (it == data.end()) {
    throw std::runtime_error(
        std::format("'{}' contains no payload to decode", origin.string()));
  }

  // JSON text always starts with '{' or '['. Our CBOR payloads are always a
  // map or an array at the top level, i.e. major type 5 or 4, which cannot
  // collide with either of those two ASCII bytes.
  auto lead = static_cast<unsigned char>(*it);
  bool isText = lead == '{' || lead == '[';
  bool isCbor = (lead & 0xE0) == 0xA0 || (lead & 0xE0) == 0x80;

  // A decompressed payload that does not parse points at the compressed file,
  // so say which stage the content came from.
  std::string_view stage = wasCompressed ? " (after decompression)" : "";

  if (!isText && !isCbor) {
    throw std::runtime_error(std::format(
        "'{}'{} is neither JSON text nor CBOR, leading byte is 0x{:02X}",
        origin.string(), stage, lead));
  }

  try {
    if (isText) {
      return nlohmann::json::parse(asStringView(data));
    }
    return nlohmann::json::from_cbor(data);
  } catch (const nlohmann::json::exception& e) {
    throw std::runtime_error(std::format(
        "Failed to decode {} payload from '{}'{}: {}", isText ? "JSON" : "CBOR",
        origin.string(), stage, e.what()));
  }
}

void Acts::detail::writeJsonFile(const std::filesystem::path& path,
                                 const nlohmann::json& payload,
                                 unsigned indentation, int compressionLevel) {
  std::vector<std::byte> encoded = encodeJson(payload, jsonFormatFromPath(path),
                                              indentation, compressionLevel);

  std::ofstream ofs{path, std::ios::binary};
  if (!ofs.good()) {
    throw std::runtime_error(
        std::format("Cannot open '{}' for writing", path.string()));
  }
  ofs.write(reinterpret_cast<const char*>(encoded.data()),
            static_cast<std::streamsize>(encoded.size()));
  if (!ofs.good()) {
    throw std::runtime_error(
        std::format("Failed to write '{}'", path.string()));
  }
}

nlohmann::json Acts::detail::readJsonFile(const std::filesystem::path& path) {
  if (!std::filesystem::exists(path)) {
    throw std::invalid_argument(
        std::format("File '{}' does not exist", path.string()));
  }

  std::ifstream ifs{path, std::ios::binary};
  if (!ifs.good()) {
    throw std::invalid_argument(std::format("Cannot open '{}'", path.string()));
  }

  ifs.seekg(0, std::ios::end);
  auto size = ifs.tellg();
  if (size < 0) {
    throw std::runtime_error(
        std::format("Cannot determine the size of '{}'", path.string()));
  }
  ifs.seekg(0, std::ios::beg);

  std::vector<std::byte> data(static_cast<std::size_t>(size));
  ifs.read(reinterpret_cast<char*>(data.data()),
           static_cast<std::streamsize>(data.size()));
  if (ifs.bad()) {
    throw std::runtime_error(std::format("Failed to read '{}'", path.string()));
  }

  return decodeJson(data, path);
}
