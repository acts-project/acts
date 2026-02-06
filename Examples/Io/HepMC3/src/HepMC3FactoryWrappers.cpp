// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3Util.hpp"

#include <stdexcept>

// Include problematic HepMC3 headers only in this translation unit
// to avoid multiple definition linker errors in HepMC3 < 3.3.0
#include <HepMC3/ReaderFactory.h>
#include <HepMC3/Writer.h>
#include <HepMC3/WriterAscii.h>

#ifdef HEPMC3_USE_COMPRESSION
#include <HepMC3/WriterGZ.h>
#endif

#ifdef HEPMC3_ROOT_SUPPORT
#include <HepMC3/WriterRootTree.h>
#endif

namespace ActsExamples::HepMC3Util {

std::shared_ptr<HepMC3::Reader> deduceReader(const std::string& filename) {
  return HepMC3::deduce_reader(filename);
}

std::unique_ptr<HepMC3::Writer> createWriter(const std::filesystem::path& path,
                                             Format format,
                                             Compression compression) {
  if (format == Format::root) {
#ifdef HEPMC3_ROOT_SUPPORT
    if (compression != Compression::none) {
      throw std::invalid_argument("Compression not supported for ROOT format");
    }
    return std::make_unique<HepMC3::WriterRootTree>(path.string());
#else
    throw std::runtime_error("ROOT support not enabled in HepMC3");
#endif
  } else {
    // ASCII format
    switch (compression) {
      using enum Compression;
      case none:
        return std::make_unique<HepMC3::WriterAscii>(path.string());
#ifdef HEPMC3_USE_COMPRESSION
      case zlib:
        return std::make_unique<
            HepMC3::WriterGZ<HepMC3::WriterAscii, HepMC3::Compression::z>>(
            path.string());
      case lzma:
        return std::make_unique<
            HepMC3::WriterGZ<HepMC3::WriterAscii, HepMC3::Compression::lzma>>(
            path.string());
      case bzip2:
        return std::make_unique<
            HepMC3::WriterGZ<HepMC3::WriterAscii, HepMC3::Compression::bz2>>(
            path.string());
      case zstd:
        return std::make_unique<
            HepMC3::WriterGZ<HepMC3::WriterAscii, HepMC3::Compression::zstd>>(
            path.string());
#endif
      default:
        throw std::invalid_argument("Unknown or unsupported compression type");
    }
  }
}

}  // namespace ActsExamples::HepMC3Util
