// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3Writer.hpp"

#include "ActsExamples/Utilities/Paths.hpp"

#include <filesystem>

#include <HepMC3/Version.h>
#include <HepMC3/WriterAscii.h>

#if HEPMC3_VERSION_CODE == 3002007
#include "./CompressedIO.h"
#endif

#ifdef HEPMC3_USE_COMPRESSION
#include <HepMC3/WriterGZ.h>
#endif

namespace ActsExamples {

HepMC3Writer::HepMC3Writer(const Config& config, Acts::Logging::Level level)
    : WriterT(config.inputEvent, "HepMC3Writer", level), m_cfg(config) {
  if (m_cfg.outputPath.empty()) {
    throw std::invalid_argument("Missing output file path");
  }

  if (std::ranges::none_of(HepMC3Util::availableCompressionModes(),
                           [this](HepMC3Util::Compression c) {
                             return c == this->m_cfg.compression;
                           })) {
    std::stringstream ss;
    ss << "Unsupported compression mode: " << m_cfg.compression;
    throw std::invalid_argument(ss.str());
  }

  if (!m_cfg.perEvent) {
    auto absolute = std::filesystem::absolute(m_cfg.outputPath);
    if (std::filesystem::exists(absolute) &&
        std::filesystem::is_directory(absolute)) {
      throw std::invalid_argument("Output path is a directory: " +
                                  absolute.string());
    }

    if (!std::filesystem::exists(absolute.parent_path())) {
      throw std::invalid_argument("Directory to write into does not exist: " +
                                  absolute.parent_path().string());
    }
    // Create a single file writer
    m_writer = createWriter(m_cfg.outputPath);
  }
}

std::unique_ptr<HepMC3::Writer> HepMC3Writer::createWriter(
    const std::filesystem::path& target) {
  std::filesystem::path path =
      target.string() +
      std::string{HepMC3Util::compressionExtension(m_cfg.compression)};

  switch (m_cfg.compression) {
    case HepMC3Util::Compression::none:
      return std::make_unique<HepMC3::WriterAscii>(path);
#ifdef HEPMC3_USE_COMPRESSION
    case HepMC3Util::Compression::zlib:
      return std::make_unique<
          HepMC3::WriterGZ<HepMC3::WriterAscii, HepMC3::Compression::z>>(path);
    case HepMC3Util::Compression::lzma:
      return std::make_unique<
          HepMC3::WriterGZ<HepMC3::WriterAscii, HepMC3::Compression::lzma>>(
          path);
    case HepMC3Util::Compression::bzip2:
      return std::make_unique<
          HepMC3::WriterGZ<HepMC3::WriterAscii, HepMC3::Compression::bz2>>(
          path);
    case HepMC3Util::Compression::zstd:
      return std::make_unique<
          HepMC3::WriterGZ<HepMC3::WriterAscii, HepMC3::Compression::zstd>>(
          path);
#endif
    default:
      throw std::invalid_argument{"Unknown compression value"};
  }
}

HepMC3Writer::~HepMC3Writer() = default;

ProcessCode HepMC3Writer::writeT(
    const AlgorithmContext& ctx,
    const std::shared_ptr<HepMC3::GenEvent>& event) {
  ACTS_VERBOSE("Writing " << event->particles().size() << " particles to "
                          << m_cfg.outputPath);

  auto write = [&event](HepMC3::Writer& writer) {
    writer.write_event(*event);
    if (writer.failed()) {
      return ProcessCode::ABORT;
    }
    return ProcessCode::SUCCESS;
  };

  if (m_cfg.perEvent) {
    std::filesystem::path perEventFile =
        perEventFilepath(m_cfg.outputPath.parent_path(),
                         m_cfg.outputPath.filename().string(), ctx.eventNumber);

    ACTS_VERBOSE("Writing per-event file " << perEventFile);
    auto writer = createWriter(perEventFile);
    auto result = write(*writer);
    writer->close();
    return result;
  } else {
    ACTS_VERBOSE("Writing to single file " << m_cfg.outputPath);
    // Take the lock until the end of the function
    std::scoped_lock lock(m_mutex);
    return write(*m_writer);
  }
}

ProcessCode HepMC3Writer::finalize() {
  ACTS_VERBOSE("Finalizing HepMC3Writer");
  if (m_writer) {
    m_writer->close();
  }
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
