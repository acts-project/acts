// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3Writer.hpp"

#include <HepMC3/WriterAscii.h>

namespace ActsExamples {

HepMC3AsciiWriter::HepMC3AsciiWriter(const Config& config,
                                     Acts::Logging::Level level)
    : WriterT(config.inputEvent, "HepMC3AsciiWriter", level), m_cfg(config) {
  if (m_cfg.outputPath.empty()) {
    throw std::invalid_argument("Missing output file path");
  }

  if (auto absolute = std::filesystem::absolute(m_cfg.outputPath);
      std::filesystem::exists(absolute.parent_path())) {
    throw std::invalid_argument("Output directory does not exist: " +
                                absolute.parent_path().string());
  }

  if (!m_cfg.perEvent) {
    // Create a single file writer
    m_writer = std::make_unique<HepMC3::WriterAscii>(m_cfg.outputPath);
  }
}

HepMC3AsciiWriter::~HepMC3AsciiWriter() = default;

ProcessCode HepMC3AsciiWriter::writeT(
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
    auto stem = m_cfg.outputPath.stem();
    auto ext = m_cfg.outputPath.extension();
    std::filesystem::path perEventFile =
        m_cfg.outputPath.parent_path() /
        (stem.string() + "_" + std::to_string(ctx.eventNumber) + ext.string());
    ACTS_VERBOSE("Writing per-event file " << perEventFile);
    HepMC3::WriterAscii writer(perEventFile);
    auto result = write(writer);
    writer.close();
    return result;
  } else {
    ACTS_VERBOSE("Writing to single file " << m_cfg.outputPath);
    // Take the lock until the end of the function
    std::scoped_lock lock(m_mutex);
    return write(*m_writer);
  }
}

ProcessCode HepMC3AsciiWriter::finalize() {
  ACTS_VERBOSE("Finalizing HepMC3AsciiWriter");
  if (m_writer) {
    m_writer->close();
  }
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
