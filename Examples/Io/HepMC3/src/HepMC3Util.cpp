// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/HepMC3/HepMC3Util.hpp"

#include "Acts/Utilities/ScopedTimer.hpp"

#include <stdexcept>

#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>

namespace ActsExamples {

namespace {
template <typename T>
void mergeEventsImpl(HepMC3::GenEvent& event, std::span<T> genEvents,
                     const Acts::Logger& logger) {
  Acts::AveragingScopedTimer mergeTimer("Merging HepMC3 events", logger(),
                                        Acts::Logging::DEBUG);

  std::vector<std::shared_ptr<HepMC3::GenParticle>> particles;

  // Loop once to find the total size we'll need
  std::size_t nParticles = 0;
  std::size_t nVertices = 0;
  for (const auto& genEvent : genEvents) {
    nParticles += genEvent->particles().size();
    nVertices += genEvent->vertices().size();
  }

  event.reserve(nParticles, nVertices);

  for (const auto& genEvent : genEvents) {
    auto sample = mergeTimer.sample();
    particles.clear();
    particles.reserve(genEvent->particles().size());

    auto copyAttributes = [&](const auto& src, auto& dst) {
      for (const auto& attr : src.attribute_names()) {
        auto value = src.attribute_as_string(attr);
        dst.add_attribute(attr,
                          std::make_shared<HepMC3::StringAttribute>(value));
      }
    };

    copyAttributes(*genEvent, event);

    // Add to combined event
    for (const auto& srcParticle : genEvent->particles()) {
      if (srcParticle->id() - 1 != static_cast<int>(particles.size())) {
        throw std::runtime_error("Particle id is not consecutive");
      }
      auto particle = std::make_shared<HepMC3::GenParticle>();
      particle->set_momentum(srcParticle->momentum());
      particle->set_generated_mass(srcParticle->generated_mass());
      particle->set_pid(srcParticle->pid());
      particle->set_status(srcParticle->status());

      particles.push_back(particle);
      event.add_particle(particle);

      copyAttributes(*srcParticle, *particle);
    }

    for (const auto& srcVertex : genEvent->vertices()) {
      auto vertex = std::make_shared<HepMC3::GenVertex>(srcVertex->position());
      vertex->set_status(srcVertex->status());

      event.add_vertex(vertex);

      copyAttributes(*srcVertex, *vertex);

      for (const auto& srcParticle : srcVertex->particles_in()) {
        const auto& particle = particles.at(srcParticle->id() - 1);
        vertex->add_particle_in(particle);
      }
      for (const auto& srcParticle : srcVertex->particles_out()) {
        const auto& particle = particles.at(srcParticle->id() - 1);
        vertex->add_particle_out(particle);
      }
    }
  }
}
}  // namespace

void HepMC3Util::mergeEvents(HepMC3::GenEvent& event,
                             std::span<const HepMC3::GenEvent*> genEvents,
                             const Acts::Logger& logger) {
  mergeEventsImpl(event, genEvents, logger);
}

void HepMC3Util::mergeEvents(
    HepMC3::GenEvent& event,
    std::span<std::shared_ptr<const HepMC3::GenEvent>> genEvents,
    const Acts::Logger& logger) {
  mergeEventsImpl(event, genEvents, logger);
}

std::string_view HepMC3Util::compressionExtension(Compression compression) {
  switch (compression) {
    using enum Compression;
    case none:
      return "";
    case zlib:
      return ".gz";
    case lzma:
      return ".xz";
    case bzip2:
      return ".bz2";
    case zstd:
      return ".zst";
    default:
      throw std::invalid_argument{"Unknown compression value"};
  }
}

std::span<const HepMC3Util::Compression>
HepMC3Util::availableCompressionModes() {
  using enum Compression;
  static const auto values = []() -> std::vector<HepMC3Util::Compression> {
    return {
        none,
#ifdef HEPMC3_Z_SUPPORT
        zlib,
#endif
#ifdef HEPMC3_LZMA_SUPPORT
        lzma,
#endif
#ifdef HEPMC3_BZ2_SUPPORT
        bzip2,
#endif
#ifdef HEPMC3_ZSTD_SUPPORT
        zstd,
#endif
    };
  }();
  return values;
}

std::ostream& HepMC3Util::operator<<(std::ostream& os,
                                     HepMC3Util::Compression compression) {
  switch (compression) {
    using enum HepMC3Util::Compression;
    case none:
      return os << "none";
    case zlib:
      return os << "zlib";
    case lzma:
      return os << "lzma";
    case bzip2:
      return os << "bzip2";
    case zstd:
      return os << "zstd";
    default:
      throw std::invalid_argument{"Unknown compression value"};
  }
}

std::ostream& HepMC3Util::operator<<(std::ostream& os,
                                     HepMC3Util::Format format) {
  switch (format) {
    using enum HepMC3Util::Format;
    case ascii:
      return os << "ascii";
    case root:
      return os << "root";
    default:
      throw std::invalid_argument{"Unknown format value"};
  }
}

std::span<const HepMC3Util::Format> HepMC3Util::availableFormats() {
  using enum Format;
  static const auto values = []() -> std::vector<HepMC3Util::Format> {
    return {
        ascii,
#ifdef HEPMC3_ROOT_SUPPORT
        root,
#endif
    };
  }();
  return values;
}

HepMC3Util::Format HepMC3Util::formatFromFilename(std::string_view filename) {
  using enum Format;

  for (auto compression : availableCompressionModes()) {
    auto ext = compressionExtension(compression);

    if (filename.ends_with(".hepmc3" + std::string(ext)) ||
        filename.ends_with(".hepmc" + std::string(ext))) {
      return ascii;
    }
  }
  if (filename.ends_with(".root")) {
    return root;
  }

  throw std::invalid_argument{"Unknown format extension: " +
                              std::string{filename}};
}

}  // namespace ActsExamples
