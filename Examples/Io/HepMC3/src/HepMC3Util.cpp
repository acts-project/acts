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

std::shared_ptr<HepMC3::GenEvent> HepMC3Util::mergeEvents(
    std::span<const HepMC3::GenEvent*> genEvents, const Acts::Logger& logger) {
  Acts::AveragingScopedTimer mergeTimer("Merging generator events", logger(),
                                        Acts::Logging::DEBUG);

  std::vector<std::shared_ptr<HepMC3::GenParticle>> particles;

  auto event = std::make_shared<HepMC3::GenEvent>();
  event->set_units(HepMC3::Units::GEV, HepMC3::Units::MM);

  for (const auto& genEvent : genEvents) {
    auto sample = mergeTimer.sample();
    particles.clear();
    particles.reserve(genEvent->particles_size());

    auto copyAttributes = [&](const auto& src, auto& dst) {
      for (const auto& attr : src.attribute_names()) {
        auto value = src.attribute_as_string(attr);
        dst.add_attribute(attr,
                          std::make_shared<HepMC3::StringAttribute>(value));
      }
    };

    copyAttributes(*genEvent, *event);

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
      event->add_particle(particle);

      copyAttributes(*srcParticle, *particle);
    }

    for (const auto& srcVertex : genEvent->vertices()) {
      auto vertex = std::make_shared<HepMC3::GenVertex>(srcVertex->position());
      vertex->set_status(srcVertex->status());

      event->add_vertex(vertex);

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

  return event;
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

}  // namespace ActsExamples
