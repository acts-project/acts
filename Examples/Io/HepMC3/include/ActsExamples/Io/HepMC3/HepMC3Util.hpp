// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <ostream>
#include <span>

namespace HepMC3 {
class GenEvent;
class GenParticle;
}  // namespace HepMC3

namespace Acts {
class Logger;
}

namespace ActsExamples::HepMC3Util {

void mergeEvents(HepMC3::GenEvent& event,
                 std::span<const HepMC3::GenEvent*> genEvents,
                 const Acts::Logger& logger);

void mergeEvents(HepMC3::GenEvent& event,
                 std::span<std::shared_ptr<const HepMC3::GenEvent>> genEvents,
                 const Acts::Logger& logger);

constexpr int kBeamParticleStatus = 4;
constexpr int kUndecayedParticleStatus = 1;
constexpr int kDecayedParticleStatus = 2;

enum class Compression { none, zlib, lzma, bzip2, zstd };

std::ostream& operator<<(std::ostream& os, HepMC3Util::Compression compression);

std::span<const Compression> availableCompressionModes();

std::string_view compressionExtension(Compression compression);

enum class Format { ascii, root };

std::ostream& operator<<(std::ostream& os, Format format);

std::span<const Format> availableFormats();

Format formatFromFilename(std::string_view filename);

static constexpr std::string_view kEventGeneratorIndexAttribute =
    "acts_gen_event_index";

int eventGeneratorIndex(const HepMC3::GenParticle& particle);

}  // namespace ActsExamples::HepMC3Util
