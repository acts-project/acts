// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Digitization/HitSmearing.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/TruthTracking/ParticleSmearing.hpp"
#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <memory>
#include <string>

#include <boost/filesystem.hpp>

/// Setup sim hit csv reader
///
/// @param variables The configuration variables
/// @param sequencer The framework sequencer
///
/// @return config for sim hits csv reader
ActsExamples::CsvSimHitReader::Config setupSimHitReading(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer);

/// Setup sim particle csv reader
///
/// @param variables The configuration variables
/// @param sequencer The framework sequencer
///
/// @return config for sim particles csv reader
ActsExamples::CsvParticleReader::Config setupParticleReading(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer);

/// Run sim hit smearing
///
/// @param variables The configuration variables
/// @param sequencer The framework sequencer
/// @param randomNumbers The random number service
/// @param trackingGeometry The TrackingGeometry for the tracking setup
/// @param inputSimHits The input sim hit collection (e.g. from sim hit reader)
///
/// @return config for hit smearing
ActsExamples::HitSmearing::Config runSimHitSmearing(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const ActsExamples::RandomNumbers> rnd,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    const std::string& inputSimHits);

/// Run particle smearing
///
/// @param variables The configuration variables
/// @param sequencer The framework sequencer
/// @param randomNumbers The random number service
/// @param inputParticles The input particle collection (e.g. from particle
/// reader or from particle selection)
/// @param inputMeasurementParticlesMap The input map between measurement and
/// particle (e.g. from hit smearing)
///
/// @return config for particle smearing
ActsExamples::ParticleSmearing::Config runParticleSmearing(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const ActsExamples::RandomNumbers> rnd,
    const std::string& inputParticles);

/// Run particle selection
///
/// @param variables The configuration variables
/// @param sequencer The framework sequencer
/// @param inputParticles The input particle collection (e.g. from particle
/// reader)
/// @param inputMeasurementParticlesMap The input map between measurement and
/// particle (e.g. from hit smearing)
///
/// @return config used for particle selection
ActsExamples::TruthSeedSelector::Config runParticleSelection(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer, const std::string& inputParticles,
    const std::string& inputMeasurementParticlesMap);
