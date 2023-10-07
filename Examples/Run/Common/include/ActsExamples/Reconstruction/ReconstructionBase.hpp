// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/CsvOptionsReader.hpp"
#include "ActsExamples/Options/ParticleSmearingOptions.hpp"
#include "ActsExamples/TruthTracking/ParticleSmearing.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <filesystem>
#include <memory>
#include <string>

/// Setup sim hit csv reader
///
/// @param vars The configuration variables
/// @param sequencer The framework sequencer
///
/// @return config for sim hits csv reader
ActsExamples::CsvSimHitReader::Config setupSimHitReading(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer);

/// Setup sim particle csv reader
///
/// @param vars The configuration variables
/// @param sequencer The framework sequencer
///
/// @return config for sim particles csv reader
ActsExamples::CsvParticleReader::Config setupParticleReading(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer);

/// Setup sim hit smearing
///
/// @param vars The configuration variables
/// @param sequencer The framework sequencer
/// @param rnd The random number service
/// @param trackingGeometry The TrackingGeometry for the tracking setup
/// @param inputSimHits The input sim hit collection (e.g. from sim hit reader)
///
/// @return config for hit smearing
ActsExamples::DigitizationConfig setupDigitization(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const ActsExamples::RandomNumbers> rnd,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    const std::string& inputSimHits);

/// Setup particle smearing
///
/// @param vars The configuration variables
/// @param sequencer The framework sequencer
/// @param rnd The random number service
/// @param inputParticles The input particle collection (e.g. from particle
/// reader or from particle selection)
///
/// @return config for particle smearing
ActsExamples::ParticleSmearing::Config setupParticleSmearing(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer,
    std::shared_ptr<const ActsExamples::RandomNumbers> rnd,
    const std::string& inputParticles);

/// Setup reading measurements
///
/// @param vars The configuration variables
/// @param sequencer The framework sequencer
///
/// @return config for reading measurements
ActsExamples::CsvMeasurementReader::Config setupMeasurementReading(
    const ActsExamples::Options::Variables& vars,
    ActsExamples::Sequencer& sequencer);
