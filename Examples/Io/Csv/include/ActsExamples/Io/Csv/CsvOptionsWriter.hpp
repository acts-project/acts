// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Io/Csv/CsvMeasurementWriter.hpp"
#include "ActsExamples/Io/Csv/CsvParticleWriter.hpp"
#include "ActsExamples/Io/Csv/CsvPlanarClusterWriter.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitWriter.hpp"
#include "ActsExamples/Io/Csv/CsvTrackingGeometryWriter.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

namespace ActsExamples {
namespace Options {

// Add common CSV writer options.
void addCsvWriterOptions(Description& desc);

/// Read the CSV particle writer options.
ActsExamples::CsvParticleWriter::Config readCsvParticleWriterConfig(
    const Variables& vm);

/// Read the CSV sim hit writer options.
ActsExamples::CsvSimHitWriter::Config readCsvSimHitWriterConfig(
    const Variables& vm);

/// Read the CSV planar cluster writer options.
ActsExamples::CsvPlanarClusterWriter::Config readCsvPlanarClusterWriterConfig(
    const Variables& vm);

/// Read the CSV measurement writer options.
ActsExamples::CsvMeasurementWriter::Config readCsvMeasurementWriterConfig(
    const Variables& vm);

/// Read the CSV tracking geometry writer config.
ActsExamples::CsvTrackingGeometryWriter::Config
readCsvTrackingGeometryWriterConfig(const Variables& vm);

}  // namespace Options
}  // namespace ActsExamples
