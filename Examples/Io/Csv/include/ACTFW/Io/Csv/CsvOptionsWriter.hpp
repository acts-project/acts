// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ACTFW/Io/Csv/CsvParticleWriter.hpp"
#include "ACTFW/Io/Csv/CsvPlanarClusterWriter.hpp"
#include "ACTFW/Io/Csv/CsvTrackingGeometryWriter.hpp"
#include "ACTFW/Utilities/OptionsFwd.hpp"

namespace FW {
namespace Options {

// Add common CSV writer options.
void addCsvWriterOptions(Description& desc);

/// Read the CSV particle writer options.
FW::CsvParticleWriter::Config readCsvParticleWriterConfig(const Variables& vm);

/// Read the CSV planar cluster writer options.
FW::CsvPlanarClusterWriter::Config readCsvPlanarClusterWriterConfig(
    const Variables& vm);

/// Read the CSV tracking geometry writer config.
FW::CsvTrackingGeometryWriter::Config readCsvTrackingGeometryWriterConfig(
    const Variables& vm);

}  // namespace Options
}  // namespace FW
