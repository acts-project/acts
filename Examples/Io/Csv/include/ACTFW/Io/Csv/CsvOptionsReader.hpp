// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ACTFW/Io/Csv/CsvParticleReader.hpp"
#include "ACTFW/Io/Csv/CsvPlanarClusterReader.hpp"
#include "ACTFW/Utilities/OptionsFwd.hpp"

namespace FW {
namespace Options {

// There are no additional CSV reader options apart from the
// format-independent, generic input option.

/// Read the CSV particle reader config.
FW::CsvParticleReader::Config readCsvParticleReaderConfig(const Variables& vm);

/// Read the CSV particle reader config.
FW::CsvPlanarClusterReader::Config readCsvPlanarClusterReaderConfig(
    const Variables& vm);

}  // namespace Options
}  // namespace FW
