// This file is part of the Acts project.
//
// Copyright (C) 2020-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <Acts/EventData/Measurement.hpp>

#include <vector>

namespace ActsExamples {

/// Variable measurement type that can contain all possible combinations.
using Measurement = ::Acts::BoundVariantMeasurement;

/// Container of measurements.
///
/// In contrast to the source links, the measurements themself must not be
/// orderable. The source links stored in the measurements are treated
/// as opaque here and no ordering is enforced on the stored measurements.
using MeasurementContainer = std::vector<Measurement>;

}  // namespace ActsExamples
