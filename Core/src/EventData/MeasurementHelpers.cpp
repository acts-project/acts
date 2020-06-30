// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/MeasurementHelpers.hpp"

bool Acts::MinimalSourceLink::operator==(
    const Acts::MinimalSourceLink& rhs) const {
  return meas == rhs.meas;
}

const Acts::Surface& Acts::MinimalSourceLink::referenceSurface() const {
  return *MeasurementHelpers::getSurface(*meas);
}

const Acts::FittableMeasurement<Acts::MinimalSourceLink>&
    Acts::MinimalSourceLink::operator*() const {
  return *meas;
}