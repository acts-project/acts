// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts {

// concept of measurement proxy:
// - access to the subset descriptor
// - access to the measurement data

// concept of measurement adapter backend:
// - access measurement range by surface

// concept of measurement adapter frontend:
// - all in one method - surface in, predicted parameters + covariance in, track
// states out

class MeasurementAdapterBase {
 public:
  // pre calibration
  // selection
  // post calibration
};

}  // namespace Acts
