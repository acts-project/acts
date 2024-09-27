// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/Conversion/MeasurementConversion.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"

#include "traccc/definitions/track_parametrization.hpp"

namespace ActsExamples::Traccc::Common::Conversion {

Acts::BoundIndices boundIndex(const traccc::bound_indices tracccBoundIndex) {
  switch (tracccBoundIndex) {
    case traccc::bound_indices::e_bound_loc0:
      return Acts::BoundIndices::eBoundLoc0;
    case traccc::bound_indices::e_bound_loc1:
      return Acts::BoundIndices::eBoundLoc1;
    case traccc::bound_indices::e_bound_phi:
      return Acts::BoundIndices::eBoundPhi;
    case traccc::bound_indices::e_bound_theta:
      return Acts::BoundIndices::eBoundTheta;
    case traccc::bound_indices::e_bound_qoverp:
      return Acts::BoundIndices::eBoundQOverP;
    case traccc::bound_indices::e_bound_time:
      return Acts::BoundIndices::eBoundTime;
    case traccc::bound_indices::e_bound_size:
      return Acts::BoundIndices::eBoundSize;
    default:
      throw std::runtime_error("Could not convert traccc bound index");
  }
}

}  // namespace ActsExamples::Traccc::Common::Conversion
