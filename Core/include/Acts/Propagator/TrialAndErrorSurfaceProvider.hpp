// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

class TrialAndErrorSurfaceProvider {
  std::vector<const Surface *> m_all_surfaces;

 public:
  TrialAndErrorSurfaceProvider(const TrackingGeometry &tgeo) {
    tgeo.visitSurfaces([&](auto s) { m_all_surfaces.push_back(s); });
  }

  template <typename propagator_state_t, typename stepper_t>
  std::vector<const Surface *> operator()(const propagator_state_t &,
                                          const stepper_t &) const {
    return m_all_surfaces;
  }
};

}  // namespace Acts
