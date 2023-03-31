// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

namespace Acts {
namespace Experimental {

template <typename axis_generator>
struct CylindricalSurfaceGridGenerator {
  // The surfaces to be indexed
  std::vector<std::shared_ptr<Surface>> surfaces = {};
};

/*
inline static SurfaceCandidatesUpdator tryAllPortalsAndSurfaces() {
  auto aps = std::make_unique<const AllPortalsAndSurfacesImpl>();
  SurfaceCandidatesUpdator nStateUpdator;
  nStateUpdator.connect<&AllPortalsAndSurfacesImpl::update>(std::move(aps));
  return nStateUpdator;
}
*/

}  // namespace Experimental

}  // namespace Acts
