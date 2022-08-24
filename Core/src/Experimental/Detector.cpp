// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Experimental/Detector.hpp"

#include "Acts/Experimental/NavigationState.hpp"
#include "Acts/Experimental/Portal.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Helpers.hpp"

Acts::Experimental::Detector::Detector(
    const std::vector<std::shared_ptr<DetectorVolume>>& volumes,
    const DetectorVolumeFinder& volumeFinder, const std::string& name) {}

std::shared_ptr<Acts::Experimental::Detector>
Acts::Experimental::Detector::getSharedPtr() {
  return shared_from_this();
}

std::shared_ptr<const Acts::Experimental::Detector>
Acts::Experimental::Detector::getSharedPtr() const {
  return shared_from_this();
}
