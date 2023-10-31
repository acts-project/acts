// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/VolumeBounds.hpp"

const std::vector<std::string> Acts::VolumeBounds::s_boundsTypeNames = {
    "Cone",          "Cuboid",    "CutoutCylinder", "Cylinder",
    "GenericCuboid", "Trapezoid", "Other"};

/**Overload of << operator for std::ostream for debug output*/
std::ostream& Acts::operator<<(std::ostream& sl, const Acts::VolumeBounds& vb) {
  return vb.toStream(sl);
}
