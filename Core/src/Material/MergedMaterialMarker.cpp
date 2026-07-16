// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/MergedMaterialMarker.hpp"

#include <ostream>

namespace Acts {

MergedMaterialMarker& MergedMaterialMarker::scale(double /*factor*/) {
  return *this;
}

const MaterialSlab& MergedMaterialMarker::materialSlab(
    const Vector2& /*lp*/) const {
  return m_slab;
}

std::vector<AxisDirection> MergedMaterialMarker::localAxisDirections() const {
  return {};
}

const MaterialSlab& MergedMaterialMarker::materialSlab(
    const Vector3& /*gp*/) const {
  return m_slab;
}

std::ostream& MergedMaterialMarker::toStream(std::ostream& sl) const {
  sl << "MergedMaterialMarker (material discarded during portal merge)";
  return sl;
}

}  // namespace Acts
