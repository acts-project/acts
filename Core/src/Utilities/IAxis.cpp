// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/IAxis.hpp"

namespace Acts {

bool operator==(const IAxis& lhs, const IAxis& rhs) {
  return lhs.getType() == rhs.getType() &&
         lhs.getBoundaryType() == rhs.getBoundaryType() &&
         lhs.getMin() == rhs.getMin() && lhs.getMax() == rhs.getMax() &&
         lhs.getNBins() == rhs.getNBins() &&
         lhs.getBinEdges() == rhs.getBinEdges();
}

std::ostream& operator<<(std::ostream& os, const IAxis& axis) {
  axis.toStream(os);
  return os;
}

}  // namespace Acts
