// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Tests/CommonHelpers/TestSourceLink.hpp"

#include <ostream>

bool Acts::Test::operator==(const TestSourceLink& lhs,
                            const TestSourceLink& rhs) {
  return (lhs.m_geometryId == rhs.m_geometryId) &&
         (lhs.sourceId == rhs.sourceId) && (lhs.indices == rhs.indices) &&
         (lhs.parameters == rhs.parameters) &&
         (lhs.covariance == rhs.covariance);
}

bool Acts::Test::operator!=(const TestSourceLink& lhs,
                            const TestSourceLink& rhs) {
  return !(lhs == rhs);
}

std::ostream& Acts::Test::operator<<(
    std::ostream& os, const Acts::Test::TestSourceLink& sourceLink) {
  os << "TestsSourceLink(geometryId=" << sourceLink.m_geometryId
     << ",sourceId=" << sourceLink.sourceId;
  if (sourceLink.indices[0] != eBoundSize) {
    os << ",index0=" << sourceLink.indices[0];
  }
  if (sourceLink.indices[1] != eBoundSize) {
    os << ",index1=" << sourceLink.indices[1];
  }
  os << ")";
  return os;
}
