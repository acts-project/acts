// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Tests/CommonHelpers/TestSourceLink.hpp"

#include <ostream>

std::ostream& Acts::Test::operator<<(
    std::ostream& os, const Acts::Test::TestSourceLink& sourceLink) {
  os << "TestsSourceLink(geometryId=" << sourceLink.geometryId()
     << ",measurement=" << sourceLink.measurement << ")";
  return os;
}
