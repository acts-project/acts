// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <utility>

namespace Acts::Test {

/// Test the constructors
BOOST_AUTO_TEST_CASE(ProtoSurfaceMaterial_construction_test) {
  BinUtility smpBU(10, -10., 10., open, binX);
  smpBU += BinUtility(10, -10., 10., open, binY);

  // Constructor from arguments
  ProtoSurfaceMaterial smp(smpBU);
  // Copy constructor
  ProtoSurfaceMaterial smpCopy(smp);
  // Copy move constructor
  ProtoSurfaceMaterial smpCopyMoved(std::move(smpCopy));
}

}  // namespace Acts::Test
