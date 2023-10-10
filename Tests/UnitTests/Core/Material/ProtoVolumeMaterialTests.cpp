// This file is part of the Acts project.
//
// Copyright (C) 2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Material/ProtoVolumeMaterial.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <utility>

namespace Acts {

namespace Test {

/// Test the constructors
BOOST_AUTO_TEST_CASE(ProtoVolumeMaterial_construction_test) {
  BinUtility vmpBU(10, -10., 10., open, binX);
  vmpBU += BinUtility(10, -10., 10., open, binY);
  vmpBU += BinUtility(10, -10., 10., open, binZ);

  // Constructor from arguments
  ProtoVolumeMaterial vmp(vmpBU);
  // Copy constructor
  ProtoVolumeMaterial vmpCopy(vmp);
  // Copy move constructor
  ProtoVolumeMaterial vmpCopyMoved(std::move(vmpCopy));
}

}  // namespace Test
}  // namespace Acts
