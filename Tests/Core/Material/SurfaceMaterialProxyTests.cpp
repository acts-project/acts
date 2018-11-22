// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE SurfaceMaterialProxy Tests
#include <boost/test/included/unit_test.hpp>
#include <climits>
#include "Acts/Material/SurfaceMaterialProxy.hpp"

namespace Acts {

namespace Test {

  /// Test the constructors
  BOOST_AUTO_TEST_CASE(SurfaceMaterialProxy_construction_test)
  {
    BinUtility smpBU(10, -10., 10., open, binX);
    smpBU += BinUtility(10, -10., 10., open, binY);

    // Constructor from arguments
    SurfaceMaterialProxy smp(smpBU);
    // Copy constructor
    SurfaceMaterialProxy smpCopy(smp);
    // Copy move constructor
    SurfaceMaterialProxy smpCopyMoved(std::move(smpCopy));
  }

}  // namespace Test
}  // namespace Acts