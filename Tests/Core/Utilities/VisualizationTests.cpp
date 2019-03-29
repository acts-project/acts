// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE VisualizationTest Test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/IVisualization.hpp"
#include "Acts/Utilities/PlyHelper.hpp"
#include "Acts/Utilities/ObjHelper.hpp"

namespace Acts {
namespace Test {

  BOOST_AUTO_TEST_SUITE(Utilities)

  BOOST_AUTO_TEST_CASE(construction_test)
  {
    // this doesn't really test anything, other than conformance to the IVisualization interface
    PlyHelper ply;
    ObjHelper obj;

    IVisualization* vis;
    vis = &ply;
    vis = &obj;

    (void)vis;
  }

  BOOST_AUTO_TEST_SUITE_END()
}}
