// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#define BOOST_TEST_MODULE Layer Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include <boost/test/data/test_case.hpp>
// leave blank line

#include <boost/test/output_test_stream.hpp>
// leave blank line

//#include <limits>
#include "ACTS/Layers/Layer.hpp"
//#include "ACTS/Utilities/Definitions.hpp"
#include "LayerStub.hpp"
#include "ACTS/Layers/GenericApproachDescriptor.hpp"

using boost::test_tools::output_test_stream;
namespace utf    = boost::unit_test;

namespace Acts {

namespace Test {
  BOOST_AUTO_TEST_SUITE(Layers);

  /// todo: make test fixture; separate out different cases

  /// Unit test for creating compliant/non-compliant Layer object
  BOOST_AUTO_TEST_CASE(LayerConstruction)
  {
    ///Minimum possible construction (default constructor is deleted)
    LayerStub minallyConstructed(nullptr);
    BOOST_TEST(minallyConstructed.layerType() == passive);
    ///Need an approach descriptor for the next level of complexity:
    //const std::vector<const Surface*> aSurfaces{new SurfaceStub(), new SurfaceStub()};
    //std::unique_ptr<ApproachDescriptor> ad(new GenericApproachDescriptor<Surface>(aSurfaces));
    //LayerStub approachDescriptorConstructed(nullptr, 1.0, ad);
    

  }

  /// Unit test for testing Layer properties
  BOOST_AUTO_TEST_CASE(LayerProperties)
  {
   
  }

  BOOST_AUTO_TEST_SUITE_END();

}  // end of namespace Test

}  // end of namespace Acts
