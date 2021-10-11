// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Experimental/PortalLinks.hpp"
#include "Acts/Experimental/VolumeLinks.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/detail/Axis.hpp"

#include <functional>

namespace Acts {
namespace Test {

GeometryContext gctx = GeometryContext();

BOOST_AUTO_TEST_SUITE(Experimental)

BOOST_AUTO_TEST_CASE(VoidVolumeLink_) {
  // The transform
  Transform3 transform = Transform3::Identity();
  Vector3 position = Vector3(0., 0., 0.);

  // Default volume link
  VolumeLink vLink = VoidVolumeLink{};
  unsigned int vIndex = vLink(transform, position);
  // Test the default link
  BOOST_CHECK(vIndex == 0u);
}

BOOST_AUTO_TEST_CASE(EquidistantVolumeLink_) {
  // The transform
  Transform3 transform = Transform3::Identity();
  transform.pretranslate(Vector3(0., 0., 5.));

  // -3 |> bin0 <| 0 |> bin1 <| 3 |> bin2 <| 6 //
  detail::EquidistantAxis eAxis(-3., 6., 3u);
  EquidistantVolumeLink evLink(std::move(eAxis));
  evLink.bvalue = Acts::binZ;

  // Assign to the volume link
  VolumeLink vLink = evLink;
  // Test underflow (0,0,0) -> (0,0,-5) : bin 0
  unsigned int vIndex = vLink(transform, Vector3(0., 0., 0.));
  BOOST_CHECK(vIndex == 0u);
  // Test (0,0,3) -> (0,0,-2) : bin 0
  vIndex = vLink(transform, Vector3(0., 0., 3.));
  BOOST_CHECK(vIndex == 0u);
  // Test (0,0,6) -> (0,0,1): bin 1
  vIndex = vLink(transform, Vector3(0., 0., 6.));
  BOOST_CHECK(vIndex == 1u);
  // Test (0,0,8.5) -> (0,0,3.5) : bin 2
  vIndex = vLink(transform, Vector3(0., 0., 8.5));
  BOOST_CHECK(vIndex == 2u);
  // Test the overflow (0,0,20) -> (0,0,15) : bin 2
  vIndex = vLink(transform, Vector3(0., 0., 20.));
  BOOST_CHECK(vIndex == 2u);
}

BOOST_AUTO_TEST_CASE(VariableVolumeLink_) {
  // The transform
  Transform3 transform = Transform3::Identity();
  transform.pretranslate(Vector3(0., 0., 5.));

  // -3.5 |> bin0 <| 2. |> bin1 <| 6.5 |> bin2 <| 8.0 //
  detail::VariableAxis vAxis({-3.5, 2., 6.5, 8.0});
  VariableVolumeLink vvLink(std::move(vAxis));
  vvLink.bvalue = Acts::binZ;

  // Assign to the volume link
  VolumeLink vLink = vvLink;
  // Test the underflow (0,0,0) -> (0,0,-5)
  unsigned int vIndex = vLink(transform, Vector3(0., 0., 0.));
  BOOST_CHECK(vIndex == 0u);
  // Test (0,0,5) -> (0,0,0) : bin 0
  vIndex = vLink(transform, Vector3(0., 0., 5.));
  BOOST_CHECK(vIndex == 0u);
  // Test (0,0,8) -> (0,0,3): bin 1
  vIndex = vLink(transform, Vector3(0., 0., 8.));
  BOOST_CHECK(vIndex == 1u);
  // Test (0,0,12) -> (0,0,7) : bin 2
  vIndex = vLink(transform, Vector3(0., 0., 12.));
  BOOST_CHECK(vIndex == 2u);
  // Test the overflow (0,0,20) -> (0,0,15) : bin 2
  vIndex = vLink(transform, Vector3(0., 0., 20.));
  BOOST_CHECK(vIndex == 2u);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Test
}  // namespace Acts