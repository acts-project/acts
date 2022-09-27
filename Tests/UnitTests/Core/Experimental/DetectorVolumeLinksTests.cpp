// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Experimental/NavigationState.hpp"
#include "Acts/Experimental/detail/DetectorVolumeLinks.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <memory>

namespace Acts {
namespace Experimental {
class DetectorVolume {};
}  // namespace Experimental
}  // namespace Acts

using namespace Acts::Experimental;

Acts::GeometryContext tContext;

BOOST_AUTO_TEST_SUITE(Experimental)

BOOST_AUTO_TEST_CASE(SingleVolumeLink) {
  auto dVolume = std::make_unique<DetectorVolume>();
  detail::SingleLinkImpl singleLink(*dVolume);

  // Plain test
  BOOST_CHECK(dVolume.get() == singleLink.dVolume);
  BOOST_CHECK(dVolume.get() ==
              singleLink.targetVolume(tContext, Acts::Vector3(0., 0., 0.),
                                      Acts::Vector3(0., 0., 0.)));

  // Test as a delegate
  DetectorVolumeLink singleLinkDel;
  singleLinkDel.connect<&detail::SingleLinkImpl::targetVolume>(&singleLink);

  BOOST_CHECK(dVolume.get() == singleLinkDel(tContext, Acts::Vector3(0., 0., 0.),
                                      Acts::Vector3(0., 0., 0.)));

  // Copy assigned of the delegate 
  DetectorVolumeLink singleLinkDelCopy = singleLinkDel;
  BOOST_CHECK(dVolume.get() == singleLinkDelCopy(tContext, Acts::Vector3(0., 0., 0.),
                                      Acts::Vector3(0., 0., 0.)));
}

BOOST_AUTO_TEST_CASE(MultiVolumeLink1DCast) {
  // Create a, b, c volumes
  auto aVolume = std::make_unique<DetectorVolume>();
  auto bVolume = std::make_unique<DetectorVolume>();
  auto cVolume = std::make_unique<DetectorVolume>();

  std::vector<const DetectorVolume*> dVolumes = {aVolume.get(), bVolume.get(),
                                                 cVolume.get()};
  std::vector<Acts::ActsScalar> zValues = {-10., 0., 100., 200.};
  detail::MultiLink1DImpl zLinks(dVolumes, zValues, Acts::binZ);

  // Check if you get the right volumes back
  Acts::Vector3 getA(1., 2., -5.);
  Acts::Vector3 getB(10., 10., 50.);
  Acts::Vector3 getC(0., 0., 170.);

  BOOST_CHECK(aVolume.get() ==
              zLinks.targetVolume(tContext, getA, Acts::Vector3(0., 0., 0.)));

  BOOST_CHECK(bVolume.get() ==
              zLinks.targetVolume(tContext, getB, Acts::Vector3(0., 0., 0.)));

  BOOST_CHECK(cVolume.get() ==
              zLinks.targetVolume(tContext, getC, Acts::Vector3(0., 0., 0.)));
}

BOOST_AUTO_TEST_SUITE_END()
