// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Experimental/detail/PortalGenerators.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <memory>

namespace Acts {
namespace Experimental {
class DetectorVolume {};
}  // namespace Experimental
}  // namespace Acts

/// Unpack to shared - simply to test the getSharedPtr mechanism
///
/// @tparam referenced_type is the type of the referenced object
///
/// @param rt is the referenced object
///
/// @returns a shared pointer
template <typename referenced_type>
std::shared_ptr<referenced_type> unpackToShared(referenced_type& rt) {
  return rt.getSharedPtr();
}

using namespace Acts::Experimental;

// A test context
Acts::GeometryContext tContext;

BOOST_AUTO_TEST_SUITE(Experimental)

DetectorVolume dVolume;
auto dTransform = Acts::Transform3::Identity();
auto pGenerator = detail::defaultPortalGenerator();

BOOST_AUTO_TEST_CASE(PortalTest) {
  // A rectangle
  auto rectangle = std::make_shared<Acts::RectangleBounds>(10., 100.);
  auto surface =
      Acts::Surface::makeShared<Acts::PlaneSurface>(dTransform, rectangle);

  // Create a portal out of it
  auto portal = Portal::makeShared(surface);

  BOOST_TEST(&(portal->surface()), surface.get());

  BOOST_CHECK(portal == unpackToShared<Portal>(*portal));
  BOOST_CHECK(portal == unpackToShared<const Portal>(*portal));
}

BOOST_AUTO_TEST_SUITE_END()
