// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

#include <iostream>
#include <vector>

using namespace std;
namespace Acts {
namespace Test {
    // make detector

    // make tracking geometry
    GeometryContext geoContext = GeometryContext();
BOOST_AUTO_TEST_SUITE(SurfaceContainer)

BOOST_AUTO_TEST_CASE(SurfaceContainerTest) {
        // Make surfaces
        int numSurfaces = 6;
        auto rBounds = std::make_shared<const RectangleBounds>(3., 4.);
        std::vector<std::shared_ptr<Acts::Surface>> surfaces;
        std::vector<Transform3> transformations;
            for (int i=0 ; i<numSurfaces ; i++) {
                Translation3 translation{i, 1., 2.};
                auto pTransform = Transform3(translation);
                transformations.push_back(pTransform);
                auto planeSurface =
                        Surface::makeShared<PlaneSurface>(pTransform, rBounds);
                surfaces.push_back(planeSurface);
                BOOST_CHECK_EQUAL(
                        surfaces.at(i)->transform(geoContext).isApprox(transformations.at(i)), true
                );
            }



        // make points

        // make spacepoints

        using namespace UnitLiterals;


}

} // namespace Test
} // namcespace Acts


}
