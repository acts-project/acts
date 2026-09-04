# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Set up the core tests.
detray_add_unit_test( svgtools
   "detectors.cpp"
   "grids.cpp"
   "groups.cpp"
   "intersections.cpp"
   "landmarks.cpp"
   "masks.cpp"
   "material.cpp"
   "surfaces.cpp"
   "trajectories.cpp"
   "volumes.cpp"
   LINK_LIBRARIES GTest::gtest_main detray::core_array detray::test_framework detray::test_common detray::validation_utils
                  actsvg::core actsvg::meta detray::svgtools
)

detray_add_unit_test( svgtools_web
   "web.cpp"
   LINK_LIBRARIES GTest::gtest_main detray::core_array detray::test_framework detray::test_common detray::validation_utils
                  actsvg::core actsvg::meta detray::svgtools actsvg::web
)
