// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <boost/test/data/test_case.hpp>
#include <cmath>
#include <vector>

#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Utilities/Units.hpp"

namespace {

namespace bdata = boost::unit_test::data;
using namespace Acts;

// reference surfaces
// this includes only those surfaces that can take unbounded local positions as
// inputs, i.e. no angles or strictly positive radii.
const auto surfaces = bdata::make(std::vector<std::shared_ptr<const Surface>>{
    Surface::makeShared<CylinderSurface>(
        std::make_shared<Transform3D>(Transform3D::Identity()), 10 /* radius */,
        100 /* half-length z */),
    // TODO perigee roundtrip local->global->local does not seem to work
    // Surface::makeShared<PerigeeSurface>(Vector3D(0, 0, -1.5)),
    Surface::makeShared<PlaneSurface>(Vector3D::Zero(), Vector3D::UnitX()),
    Surface::makeShared<PlaneSurface>(Vector3D::Zero(), Vector3D::UnitY()),
    Surface::makeShared<PlaneSurface>(Vector3D::Zero(), Vector3D::UnitZ()),
});
// positions
const auto posAngle = bdata::xrange(-M_PI, M_PI, 0.25);
const auto posPositiveNonzero = bdata::xrange(0.25, 1.0, 0.25);
const auto posPositive = bdata::make(0.0) + posPositiveNonzero;
const auto posSymmetric = bdata::xrange(-1.0, 1.0, 0.25);
// time
const auto ts = bdata::make(1.0);
// direction angles
const auto phis = bdata::make({0.0, M_PI, -M_PI, M_PI_2, -M_PI_2});
const auto thetasNoForwardBackward = bdata::xrange(M_PI_4, M_PI, M_PI_4);
const auto thetas = bdata::make({0.0, M_PI}) + thetasNoForwardBackward;
// p,q,q/p are intened to be zipped not to construct a cartesian product
const auto ps = bdata::make(1.0);
const auto qs = bdata::make(-UnitConstants::e);
const auto qOverPs = bdata::make(-UnitConstants::e / 1.0);

}  // namespace
