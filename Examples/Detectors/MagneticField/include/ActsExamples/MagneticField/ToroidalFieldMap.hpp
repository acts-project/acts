// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/Utilities/Grid.hpp"

#include <tuple>
#include <utility>

namespace Acts {

// Forward declaration
class ToroidalField;

// 3D grid in cylindrical (r,phi,z), table stores (Br,Bphi,Bz)
InterpolatedBFieldMap<
    Grid<Vector3, Axis<AxisType::Equidistant>, Axis<AxisType::Equidistant>,
         Axis<AxisType::Equidistant>>>
toroidalFieldMapCyl(
    const std::pair<double, double>& rLim,
    const std::pair<double, double>& phiLim,
    const std::pair<double, double>& zLim,
    const std::tuple<std::size_t, std::size_t, std::size_t>& nBins,
    const ToroidalField& field);

// Optional: a Cartesian (x,y,z) variant that stores (Bx,By,Bz) directly.
InterpolatedBFieldMap<
    Grid<Vector3, Axis<AxisType::Equidistant>, Axis<AxisType::Equidistant>,
         Axis<AxisType::Equidistant>>>
toroidalFieldMapXYZ(
    const std::pair<double, double>& xLim,
    const std::pair<double, double>& yLim,
    const std::pair<double, double>& zLim,
    const std::tuple<std::size_t, std::size_t, std::size_t>& nBins,
    const ToroidalField& field);

}  // namespace Acts
