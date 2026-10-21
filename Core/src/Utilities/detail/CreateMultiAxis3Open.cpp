// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Explicit instantiations of the 3D multi-axis factory for the
// AxisBoundaryType::Open first axis. Splitting the instantiations across one
// TU per boundary type bounds the per-TU compiler memory of the otherwise
// 216-way combinatorial `MultiAxis` instantiation.

#include "CreateMultiAxis3Impl.hpp"

namespace Acts::detail {

template std::unique_ptr<IMultiAxis3D>
createMultiAxis3<Axis<AxisType::Equidistant, AxisBoundaryType::Open>>(
    const Axis<AxisType::Equidistant, AxisBoundaryType::Open>& first,
    const IAxis& second, const IAxis& third);

template std::unique_ptr<IMultiAxis3D>
createMultiAxis3<Axis<AxisType::Variable, AxisBoundaryType::Open>>(
    const Axis<AxisType::Variable, AxisBoundaryType::Open>& first,
    const IAxis& second, const IAxis& third);

}  // namespace Acts::detail
