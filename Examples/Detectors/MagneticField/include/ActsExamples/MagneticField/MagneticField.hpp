// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/AxisFwd.hpp"
#include "Acts/Utilities/detail/grid_helper.hpp"
#include "ActsExamples/MagneticField/ScalableBField.hpp"

#include <memory>
#include <variant>
#include <vector>

namespace ActsExamples::detail {

using InterpolatedMagneticField2 = Acts::InterpolatedBFieldMap<
    Acts::Grid<Acts::Vector2, Acts::detail::EquidistantAxis,
               Acts::detail::EquidistantAxis>>;

using InterpolatedMagneticField3 = Acts::InterpolatedBFieldMap<
    Acts::Grid<Acts::Vector3, Acts::detail::EquidistantAxis,
               Acts::detail::EquidistantAxis, Acts::detail::EquidistantAxis>>;

}  // namespace ActsExamples::detail
