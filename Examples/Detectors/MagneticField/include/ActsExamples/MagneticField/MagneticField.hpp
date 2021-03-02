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
#include "Acts/MagneticField/NullBField.hpp"
#include "ActsExamples/MagneticField/ScalableBField.hpp"

#include <memory>
#include <variant>

namespace ActsExamples {
namespace detail {

using InterpolatedMagneticFieldMapper2 = Acts::InterpolatedBFieldMapper<
    Acts::detail::Grid<Acts::Vector2, Acts::detail::EquidistantAxis,
                       Acts::detail::EquidistantAxis>>;
using InterpolatedMagneticField2 =
    Acts::InterpolatedBFieldMap<InterpolatedMagneticFieldMapper2>;

using InterpolatedMagneticFieldMapper3 =
    Acts::InterpolatedBFieldMapper<Acts::detail::Grid<
        Acts::Vector3, Acts::detail::EquidistantAxis,
        Acts::detail::EquidistantAxis, Acts::detail::EquidistantAxis>>;
using InterpolatedMagneticField3 =
    Acts::InterpolatedBFieldMap<InterpolatedMagneticFieldMapper3>;

}  // namespace detail

}  // namespace ActsExamples
