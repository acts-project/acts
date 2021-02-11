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

/// Magnetic field variant with all supported fields.
///
/// This is a value-like type, i.e. can be copied and stored by-value, that can
/// be used wherever magnetic field information is needed. The examples support
/// only the finite set of magnetic field type contained in this variant. This
/// enables the use of a single concrete value-like type that can be used in
/// interfaces.
using MagneticField =
    std::variant<std::shared_ptr<Acts::NullBField>,
                 std::shared_ptr<Acts::ConstantBField>,
                 std::shared_ptr<ScalableBField>,
                 std::shared_ptr<detail::InterpolatedMagneticField2>,
                 std::shared_ptr<detail::InterpolatedMagneticField3>>;

}  // namespace ActsExamples
