// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/detail/AxisFwd.hpp"
#include "Acts/Utilities/detail/GridFwd.hpp"
#include "ActsExamples/Utilities/OptionsFwd.hpp"

#include <memory>
#include <tuple>
#include <variant>

// Forward declarations
namespace Acts {
template <typename G>
struct InterpolatedBFieldMapper;

template <typename M>
class InterpolatedBFieldMap;

class ConstantBField;
}  // namespace Acts

namespace ActsExamples {
namespace BField {
class ScalableBField;
}
}  // namespace ActsExamples

using InterpolatedMapper2D = Acts::InterpolatedBFieldMapper<
    Acts::detail::Grid<Acts::Vector2D, Acts::detail::EquidistantAxis,
                       Acts::detail::EquidistantAxis>>;

using InterpolatedMapper3D = Acts::InterpolatedBFieldMapper<Acts::detail::Grid<
    Acts::Vector3D, Acts::detail::EquidistantAxis,
    Acts::detail::EquidistantAxis, Acts::detail::EquidistantAxis>>;

using InterpolatedBFieldMap2D =
    Acts::InterpolatedBFieldMap<InterpolatedMapper2D>;
using InterpolatedBFieldMap3D =
    Acts::InterpolatedBFieldMap<InterpolatedMapper3D>;

namespace ActsExamples {

namespace Options {

using BFieldVariant =
    std::variant<std::shared_ptr<InterpolatedBFieldMap2D>,
                 std::shared_ptr<InterpolatedBFieldMap3D>,
                 std::shared_ptr<Acts::ConstantBField>,
                 std::shared_ptr<ActsExamples::BField::ScalableBField>>;

// common bfield options, with a bf prefix
void addBFieldOptions(boost::program_options::options_description& opt);

// create the bfield maps
BFieldVariant readBField(const boost::program_options::variables_map& vm);

}  // namespace Options
}  // namespace ActsExamples
