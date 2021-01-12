// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/Plugins/BField/ScalableBField.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"

namespace {

template <typename bfield_t>
struct BFieldGetterImpl {
  bfield_t bField;

  BFieldGetterImpl(bfield_t&& b) : bField(std::move(b)) {}

  Acts::Vector3 operator()(const Acts::Vector3& position) const {
    return bField.getField(position);
  };
};
}  // namespace

ActsExamples::TrackParamsEstimationAlgorithm::BFieldGetter
ActsExamples::TrackParamsEstimationAlgorithm::makeBFieldGetter(
    Options::BFieldVariant magneticField) {
  return std::visit(
      [](auto&& inputField) -> BFieldGetter {
        using InputMagneticField =
            typename std::decay_t<decltype(inputField)>::element_type;
        using MagneticField = Acts::SharedBField<InputMagneticField>;
        MagneticField field(std::move(inputField));

        return BFieldGetterImpl<MagneticField>(std::move(field));
      },
      std::move(magneticField));
}
