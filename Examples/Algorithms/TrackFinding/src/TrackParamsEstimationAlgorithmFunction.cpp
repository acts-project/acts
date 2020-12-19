// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
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
  // unpack the magnetic field variant and instantiate the corresponding fitter.
  return std::visit(
      [](auto&& inputField) -> BFieldGetter {
        // each entry in the variant is already a shared_ptr
        // need ::element_type to get the real magnetic field type
        using InputMagneticField =
            typename std::decay_t<decltype(inputField)>::element_type;
        using MagneticField = Acts::SharedBField<InputMagneticField>;

        MagneticField field(std::move(inputField));

        // build the B field getter functions. owns the fitter object.
        return BFieldGetterImpl<MagneticField>(std::move(field));
      },
      std::move(magneticField));
}
