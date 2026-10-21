// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// One instantiation per translation unit; see the header.

#include "SympyStepperDenseStepImpl.hpp"
#include "SympyStepperStepImpl.hpp"

namespace Acts {

template detail::Rk4Status detail::sympyDenseStep<true>(
    const SympyStepper&, SympyStepper::State&, const IVolumeMaterial&, double,
    double, double&, Vector3&, std::error_code&, std::span<double>);

template Result<double> detail::sympyStep<detail::SympyStepMode::Dense, true>(
    const SympyStepper&, SympyStepper::State&, Direction,
    const IVolumeMaterial*);

}  // namespace Acts
