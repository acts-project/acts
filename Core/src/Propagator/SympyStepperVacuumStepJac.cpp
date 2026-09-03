// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// One specialisation per translation unit on purpose; see the header.

#include "Acts/Propagator/detail/SympyStepperVacuumStep.hpp"
#include "Acts/Propagator/detail/SympyStepperVacuumStepDecl.hpp"

namespace Acts::detail {

template Result<double> sympyVacuumStep<true>(const SympyStepper&,
                                              SympyStepper::State&, Direction);

}  // namespace Acts::detail
