// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// One specialisation per translation unit; see the vacuum step's header.

#include "Acts/Propagator/detail/SympyStepperDenseStep.hpp"
#include "Acts/Propagator/detail/SympyStepperDenseStepImpl.hpp"

namespace Acts::detail {

template Result<double> sympyDenseStepFull<false>(const SympyStepper&,
                                                  SympyStepper::State&,
                                                  Direction,
                                                  const IVolumeMaterial*);

}  // namespace Acts::detail
