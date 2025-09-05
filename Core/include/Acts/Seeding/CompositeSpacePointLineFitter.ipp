// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/detail/CompositeSpacePointLineFitter.hpp"
namespace Acts::Experimental {

template <CompositeSpacePointContainer Cont_t,
            CompositeSpacePointCalibrator<Cont_t, Cont_t> Calibrator_t>
  CompositeSpacePointLineFitter::FitResult<Cont_t> CompositeSpacePointLineFitter::fit(FitOptions<Cont_t, Calibrator_t>&& fitOpts) const {
      FitResult<Cont_t> result{};
      return result;
  }




}  // namespace Acts::Experimental

