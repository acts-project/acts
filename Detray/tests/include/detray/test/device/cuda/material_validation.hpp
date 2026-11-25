// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/core/concepts.hpp"
#include "detray/core/detector.hpp"

// Detray test include(s)
#include "detray/test/common/material_validation_fixture.hpp"
#include "detray/test/validation/cuda/material_recorder.hpp"

namespace detray::test::cuda {

// The material scan is only available on the host, so the host material
// validation needs to be included for this

/// Material validator using the CUDA material recording
template <concepts::detector detector_t>
using material_validation = detray::test::material_validation_fixture<
    detector_t, detray::validation::cuda::material_recorder>;

}  // namespace detray::test::cuda
