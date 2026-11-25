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

// Detray test include(s)
#include "detray/test/common/material_scan_fixture.hpp"
#include "detray/test/common/material_validation_fixture.hpp"
#include "detray/test/validation/cpu/material_recorder.hpp"

namespace detray::test::host {

template <concepts::detector detector_t>
using material_scan = test::material_scan_fixture<detector_t>;

template <concepts::detector detector_t>
using material_validation = test::material_validation_fixture<
    detector_t, detray::validation::host::material_recorder>;

}  // namespace detray::test::host
