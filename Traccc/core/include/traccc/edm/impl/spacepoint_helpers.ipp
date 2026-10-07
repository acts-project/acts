// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/utils/detray_conversion.hpp"

namespace traccc::edm {

template <detray::concepts::algebra algebra_t, typename spacepoint_backend_t>
TRACCC_HOST_DEVICE detray::dpoint3D<algebra_t> get_spacepoint_global(
    const edm::spacepoint<spacepoint_backend_t>& sp) {
  return utils::to_dpoint3D<algebra_t>(sp.global());
}

}  // namespace traccc::edm
