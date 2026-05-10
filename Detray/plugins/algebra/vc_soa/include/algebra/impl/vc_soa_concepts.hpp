// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "algebra/impl/vc_aos_concepts.hpp"

// System include(s).
#include <concepts>

namespace detray::concepts {

/// Vc SoA vector
template <typename T>
concept vc_soa_vector =
    (simd_storage_vector<T> || vc_simd_vector<typename T::value_type>);

}  // namespace detray::concepts
