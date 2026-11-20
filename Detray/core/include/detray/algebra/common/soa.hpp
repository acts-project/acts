// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

/// Vector and matrix functions for the SoA plugins (Vc SoA and array SoA)
///
/// The functions operate on the @c algebra::storage types with a SIMD scalar.
///
/// @note The plugin has to declare the overloads of the math functions and
/// constants for its SIMD type before this header is included (see
/// @c algebra::soa::dot).

// Impl include(s).
#include "detray/algebra/common/impl/soa_matrix.hpp"
#include "detray/algebra/common/impl/soa_vector.hpp"
