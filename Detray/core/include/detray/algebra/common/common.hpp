// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Common implementations that are based on the common storage types
//
// The functions operate on the algebra::storage types with a SIMD scalar.
// The vector functions use the generic implementation (see
// detray/algebra/generic/impl/generic_vector.hpp).

// Impl include(s).
#include "detray/algebra/common/impl/matrix.hpp"
#include "detray/algebra/common/impl/transform3.hpp"
#include "detray/algebra/common/storage/matrix.hpp"
#include "detray/algebra/common/storage/vector.hpp"
