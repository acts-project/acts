// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Vector, matrix and transform types on top of a storage array
//
// The types are used by the plugins that are based on the
// algebra::storage types (Vc AoS, Vc SoA and array SoA). Their scalar type
// can be a plain value type (AoS) or a SIMD type (SoA).
//
// Note: The overloads of approx_equal for the scalar type of a plugin
// have to be declared before this header is included (see
// algebra::storage::transform3).

// Impl include(s).
#include "detray/algebra/common/impl/storage_array_operators.hpp"
#include "detray/algebra/common/impl/storage_getter.hpp"
#include "detray/algebra/common/impl/storage_matrix.hpp"
#include "detray/algebra/common/impl/storage_transform3.hpp"
#include "detray/algebra/common/impl/storage_vector.hpp"
