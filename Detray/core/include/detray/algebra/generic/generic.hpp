// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Impl include(s).
#include "detray/algebra/common/boolean.hpp"
#include "detray/algebra/common/math.hpp"
#include "detray/algebra/generic/impl/generic_matrix.hpp"
#include "detray/algebra/generic/impl/generic_transform3.hpp"
#include "detray/algebra/generic/impl/generic_vector.hpp"

// Algorithms include(s).
#include "detray/algebra/generic/algorithms/matrix/decomposition/partial_pivot_lud.hpp"
#include "detray/algebra/generic/algorithms/matrix/determinant/cofactor.hpp"
#include "detray/algebra/generic/algorithms/matrix/determinant/hard_coded.hpp"
#include "detray/algebra/generic/algorithms/matrix/determinant/partial_pivot_lud.hpp"
#include "detray/algebra/generic/algorithms/matrix/inverse/cofactor.hpp"
#include "detray/algebra/generic/algorithms/matrix/inverse/hard_coded.hpp"
#include "detray/algebra/generic/algorithms/matrix/inverse/partial_pivot_lud.hpp"
#include "detray/algebra/generic/algorithms/utils/algorithm_selector.hpp"
