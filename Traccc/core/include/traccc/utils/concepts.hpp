/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include "detray/algebra/concepts.hpp"

namespace traccc::concepts {
using detray::concepts::column_matrix;
using detray::concepts::diagonal_matrix;
using detray::concepts::matrix;
using detray::concepts::row_matrix;
using detray::concepts::sized_matrix;
using detray::concepts::square_matrix;
using detray::concepts::symmetric_matrix;
}  // namespace traccc::concepts
