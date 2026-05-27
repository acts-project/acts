/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/bfield/construct_const_bfield.hpp"

namespace traccc {

magnetic_field construct_const_bfield(const vector3& v) {

    return construct_const_bfield(v[0], v[1], v[2]);
}

}  // namespace traccc
