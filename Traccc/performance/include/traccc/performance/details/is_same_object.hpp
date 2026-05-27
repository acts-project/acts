/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"

// System include(s).
#include <functional>

namespace traccc::details {

/// Generic unary functor for checking equality with a certain object
///
/// The idea is to use this class as-is for types that have a well-defined
/// equality operator, and to specialise it for cases where the uncertainty
/// of the equality can be set by the user. (To account for floating point
/// uncertainties in the algorithms.)
///
template <typename T>
class is_same_object {

    public:
    /// Constructor with a reference object, and an allowed uncertainty
    is_same_object(const T& ref, scalar /*unc*/ = float_epsilon);

    /// Default implementation relying on @c operator== for the given type
    bool operator()(const T& obj) const;

    private:
    /// The reference object
    std::reference_wrapper<const T> m_ref;

};  // class is_same_object

}  // namespace traccc::details

// Include the generic implementation.
#include "traccc/performance/impl/is_same_object.ipp"

// Include specialisations for the core library types
#include "traccc/performance/impl/is_same_measurement.ipp"
#include "traccc/performance/impl/is_same_seed.ipp"
#include "traccc/performance/impl/is_same_spacepoint.ipp"
#include "traccc/performance/impl/is_same_track.ipp"
#include "traccc/performance/impl/is_same_track_parameters.ipp"
#include "traccc/performance/impl/is_same_track_state.ipp"
