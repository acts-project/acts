/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/performance/details/is_same_object.hpp"

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"

namespace traccc::details {

/// Factory creating instances of "comparator objects" for a given type
///
/// This level of abstraction is necessary to be able to construct comparator
/// objects that would have extra configuration parameters over the reference
/// object and the comparison uncertainty.
///
/// @tparam TYPE The type for which a comparator object should be generated
///
template <typename TYPE>
class comparator_factory {

    public:
    /// Instantiate an instance of a comparator object
    is_same_object<TYPE> make_comparator(const TYPE& ref,
                                         scalar unc = float_epsilon) const;

};  // class comparator_factory

}  // namespace traccc::details

// Include the generic implementation.
#include "traccc/performance/impl/comparator_factory.ipp"

// Include the specialised implementation(s).
#include "traccc/performance/impl/seed_comparator_factory.ipp"
#include "traccc/performance/impl/track_comparator_factory.ipp"
