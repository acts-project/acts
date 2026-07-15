/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/performance/details/is_same_scalar.hpp"

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"

namespace traccc::details {

/// @c traccc::is_same_object specialisation for @c traccc::measurement
template <typename T>
class is_same_object<edm::measurement<T>> {

    public:
    /// Constructor with a reference object, and an allowed uncertainty
    is_same_object(const edm::measurement<T>& ref, scalar unc = float_epsilon)
        : m_ref(ref), m_unc(unc) {}

    /// Specialised implementation for @c traccc::measurement
    bool operator()(const edm::measurement<T>& obj) const {

        return ((obj.surface_link() == m_ref.surface_link()) &&
                is_same_scalar(obj.local_position()[0],
                               m_ref.local_position()[0], m_unc) &&
                is_same_scalar(obj.local_position()[1],
                               m_ref.local_position()[1], m_unc) &&
                is_same_scalar(obj.local_variance()[0],
                               m_ref.local_variance()[0], m_unc) &&
                is_same_scalar(obj.local_variance()[1],
                               m_ref.local_variance()[1], m_unc));
    }

    private:
    /// The reference object
    const edm::measurement<T> m_ref;
    /// The uncertainty
    scalar m_unc;

};  // class is_same_object<measurement>

}  // namespace traccc::details
