/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "traccc/performance/details/is_same_object.hpp"

#include "traccc/performance/details/is_same_angle.hpp"
#include "traccc/performance/details/is_same_scalar.hpp"

namespace traccc::details {

/// @name Implementation for
///       @c traccc::details::is_same_object<bound_track_parameters>
/// @{
is_same_object<bound_track_parameters<>>::is_same_object(
    const bound_track_parameters<>& ref, scalar unc)
    : m_ref(ref), m_unc(unc) {}

bool is_same_object<bound_track_parameters<>>::operator()(
    const bound_track_parameters<>& obj) const {

    return ((obj.surface_link() == m_ref.get().surface_link()) &&
            is_same_scalar(obj.bound_local()[0], m_ref.get().bound_local()[0],
                           m_unc) &&
            is_same_scalar(obj.bound_local()[1], m_ref.get().bound_local()[1],
                           m_unc) &&
            is_same_angle(obj.phi(), m_ref.get().phi(), m_unc) &&
            is_same_scalar(obj.theta(), m_ref.get().theta(), m_unc) &&
            is_same_scalar(obj.time(), m_ref.get().time(), m_unc) &&
            is_same_scalar(obj.qop(), m_ref.get().qop(), m_unc));
}

/// @}

}  // namespace traccc::details
