/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/performance/details/is_same_object.hpp"
#include "traccc/performance/details/is_same_scalar.hpp"

// Project include(s).
#include "traccc/edm/seed_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"

namespace traccc::details {

/// @c traccc::is_same_object specialisation for @c traccc::edm::seed
template <typename T>
class is_same_object<edm::seed<T>> {

    public:
    /// Constructor with all necessary arguments
    is_same_object(
        const edm::spacepoint_collection::const_view& ref_spacepoints,
        const edm::spacepoint_collection::const_view& test_spacepoints,
        const edm::seed<T>& ref, scalar unc = float_epsilon)
        : m_ref_spacepoints(ref_spacepoints),
          m_spacepoints(test_spacepoints),
          m_ref(ref),
          m_unc(unc) {}

    /// Specialised implementation for @c traccc::seed
    bool operator()(const edm::seed<T>& obj) const {

        // Access the reference and test spacepoints.
        const edm::spacepoint_collection::const_device ref_spacepoints{
            m_ref_spacepoints};
        const edm::spacepoint ref_bottom_spacepoint =
            ref_spacepoints.at(m_ref.bottom_index());
        const edm::spacepoint ref_middle_spacepoint =
            ref_spacepoints.at(m_ref.middle_index());
        const edm::spacepoint ref_top_spacepoint =
            ref_spacepoints.at(m_ref.top_index());

        const edm::spacepoint_collection::const_device test_spacepoints{
            m_spacepoints};
        const edm::spacepoint test_bottom_spacepoint =
            test_spacepoints.at(obj.bottom_index());
        const edm::spacepoint test_middle_spacepoint =
            test_spacepoints.at(obj.middle_index());
        const edm::spacepoint test_top_spacepoint =
            test_spacepoints.at(obj.top_index());

        // Compare the spacepoints.
        return (is_same_object<
                    edm::spacepoint_collection::const_device::const_proxy_type>(
                    ref_bottom_spacepoint, m_unc)(test_bottom_spacepoint) &&
                is_same_object<
                    edm::spacepoint_collection::const_device::const_proxy_type>(
                    ref_middle_spacepoint, m_unc)(test_middle_spacepoint) &&
                is_same_object<
                    edm::spacepoint_collection::const_device::const_proxy_type>(
                    ref_top_spacepoint, m_unc)(test_top_spacepoint));
    }

    private:
    /// Spacepoints for the reference object
    const edm::spacepoint_collection::const_view m_ref_spacepoints;
    /// Spacepoint container for the test seeds
    const edm::spacepoint_collection::const_view m_spacepoints;

    /// The reference object
    const edm::seed<T> m_ref;
    /// The uncertainty
    const scalar m_unc;

};  // class is_same_object<seed>

}  // namespace traccc::details
