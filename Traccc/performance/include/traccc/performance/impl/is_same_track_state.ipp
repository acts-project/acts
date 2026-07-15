/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_state_collection.hpp"

// Local include(s).
#include "traccc/performance/details/is_same_scalar.hpp"
#include "traccc/performance/impl/is_same_measurement.ipp"
#include "traccc/performance/impl/is_same_track_parameters.ipp"

namespace traccc::details {

/// @c traccc::is_same_object specialisation for
/// @c traccc::edm::track_state_collection
template <typename T>
class is_same_object<edm::track_state<T>> {

    public:
    /// Constructor with a reference object, and an allowed uncertainty
    is_same_object(const edm::measurement_collection::const_view& ref_meas,
                   const edm::measurement_collection::const_view& test_meas,
                   const edm::track_state<T>& ref, scalar unc = float_epsilon)
        : m_ref_meas(ref_meas),
          m_test_meas(test_meas),
          m_ref(ref),
          m_unc(unc) {}

    /// Specialised implementation for @c traccc::edm::track_state<T>
    bool operator()(const edm::track_state<T>& obj) const {

        // Compare the state words.
        if (obj.state() != m_ref.state()) {
            return false;
        }
        // Compare the chi2 values.
        if (!is_same_scalar(obj.filtered_chi2(), m_ref.filtered_chi2(),
                            m_unc) ||
            !is_same_scalar(obj.smoothed_chi2(), m_ref.smoothed_chi2(),
                            m_unc) ||
            !is_same_scalar(obj.backward_chi2(), m_ref.backward_chi2(),
                            m_unc)) {
            return false;
        }
        // Compare the fitted parameters.
        if (!is_same_object<bound_track_parameters<>>(
                m_ref.filtered_params(), m_unc)(obj.filtered_params()) ||
            !is_same_object<bound_track_parameters<>>(
                m_ref.smoothed_params(), m_unc)(obj.smoothed_params())) {
            return false;
        }
        // Compare the measurements that they point at.
        const edm::measurement_collection::const_device ref_meas{m_ref_meas};
        const edm::measurement_collection::const_device test_meas{m_test_meas};
        if (!is_same_object<
                edm::measurement_collection::const_device::const_proxy_type>(
                ref_meas.at(m_ref.measurement_index()),
                m_unc)(test_meas.at(obj.measurement_index()))) {
            return false;
        }

        // If we got here, the two track states are the same.
        return true;
    }

    private:
    /// Measurements for the reference object
    const edm::measurement_collection::const_view m_ref_meas;
    /// Measurements for the test object
    const edm::measurement_collection::const_view m_test_meas;
    /// The reference object
    const edm::track_state<T> m_ref;
    /// The uncertainty
    scalar m_unc;

};  // class is_same_object<edm::track_state<T>>

}  // namespace traccc::details
