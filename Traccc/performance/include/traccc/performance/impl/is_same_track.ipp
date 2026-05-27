/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/track_collection.hpp"
#include "traccc/edm/track_state_collection.hpp"

// Local include(s).
#include "traccc/performance/impl/is_same_measurement.ipp"
#include "traccc/performance/impl/is_same_track_state.ipp"

// System include(s).
#include <stdexcept>

namespace traccc::details {

/// @c traccc::is_same_object specialisation for @c traccc::edm::track<T>
template <typename T>
class is_same_object<edm::track<T>> {

    public:
    /// Constructor with a reference object, and an allowed uncertainty
    is_same_object(
        const edm::measurement_collection::const_view& ref_meas,
        const edm::measurement_collection::const_view& test_meas,
        const edm::track_state_collection<default_algebra>::const_view&
            ref_states,
        const edm::track_state_collection<default_algebra>::const_view&
            test_states,
        const edm::track<T>& ref, scalar unc = float_epsilon)
        : m_ref_meas(ref_meas),
          m_test_meas(test_meas),
          m_ref_states(ref_states),
          m_test_states(test_states),
          m_ref(ref),
          m_unc(unc) {}

    /// Specialised implementation for @c traccc::edm::track<T>
    bool operator()(const edm::track<T>& obj) const {

        // Compare the fit outcomes.
        if (obj.fit_outcome() != m_ref.fit_outcome()) {
            return false;
        }
        // Compare the parameters.
        if (!is_same_object<bound_track_parameters<>>(m_ref.params(),
                                                      m_unc)(obj.params()) ||
            !is_same_scalar(obj.ndf(), m_ref.ndf(), m_unc)) {
            return false;
        }
        // Compare the fitted parameters. But only in case of a successful fit.
        if ((m_ref.fit_outcome() == track_fit_outcome::SUCCESS) &&
            (!is_same_scalar(obj.chi2(), m_ref.chi2(), m_unc) ||
             !is_same_scalar(obj.pval(), m_ref.pval(), m_unc))) {
            return false;
        }
        // Compare the number of holes.
        if (obj.nholes() != m_ref.nholes()) {
            return false;
        }

        // The two tracks need to have the same number of constituents.
        if (obj.constituent_links().size() !=
            m_ref.constituent_links().size()) {
            return false;
        }

        // Now compare the constituents one by one.
        const edm::measurement_collection::const_device ref_meas{m_ref_meas};
        const edm::measurement_collection::const_device test_meas{m_test_meas};
        const edm::track_state_collection<default_algebra>::const_device
            ref_states{m_ref_states};
        const edm::track_state_collection<default_algebra>::const_device
            test_states{m_test_states};
        for (unsigned int i = 0; i < obj.constituent_links().size(); ++i) {

            // Make sure that both links point at the same type of object.
            if (obj.constituent_links()[i].type !=
                m_ref.constituent_links()[i].type) {
                return false;
            }
            // Now compare the constituents themselves.
            if (obj.constituent_links()[i].type ==
                edm::track_constituent_link::measurement) {

                if (!is_same_object<edm::measurement_collection::const_device::
                                        const_proxy_type>(
                        ref_meas.at(m_ref.constituent_links()[i].index), m_unc)(
                        test_meas.at(obj.constituent_links()[i].index))) {
                    return false;
                }

            } else if (obj.constituent_links()[i].type ==
                       edm::track_constituent_link::track_state) {

                if (!is_same_object<edm::track_state_collection<
                        default_algebra>::const_device::const_proxy_type>(
                        m_ref_meas, m_test_meas,
                        ref_states.at(m_ref.constituent_links()[i].index),
                        m_unc)(
                        test_states.at(obj.constituent_links()[i].index))) {
                    return false;
                }

            } else {
                // Unknown constituent type!
                throw std::runtime_error(
                    "Unknown track constituent link type encountered!");
            }
        }

        // If we got here, the two track fits are the same.
        return true;
    }

    private:
    /// Measurements for the reference object
    const edm::measurement_collection::const_view m_ref_meas;
    /// Measurements for the test object
    const edm::measurement_collection::const_view m_test_meas;
    /// States for the reference object
    const edm::track_state_collection<default_algebra>::const_view m_ref_states;
    /// States for the test object
    const edm::track_state_collection<default_algebra>::const_view
        m_test_states;
    /// The reference object
    const edm::track<T> m_ref;
    /// The uncertainty
    scalar m_unc;

};  // class is_same_object<edm::track<T>>

}  // namespace traccc::details
