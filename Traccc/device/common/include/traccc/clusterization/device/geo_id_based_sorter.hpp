/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"

// Detray include(s).
#include <detray/geometry/identifier.hpp>

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace traccc::device {

/// Functor to help with sorting measurements based on their surface identifiers
class geo_id_based_sorter {

    public:
    /// Constructor, capturing the geo IDs to sort indices by
    explicit geo_id_based_sorter(
        const vecmem::device_vector<const detray::geometry::identifier>&
            geo_ids)
        : m_geo_ids(geo_ids) {}

    /// Index comparison operator
    ///
    /// The logic is a bit convoluted here. :-( The index sequence/vector that
    /// we sort may be larger than the geo ID vector. (Whenever we are dealing)
    /// with resizable measurement collections. Which is often.) In this case
    /// the behaviour should be that for the non-existent geo IDs the indices
    /// should be left unchanged. Which the following logic should do...
    ///
    TRACCC_HOST_DEVICE bool operator()(unsigned int lhs,
                                       unsigned int rhs) const {

        if (lhs >= m_geo_ids.size()) {
            return false;
        }
        if (rhs >= m_geo_ids.size()) {
            return true;
        }
        return m_geo_ids.at(lhs) < m_geo_ids.at(rhs);
    }

    private:
    /// The geo IDs to sort indices by
    vecmem::device_vector<const detray::geometry::identifier> m_geo_ids;

};  // class geo_id_based_sorter

}  // namespace traccc::device
