/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/performance/details/is_same_object.hpp"
#include "traccc/performance/details/projector.hpp"

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"

// System include(s).
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <string_view>
#include <vector>

namespace traccc {

template <typename TYPE>
soa_comparator<TYPE>::soa_comparator(
    std::string_view type_name,
    details::comparator_factory<typename TYPE::const_device::const_proxy_type>
        comp_factory,
    std::string_view lhs_type, std::string_view rhs_type, std::ostream& out,
    const std::vector<float>& uncertainties)
    : m_type_name(type_name),
      m_lhs_type(lhs_type),
      m_rhs_type(rhs_type),
      m_comp_factory(comp_factory),
      m_out(out),
      m_uncertainties(uncertainties) {

    std::sort(m_uncertainties.begin(), m_uncertainties.end());
}

template <typename TYPE>
void soa_comparator<TYPE>::operator()(
    const typename TYPE::const_view& lhs_view,
    const typename TYPE::const_view& rhs_view) const {

    // Create device containers on top of the views.
    const typename TYPE::const_device lhs{lhs_view}, rhs{rhs_view};

    // Print some basic output.
    m_out.get() << "Number of " << m_type_name << ": " << lhs.size() << " ("
                << m_lhs_type << "), " << rhs.size() << " (" << m_rhs_type
                << ")\n";

    // Calculate the agreements at various uncertainties.
    std::vector<float> agreements;
    agreements.reserve(m_uncertainties.size());

    std::vector<bool> is_matched(lhs.size(), false);
    for (float uncertainty : m_uncertainties) {
        // The number of matched items between the containers.
        std::size_t matched = 0;
        // Iterate over all elements of the LHS collection.
        for (typename TYPE::const_device::size_type i = 0; i < lhs.size();
             ++i) {
            if (is_matched[i]) {
                ++matched;
            } else {
                // Check if there's an equivalent element in the RHS collection.
                const auto comparator =
                    m_comp_factory.make_comparator(lhs[i], uncertainty);
                for (typename TYPE::const_device::size_type j = 0;
                     j < rhs.size(); ++j) {
                    if (comparator(rhs[j])) {
                        ++matched;
                        is_matched[i] = true;
                        break;
                    }
                }
            }
        }
        // Calculate the agreement value.
        agreements.push_back(
            static_cast<float>(matched) /
            static_cast<float>(std::max(lhs.size(), rhs.size())) * 100.f);
    }
    assert(agreements.size() == m_uncertainties.size());

    // Now print them.
    m_out.get() << "  Matching rate(s):\n";
    for (std::size_t i = 0; i < m_uncertainties.size(); ++i) {
        m_out.get() << "    - " << agreements.at(i) << "% at "
                    << m_uncertainties.at(i) * 100. << "% uncertainty\n";
    }
    m_out.get() << std::flush;
}

}  // namespace traccc
