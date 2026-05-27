/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/performance/details/is_same_object.hpp"

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/container.hpp"

// System include(s).
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <string_view>
#include <vector>

namespace traccc {

template <typename HEADER_TYPE, typename ITEM_TYPE>
container_comparator<HEADER_TYPE, ITEM_TYPE>::container_comparator(
    std::string_view type_name,
    details::comparator_factory<HEADER_TYPE> header_comp_factory,
    details::comparator_factory<ITEM_TYPE> item_comp_factory,
    std::string_view lhs_type, std::string_view rhs_type, std::ostream& out,
    const std::vector<scalar>& uncertainties)
    : m_type_name(type_name),
      m_lhs_type(lhs_type),
      m_rhs_type(rhs_type),
      m_header_comp_factory(header_comp_factory),
      m_item_comp_factory(item_comp_factory),
      m_out(out),
      m_uncertainties(uncertainties) {}

template <typename HEADER_TYPE, typename ITEM_TYPE>
void container_comparator<HEADER_TYPE, ITEM_TYPE>::operator()(
    const typename container_types<HEADER_TYPE, ITEM_TYPE>::const_view& lhs,
    const typename container_types<HEADER_TYPE, ITEM_TYPE>::const_view& rhs)
    const {

    // Create device containers on top of the views.
    const typename container_types<HEADER_TYPE, ITEM_TYPE>::const_device
        lhs_cont{lhs},
        rhs_cont{rhs};

    // Print some basic output.
    m_out.get() << "Number of " << m_type_name << ": " << lhs_cont.total_size()
                << " (" << m_lhs_type << "), " << rhs_cont.total_size() << " ("
                << m_rhs_type << ")\n";

    // Calculate the agreements at various uncertainties.
    std::vector<scalar> agreements;
    agreements.reserve(m_uncertainties.size());
    for (scalar uncertainty : m_uncertainties) {
        // The number of matched items between the containers.
        std::size_t matched = 0;
        // If the two containers differ in size, only compare the first N
        // "inner vectors", that both of them have.
        const std::size_t cont_size =
            std::min(lhs_cont.get_items().size(), rhs_cont.get_items().size());
        // Iterate over the "outer vectors".
        for (std::size_t i = 0; i < cont_size; ++i) {
            // If the headers don't match, don't even compare the items.
            if (m_header_comp_factory.make_comparator(
                    lhs_cont.get_headers().at(i),
                    uncertainty)(rhs_cont.get_headers().at(i)) == false) {
                continue;
            }
            // Compare the items.
            const typename collection_types<ITEM_TYPE>::const_device lhs_items =
                lhs_cont.get_items().at(i);
            const typename collection_types<ITEM_TYPE>::const_device rhs_items =
                rhs_cont.get_items().at(i);
            for (const ITEM_TYPE& obj : lhs_items) {
                if (std::find_if(rhs_items.begin(), rhs_items.end(),
                                 m_item_comp_factory.make_comparator(
                                     obj, uncertainty)) != rhs_items.end()) {
                    ++matched;
                }
            }
        }
        // Calculate the agreement value.
        agreements.push_back(
            static_cast<scalar>(matched) /
            static_cast<scalar>(
                std::max(lhs_cont.total_size(), rhs_cont.total_size())) *
            100.);
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
