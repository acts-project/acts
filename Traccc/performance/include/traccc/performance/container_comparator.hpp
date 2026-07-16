/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/performance/details/comparator_factory.hpp"

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/edm/container.hpp"

// System include(s).
#include <functional>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

namespace traccc {

/// Functor comparing two containers, and printing the results to the output
///
/// This is meant to be used in the example applications for nicely comparing
/// the results made on the host and on a device. Though the code actually
/// allows comparisons between any two containers.
///
/// @tparam HEADER_TYPE The header type in the container
/// @tparam ITEM_TYPE The item type in the container
///
template <typename HEADER_TYPE, typename ITEM_TYPE>
class container_comparator {

    public:
    /// Constructor with all configurable options
    container_comparator(
        std::string_view type_name,
        details::comparator_factory<HEADER_TYPE> header_comp_factory = {},
        details::comparator_factory<ITEM_TYPE> item_comp_factory = {},
        std::string_view lhs_type = "host",
        std::string_view rhs_type = "device", std::ostream& out = std::cout,
        const std::vector<scalar>& uncertainties = {0.0001f, 0.001f, 0.01f,
                                                    0.05f});

    /// Function comparing two collections, and printing the results
    void operator()(
        const typename container_types<HEADER_TYPE, ITEM_TYPE>::const_view& lhs,
        const typename container_types<HEADER_TYPE, ITEM_TYPE>::const_view& rhs)
        const;

    private:
    /// Container type name to print
    std::string m_type_name;
    /// Type of the "Left Hand Side" collection
    std::string m_lhs_type;
    /// Type of the "Right Hand Side" collection
    std::string m_rhs_type;

    /// Factory for making comparator objects for header objects
    details::comparator_factory<HEADER_TYPE> m_header_comp_factory;
    /// Factory for making comparator objects for item objects
    details::comparator_factory<ITEM_TYPE> m_item_comp_factory;

    /// Output stream to print the results to
    std::reference_wrapper<std::ostream> m_out;

    /// Uncertainties to evaluate the comparison for
    std::vector<scalar> m_uncertainties;

};  // class container_comparator

}  // namespace traccc

// Include the implementation.
#include "traccc/performance/impl/container_comparator.ipp"
