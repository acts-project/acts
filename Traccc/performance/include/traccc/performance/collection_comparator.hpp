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

/// Functor comparing two collections, and printing the results to the output
///
/// This is meant to be used in the example applications for nicely comparing
/// the results made on the host and on a device. Though the code actually
/// allows comparisons between any two containers.
///
/// @tparam TYPE The type in the collection
///
template <typename TYPE>
class collection_comparator {

    public:
    /// Constructor with all configurable options
    collection_comparator(std::string_view type_name,
                          details::comparator_factory<TYPE> comp_factory = {},
                          std::string_view lhs_type = "host",
                          std::string_view rhs_type = "device",
                          std::ostream& out = std::cout,
                          const std::vector<scalar>& uncertainties = {
                              0.0001f, 0.001f, 0.01f, 0.05f});

    /// Function comparing two collections, and printing the results
    void operator()(
        const typename collection_types<TYPE>::const_view& lhs,
        const typename collection_types<TYPE>::const_view& rhs) const;

    private:
    /// Container type name to print
    std::string m_type_name;
    /// Type of the "Left Hand Side" collection
    std::string m_lhs_type;
    /// Type of the "Right Hand Side" collection
    std::string m_rhs_type;

    /// Factory for making comparator objects
    details::comparator_factory<TYPE> m_comp_factory;

    /// Output stream to print the results to
    std::reference_wrapper<std::ostream> m_out;

    /// Uncertainties to evaluate the comparison for
    std::vector<scalar> m_uncertainties;

};  // class collection_comparator

}  // namespace traccc

// Include the implementation.
#include "traccc/performance/impl/collection_comparator.ipp"
