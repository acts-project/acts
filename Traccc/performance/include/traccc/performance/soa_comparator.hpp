/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/performance/details/comparator_factory.hpp"

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"

// System include(s).
#include <functional>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

namespace traccc {

/// Functor comparing two SoA containers, and printing the results
///
/// This is meant to be used in the example applications for nicely comparing
/// the results made on the host and on a device. Though the code actually
/// allows comparisons between any two SoA containers.
///
/// @tparam TYPE The trait/meta type describing the SoA container
///
template <typename TYPE>
class soa_comparator {

    public:
    /// Constructor with all configurable options
    soa_comparator(std::string_view type_name,
                   details::comparator_factory<
                       typename TYPE::const_device::const_proxy_type>
                       comp_factory = {},
                   std::string_view lhs_type = "host",
                   std::string_view rhs_type = "device",
                   std::ostream& out = std::cout,
                   const std::vector<float>& uncertainties = {0.0001f, 0.001f,
                                                              0.01f, 0.05f});

    /// Function comparing two collections, and printing the results
    void operator()(const typename TYPE::const_view& lhs,
                    const typename TYPE::const_view& rhs) const;

    private:
    /// Container type name to print
    std::string m_type_name;
    /// Type of the "Left Hand Side" collection
    std::string m_lhs_type;
    /// Type of the "Right Hand Side" collection
    std::string m_rhs_type;

    /// Factory for making comparator objects
    details::comparator_factory<typename TYPE::const_device::const_proxy_type>
        m_comp_factory;

    /// Output stream to print the results to
    std::reference_wrapper<std::ostream> m_out;

    /// Uncertainties to evaluate the comparison for
    std::vector<float> m_uncertainties;

};  // class collection_comparator

}  // namespace traccc

// Include the implementation.
#include "traccc/performance/impl/soa_comparator.ipp"
