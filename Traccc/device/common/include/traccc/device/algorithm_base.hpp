/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/utils/memory_resource.hpp"

// VecMem include(s).
#include <vecmem/utils/copy.hpp>

// System include(s).
#include <functional>

namespace traccc::device {

/// Base class for all of the device algorithms
///
/// It would hold on to members that all device algorithms would need.
///
class algorithm_base {

    public:
    /// Constructor for the device algorithm base class
    ///
    /// @param mr     The memory resource(s) to use
    /// @param copy   The copy object to use
    ///
    algorithm_base(const memory_resource& mr, vecmem::copy& copy);

    /// The memory resource(s) to use in the algorithm (non-const)
    memory_resource& mr();
    /// The memory resource(s) to use in the algorithm (const)
    const memory_resource& mr() const;

    /// The copy object to use in the algorithm (non-const)
    vecmem::copy& copy();
    /// The copy object to use in the algorithm (const)
    const vecmem::copy& copy() const;

    private:
    /// Memory resource(s) to use in the algorithm
    memory_resource m_mr;
    /// Copy object to use in the algorithm
    std::reference_wrapper<vecmem::copy> m_copy;

};  // class algorithm_base

}  // namespace traccc::device
