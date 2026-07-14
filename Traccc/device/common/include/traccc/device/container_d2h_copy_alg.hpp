/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/memory_resource.hpp"
#include "traccc/utils/messaging.hpp"

// VecMem include(s).
#include <vecmem/utils/copy.hpp>

namespace traccc::device {

/// Algorithm to copy a container from a device to the host
///
/// Specialisations of this algorithm can be used to get a container back
/// from a device to the host. Ideally at the end of an algorithm chain that
/// would've produced something in device memory.
///
/// Note that the algorithm needs to treat constant views and buffers a little
/// differently, so the user should make sure to pass in the correct type.
///
/// @tparam CONTAINER_TYPES One of the "container types" traits
///
template <typename CONTAINER_TYPES>
class container_d2h_copy_alg
    : public algorithm<typename CONTAINER_TYPES::host(
          const typename CONTAINER_TYPES::const_view&)>,
      public messaging {

    public:
    /// Helper type declaration for the input type
    typedef const typename CONTAINER_TYPES::const_view& input_type;
    /// Help the compiler understand what @c output_type is
    using output_type = typename algorithm<typename CONTAINER_TYPES::host(
        const typename CONTAINER_TYPES::const_view&)>::output_type;

    /// Constructor with the needed resources
    container_d2h_copy_alg(
        const memory_resource& mr, vecmem::copy& deviceCopy,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());

    /// Function executing the copy to the host
    virtual output_type operator()(input_type input) const override;

    private:
    /// The memory resource(s) to use
    memory_resource m_mr;
    /// The D->H copy object to use
    vecmem::copy& m_deviceCopy;
    /// The H->H copy object to use
    vecmem::copy m_hostCopy;
};  // class container_d2h_copy_alg

}  // namespace traccc::device

// Include the implementation.
#include "traccc/device/impl/container_d2h_copy_alg.ipp"
