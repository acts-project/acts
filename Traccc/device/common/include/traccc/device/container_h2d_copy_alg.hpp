/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/utils/memory_resource.hpp"
#include "traccc/utils/messaging.hpp"

// VecMem include(s).
#include <vecmem/utils/copy.hpp>

// System include(s).
#include <type_traits>

namespace traccc::device {

/// Algorithm to copy a container from the host to the device
///
/// Specialisations of this algorithm can be used to copy a container to a
/// device for a subsequent algorithm that expects its input to already be on
/// the device that it wants to run on.
///
/// Note that the input type is not a "host type", but rather a constant view.
/// This is meant to allow using this algorithm on top of more complicated
/// data objects.
///
/// @tparam CONTAINER_TYPES One of the "container types" traits
///
template <typename CONTAINER_TYPES>
class container_h2d_copy_alg : public messaging {

    public:
    /// Helper type declaration for the input type
    typedef const typename CONTAINER_TYPES::const_view& input_type;
    /// Helper type declaration for the output type
    typedef typename CONTAINER_TYPES::buffer output_type;

    /// Constructor with the needed resources
    container_h2d_copy_alg(
        const memory_resource& mr, vecmem::copy& deviceCopy,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone());

    /// Function executing a simple copy to the device
    output_type operator()(input_type input) const;
    /// Function executing an optimised copy to the device
    output_type operator()(input_type input,
                           typename CONTAINER_TYPES::buffer& hostBuffer) const;

    private:
    /// Size type for the handled container's header vector
    using header_size_type =
        const typename std::remove_reference<typename std::remove_cv<
            input_type>::type>::type::header_vector::size_type;

    /// Helper function calculating the size(s) of the input container
    std::vector<std::size_t> get_sizes(input_type input) const;

    /// The memory resource(s) to use
    memory_resource m_mr;
    /// The H->D copy object to use
    vecmem::copy& m_deviceCopy;
    /// The H->H copy object to use
    vecmem::copy m_hostCopy;
};  // class container_h2d_copy_alg

}  // namespace traccc::device

// Include the implementation.
#include "traccc/device/impl/container_h2d_copy_alg.ipp"
