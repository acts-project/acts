/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <algorithm>
#include <cassert>
#include <type_traits>
#include <vector>

namespace traccc::device {

template <typename CONTAINER_TYPES>
container_d2h_copy_alg<CONTAINER_TYPES>::container_d2h_copy_alg(
    const memory_resource& mr, vecmem::copy& deviceCopy,
    std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)), m_mr(mr), m_deviceCopy(deviceCopy) {}

template <typename CONTAINER_TYPES>
typename container_d2h_copy_alg<CONTAINER_TYPES>::output_type
container_d2h_copy_alg<CONTAINER_TYPES>::operator()(input_type input) const {

    // A sanity check.
    assert(input.headers.size() == input.items.size());

    // Decide what memory resource to use for the host container.
    vecmem::memory_resource* host_mr =
        (m_mr.host != nullptr) ? m_mr.host : &(m_mr.main);

    // Copy the device container into 2 temporary host buffers.
    auto header_buffer = m_deviceCopy.to(input.headers, *host_mr,
                                         vecmem::copy::type::device_to_host);
    auto item_buffer = m_deviceCopy.to(input.items, *host_mr, nullptr,
                                       vecmem::copy::type::device_to_host);

    // Create the result object, giving it the appropriate memory resource for
    // all of its elements.
    output_type result{header_buffer.size(), host_mr};
    for (std::size_t i = 0; i < result.size(); ++i) {
        result[i].items =
            typename CONTAINER_TYPES::host::item_vector::value_type{host_mr};
    }

    // Perform the H->H copy.
    vecmem::copy::event_type host_header_copy_event =
        m_hostCopy(header_buffer, result.get_headers());
    vecmem::copy::event_type host_item_copy_event =
        m_hostCopy(item_buffer, result.get_items());
    host_header_copy_event->wait();
    host_item_copy_event->wait();

    // Return the host object.
    return result;
}

}  // namespace traccc::device
