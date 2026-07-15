/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "traccc/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"

#include "silicon_pixel_spacepoint_formation.hpp"

namespace traccc::host {

silicon_pixel_spacepoint_formation_algorithm::
    silicon_pixel_spacepoint_formation_algorithm(
        vecmem::memory_resource& mr, std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)), m_mr(mr) {}

silicon_pixel_spacepoint_formation_algorithm::output_type
silicon_pixel_spacepoint_formation_algorithm::operator()(
    const host_detector& det,
    const edm::measurement_collection::const_view& meas) const {

    return host_detector_visitor<detector_type_list>(
        det, [&]<typename detector_traits_t>(
                 const typename detector_traits_t::host& detector) {
            return details::silicon_pixel_spacepoint_formation(detector, meas,
                                                               m_mr);
        });
}

}  // namespace traccc::host
