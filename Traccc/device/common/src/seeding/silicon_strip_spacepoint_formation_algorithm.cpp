/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/seeding/device/silicon_strip_spacepoint_formation_algorithm.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_buffer.hpp>
#include <vecmem/containers/vector.hpp>

namespace traccc::device {

silicon_strip_spacepoint_formation_algorithm::
    silicon_strip_spacepoint_formation_algorithm(
        const traccc::memory_resource& mr, const vecmem::copy& copy,
        std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)), algorithm_base(mr, copy) {}

auto silicon_strip_spacepoint_formation_algorithm::operator()(
    const detector_buffer& det,
    const edm::measurement_collection::const_view& measurements,
    const strip_measurement_surface_info_collection_types::const_view&
        surface_infos,
    const strip_pairing_rule_collection_types::const_view& pairing_rules,
    const point3& beam_spot) const -> output_type {
  edm::measurement_collection::const_view::size_type n_measurements = 0u;
  if (mr().host) {
    vecmem::async_size size = copy().get_size(measurements, *(mr().host));
    n_measurements = size.get();
  } else {
    n_measurements = copy().get_size(measurements);
  }

  if (n_measurements == 0u) {
    return {};
  }

  vecmem::data::vector_buffer<unsigned int> pair_counter_buffer(2u, mr().main);
  copy().setup(pair_counter_buffer)->ignore();
  copy().memset(pair_counter_buffer, 0)->ignore();
  count_strip_pairs_kernel(
      {n_measurements, det, measurements, surface_infos, pairing_rules,
       beam_spot, pair_counter_buffer.ptr()[0], pair_counter_buffer.ptr()[1]});

  vecmem::vector<unsigned int> pair_counter_host(mr().host ? mr().host
                                                           : &(mr().main));
  copy()(pair_counter_buffer, pair_counter_host)->wait();
  const unsigned int n_opposite_pairs = pair_counter_host.at(0);
  const unsigned int n_overlap_pairs = pair_counter_host.at(1);

  strip_pair_collection_types::buffer opposite_pairs_buffer(n_opposite_pairs,
                                                            mr().main);
  strip_pair_collection_types::buffer overlap_pairs_buffer(n_overlap_pairs,
                                                           mr().main);
  copy().setup(opposite_pairs_buffer)->ignore();
  copy().setup(overlap_pairs_buffer)->ignore();

  if ((n_opposite_pairs + n_overlap_pairs) > 0u) {
    copy().memset(pair_counter_buffer, 0)->ignore();
    find_strip_pairs_kernel(
        {n_measurements, det, measurements, surface_infos, pairing_rules,
         beam_spot, pair_counter_buffer.ptr()[0], pair_counter_buffer.ptr()[1],
         opposite_pairs_buffer, overlap_pairs_buffer});
  }

  edm::spacepoint_collection::buffer opposite_spacepoints(
      n_opposite_pairs, mr().main, vecmem::data::buffer_type::resizable);
  edm::spacepoint_collection::buffer overlap_spacepoints(
      n_overlap_pairs, mr().main, vecmem::data::buffer_type::resizable);
  copy().setup(opposite_spacepoints)->ignore();
  copy().setup(overlap_spacepoints)->ignore();

  if (n_opposite_pairs > 0u) {
    form_spacepoints_kernel({n_opposite_pairs, det, measurements,
                             opposite_pairs_buffer, surface_infos, beam_spot,
                             opposite_spacepoints});
  }
  if (n_overlap_pairs > 0u) {
    form_spacepoints_kernel({n_overlap_pairs, det, measurements,
                             overlap_pairs_buffer, surface_infos, beam_spot,
                             overlap_spacepoints});
  }

  return {std::move(opposite_spacepoints), std::move(overlap_spacepoints)};
}

}  // namespace traccc::device
