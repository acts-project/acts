// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/seeding/device/silicon_strip_spacepoint_formation_algorithm.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_buffer.hpp>
#include <vecmem/containers/vector.hpp>

// System include(s).
#include <cassert>

namespace traccc::device {

silicon_strip_spacepoint_formation_algorithm::
    silicon_strip_spacepoint_formation_algorithm(
        const traccc::memory_resource& mr, const vecmem::copy& copy,
        std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)), algorithm_base(mr, copy) {}

auto silicon_strip_spacepoint_formation_algorithm::operator()(
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

  assert(input_is_sorted(measurements));

  // Write one count per measurement, then scan to give each measurement its
  // own deterministic range. Pair order follows measurement, rule and candidate
  // traversal order, independently of GPU scheduling.
  vecmem::data::vector_buffer<unsigned int> standard_offsets_buffer(
      n_measurements, mr().main);
  vecmem::data::vector_buffer<unsigned int> overlap_offsets_buffer(
      n_measurements, mr().main);
  copy().setup(standard_offsets_buffer)->ignore();
  copy().setup(overlap_offsets_buffer)->ignore();
  vecmem::data::vector_view<unsigned int> standard_offsets(
      standard_offsets_buffer);
  vecmem::data::vector_view<unsigned int> overlap_offsets(
      overlap_offsets_buffer);
  count_strip_pairs_kernel({n_measurements, measurements, surface_infos,
                            pairing_rules, beam_spot, standard_offsets,
                            overlap_offsets});
  scan_offsets(standard_offsets);
  scan_offsets(overlap_offsets);

  // Read just the last inclusive offset to size the corresponding output.
  const auto read_total = [&](vecmem::data::vector_view<unsigned int> offsets) {
    vecmem::vector<unsigned int> total(mr().host ? mr().host : &mr().main);
    copy()(vecmem::data::vector_view<const unsigned int>(
               1u, offsets.ptr() + offsets.capacity() - 1u),
           total)
        ->wait();
    return total.at(0);
  };
  const unsigned int n_opposite_pairs = read_total(standard_offsets);
  const unsigned int n_overlap_pairs = read_total(overlap_offsets);

  strip_pair_collection_types::buffer opposite_pairs_buffer(n_opposite_pairs,
                                                            mr().main);
  strip_pair_collection_types::buffer overlap_pairs_buffer(n_overlap_pairs,
                                                           mr().main);
  copy().setup(opposite_pairs_buffer)->ignore();
  copy().setup(overlap_pairs_buffer)->ignore();

  // Count and find may disagree at floating-point cut boundaries. Initialise
  // unused slots to invalid indices; formation rejects them before access.
  // Find is bounded by each measurement's scanned range, so an excess pair
  // cannot overwrite the range assigned to another measurement.
  copy().memset(opposite_pairs_buffer, 0xff)->ignore();
  copy().memset(overlap_pairs_buffer, 0xff)->ignore();
  if ((n_opposite_pairs > 0u) || (n_overlap_pairs > 0u)) {
    find_strip_pairs_kernel({n_measurements, measurements, surface_infos,
                             pairing_rules, beam_spot, standard_offsets,
                             overlap_offsets, opposite_pairs_buffer,
                             overlap_pairs_buffer});
  }

  const auto form_and_gather =
      [&](unsigned int n_pairs,
          const strip_pair_collection_types::const_view& pairs) {
        if (n_pairs == 0u) {
          return edm::spacepoint_collection::buffer{};
        }
        // Formation writes to the pair's own slot. Scan its acceptance flag and
        // gather surviving candidates into a fixed-size, ordered output buffer.
        edm::spacepoint_collection::buffer candidates(n_pairs, mr().main);
        vecmem::data::vector_buffer<unsigned int> accepted_buffer(n_pairs,
                                                                  mr().main);
        copy().setup(candidates)->ignore();
        copy().setup(accepted_buffer)->ignore();
        vecmem::data::vector_view<unsigned int> accepted(accepted_buffer);
        form_spacepoints_kernel({n_pairs, measurements, pairs, surface_infos,
                                 beam_spot, accepted, candidates});
        scan_offsets(accepted);
        const unsigned int n_spacepoints = read_total(accepted);
        if (n_spacepoints == 0u) {
          return edm::spacepoint_collection::buffer{};
        }
        edm::spacepoint_collection::buffer result(n_spacepoints, mr().main);
        copy().setup(result)->ignore();
        const edm::spacepoint_collection::const_view candidates_view(
            candidates);
        gather_spacepoints_kernel({n_pairs, candidates_view, accepted, result});
        return result;
      };
  auto opposite_spacepoints =
      form_and_gather(n_opposite_pairs, opposite_pairs_buffer);
  auto overlap_spacepoints =
      form_and_gather(n_overlap_pairs, overlap_pairs_buffer);

  return {std::move(opposite_spacepoints), std::move(overlap_spacepoints)};
}

}  // namespace traccc::device
