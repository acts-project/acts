// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "../utils/cuda_error_handling.hpp"
#include "../utils/magnetic_field_types.hpp"
#include "../utils/utils.hpp"
#include "./kernels/collect_expected_layer_patterns.cuh"
#include "traccc/cuda/finding/expected_layer_pattern_algorithm.hpp"

// Project include(s).
#include "traccc/utils/detector_buffer_bfield_visitor.hpp"

namespace traccc::cuda {

expected_layer_pattern_algorithm::expected_layer_pattern_algorithm(
    const finding_config& config, const traccc::memory_resource& mr,
    const vecmem::copy& copy, const stream_wrapper& str)
    : algorithm_base(str), m_config(config), m_mr(mr), m_copy(copy) {}

expected_layer_pattern_algorithm::output_type
expected_layer_pattern_algorithm::operator()(
    const detector_buffer& det, const magnetic_field& bfield,
    const edm::track_container<default_algebra>::buffer& tracks,
    vecmem::data::vector_view<const expected_layer_mapping_entry>
        expected_layer_map) const {
  const auto n_tracks = m_copy.get().get_size(tracks.tracks);
  output_type output_patterns(n_tracks, m_mr.main);
  m_copy.get().setup(output_patterns)->wait();

  if (n_tracks == 0u || expected_layer_map.size() == 0u) {
    if (n_tracks != 0u) {
      TRACCC_CUDA_ERROR_CHECK(
          cudaMemsetAsync(output_patterns.ptr(), 0,
                          n_tracks * sizeof(expected_layer_pattern_type),
                          details::get_stream(stream())));
      stream().synchronize();
    }
    return output_patterns;
  }

  const auto map_size = static_cast<unsigned int>(expected_layer_map.size());
  vecmem::data::vector_buffer<expected_layer_mapping_entry> device_map(
      map_size, m_mr.main);
  m_copy.get().setup(device_map)->wait();
  m_copy
      .get()(vecmem::data::vector_view<const expected_layer_mapping_entry>(
                 map_size, expected_layer_map.ptr()),
             device_map)
      ->wait();

  const vecmem::data::vector_view<const expected_layer_mapping_entry>
      device_map_view(map_size, device_map.ptr());
  detector_buffer_magnetic_field_visitor<detector_type_list,
                                         cuda::bfield_type_list<scalar>>(
      det, bfield,
      [&]<detray::concepts::detector detector_t, typename bfield_view_t>(
          const detray::detector_view_t<detector_t>& detector,
          const bfield_view_t& field) {
        constexpr unsigned int n_threads = 128u;
        const unsigned int n_blocks =
            static_cast<unsigned int>((n_tracks + n_threads - 1u) / n_threads);
        collect_expected_layer_patterns<detray::detector_device_t<detector_t>,
                                        bfield_view_t>(
            n_blocks, n_threads, 0u, details::get_stream(stream()), detector,
            field, {tracks}, m_config, device_map_view, output_patterns);
      });

  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
  stream().synchronize();
  return output_patterns;
}

}  // namespace traccc::cuda
