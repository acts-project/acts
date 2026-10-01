// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/cuda/utils/algorithm_base.hpp"

// Project include(s).
#include "traccc/bfield/magnetic_field.hpp"
#include "traccc/edm/track_container.hpp"
#include "traccc/finding/actors/expected_layer_pattern_collector.hpp"
#include "traccc/finding/finding_config.hpp"
#include "traccc/geometry/detector_buffer.hpp"
#include "traccc/utils/memory_resource.hpp"

// Standard library include(s).
#include <functional>

namespace traccc::cuda {

/// Collect expected detector-layer patterns from finalized tracks.
class expected_layer_pattern_algorithm : public cuda::algorithm_base {
 public:
  using output_type = vecmem::data::vector_buffer<expected_layer_pattern_type>;

  expected_layer_pattern_algorithm(const finding_config& config,
                                   const traccc::memory_resource& mr,
                                   const vecmem::copy& copy,
                                   const stream_wrapper& str);

  /// Return one device-resident pattern per input track.
  output_type operator()(
      const detector_buffer& det, const magnetic_field& bfield,
      const edm::track_container<default_algebra>::buffer& tracks,
      vecmem::data::vector_view<const expected_layer_mapping_entry>
          expected_layer_map) const;

 private:
  finding_config m_config;
  traccc::memory_resource m_mr;
  std::reference_wrapper<const vecmem::copy> m_copy;
};

}  // namespace traccc::cuda
