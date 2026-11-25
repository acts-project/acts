// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/core/concepts.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/ranges.hpp"

// Detray test include(s)
#include "detray/test/utils/data_record.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

namespace detray::validation {

/// Compare the result of two material tracers
///
/// @param reference the reference trace to compare against
/// @param ref_record the reference total material to compare against (s/X0)
/// @param mat_trace the recoded material trace
/// @param mat_record the recoded total material along the track (s/X0)
/// @param trk_i the track index/identifier
/// @param rel_tol the maximum allowed relative error for any parameters
/// @param verbose verbose debug output
///
/// @returns whether the traces contain matching material steps and whether
/// matching traces contain the same total material
inline auto compare_material_traces(
    const dvector<material_record> &reference, const track_material &ref_record,
    const dvector<material_record> &mat_trace, const track_material &mat_record,
    std::size_t trk_i = detail::invalid_value<std::size_t>(),
    const double rel_tol = 0.01, const bool verbose = true) {
  // Material traces contain records of different surfaces/material
  bool is_bad_steps{reference.size() != mat_trace.size()};
  // Overall recorded material is different (in case small errors accumulate)
  bool is_diff_mat{false};

  std::stringstream debug_msg{};
  debug_msg << "Track No. " << trk_i << ":\n----------------" << std::endl;

  if (is_bad_steps) {
    debug_msg << "-> Different no. of surfaces: " << mat_trace.size()
              << " (ref.: " << reference.size() << ")\n"
              << std::endl;
  } else {
    for (std::size_t j = 0u; j < math::min(reference.size(), mat_trace.size());
         ++j) {
      if (reference[j].geo_id != mat_trace[j].geo_id) {
        is_bad_steps = true;
        debug_msg << "-> Surfaces don't match: " << mat_trace[j].geo_id
                  << " (ref.: " << reference[j].geo_id << ")" << std::endl;
        continue;
      }

      // TODO: Use approx_equal from algebra-plugins
      // Compare thickness of the surface material
      if ((reference[j].thickness - mat_trace[j].thickness) /
              reference[j].thickness >
          rel_tol) {
        is_bad_steps = true;
        debug_msg << "-> On surface: " << reference[j].geo_id << ":"
                  << std::endl;
        debug_msg << "-> thickness: " << mat_trace[j].thickness
                  << ", thickness ref.: " << reference[j].thickness
                  << std::endl;
      }

      // Compare radiation length of the surface material
      if ((reference[j].mat_X0 - mat_trace[j].mat_X0) / reference[j].mat_X0 >
          rel_tol) {
        is_bad_steps = true;
        debug_msg << "-> On surface: " << reference[j].geo_id << ":"
                  << std::endl;
        debug_msg << "-> X0: " << mat_trace[j].mat_X0
                  << ", X0 ref.: " << reference[j].mat_X0 << std::endl;
      }

      // Compare interaction length of the surface material
      if ((reference[j].mat_L0 - mat_trace[j].mat_L0) / reference[j].mat_L0 >
          rel_tol) {
        is_bad_steps = true;
        debug_msg << "-> On surface: " << reference[j].geo_id << ":"
                  << std::endl;
        debug_msg << "-> L0: " << mat_trace[j].mat_L0
                  << ", L0 ref.: " << reference[j].mat_L0 << std::endl;
      }

      // Compare path of the track through the surface material
      if ((reference[j].path - mat_trace[j].path) / reference[j].path >
          rel_tol) {
        is_bad_steps = true;
        debug_msg << "-> On surface: " << reference[j].geo_id << ":"
                  << std::endl;
        debug_msg << "-> Mat. path: " << mat_trace[j].path / unit<double>::mm
                  << " mm, mat. path ref.: "
                  << reference[j].path / unit<double>::mm << " mm" << std::endl;
      }
    }
  }

  // Check the accumulated material data
  auto check_discrepancy = [rel_tol, trk_i, &debug_msg](
                               const double truth, const double rec,
                               const std::string_view what) {
    constexpr auto e{std::numeric_limits<double>::epsilon()};
    auto rel_error{detail::invalid_value<double>()};

    if (truth <= e && rec <= e) {
      // No material for this ray => valid
      rel_error = 0.;
    } else if (truth <= e) {
      // Material found where none should be
      DETRAY_ERROR_HOST(
          "Found material on recorded material step, but none in truth trace");
    } else {
      rel_error = math::fabs(truth - rec) / truth;
    }
    if (!(rel_error < rel_tol)) {
      debug_msg << "Track " << trk_i << " (" << what << "): Truth " << truth
                << ", Nav. " << rec;
      debug_msg << "\n -> Total discrepancy of " << 100. * rel_error << "%\n"
                << std::endl;
      return false;
    }
    return true;
  };

  // Already know that material cannot add up: Only flag the additional cases
  if (is_bad_steps) {
    is_diff_mat = true;
  }

  is_diff_mat = is_diff_mat ||
                check_discrepancy(ref_record.sX0, mat_record.sX0, "X0 / path");
  is_diff_mat = is_diff_mat || check_discrepancy(ref_record.tX0, mat_record.tX0,
                                                 "X0 / thickness");
  is_diff_mat = is_diff_mat ||
                check_discrepancy(ref_record.sL0, mat_record.sL0, "L0 / path");
  is_diff_mat = is_diff_mat || check_discrepancy(ref_record.tL0, mat_record.tL0,
                                                 "L0 / thickness");

  if (verbose && (is_bad_steps || is_diff_mat)) {
    DETRAY_INFO_HOST(debug_msg.str());
  }

  return std::make_tuple(is_bad_steps, is_diff_mat);
}

}  // namespace detray::validation
