// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/device/algorithm_base.hpp"

// Project include(s).
#include <vecmem/containers/data/vector_view.hpp>

#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/seeding/detail/strip_pair.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/memory_resource.hpp"
#include "traccc/utils/messaging.hpp"

namespace traccc::device {

/// Separate standard (opposite-side) and overlap strip spacepoints.
struct strip_spacepoint_formation_output {
  edm::spacepoint_collection::buffer spacepoints;
  edm::spacepoint_collection::buffer overlap_spacepoints;
};

/// Algorithm forming space points out of measurements
///
/// This algorithm forms 3D spacepoints from compatible pairs of 1D
/// measurements. Static surface descriptors and directed pairing rules are
/// adapter inputs.
///
class silicon_strip_spacepoint_formation_algorithm
    : public algorithm<strip_spacepoint_formation_output(
          const edm::measurement_collection::const_view&,
          const strip_measurement_surface_info_collection_types::const_view&,
          const strip_pairing_rule_collection_types::const_view&,
          const point3&)>,
      public messaging,
      public algorithm_base {
 public:
  /// Constructor for spacepoint_formation algorithm
  ///
  /// @param mr The memory resource(s) to use in the algorithm
  /// @param copy The copy object to use for copying data between device
  ///             and host memory blocks
  /// @param logger The logger instance to use
  ///
  silicon_strip_spacepoint_formation_algorithm(
      const traccc::memory_resource& mr, const vecmem::copy& copy,
      std::unique_ptr<const Logger> logger = getDummyLogger().clone());

  /// Construct spacepoints from pairs of 1D silicon strip measurements
  ///
  /// @param measurements A collection of measurements
  /// @param surface_infos Static strip surface information
  /// @param pairing_rules Sorted, unique directed surface-pair rules
  /// @param beam_spot Beam-spot position for strip-plane construction
  /// @return Separate standard and overlap spacepoint buffers
  ///
  output_type operator()(
      const edm::measurement_collection::const_view& measurements,
      const strip_measurement_surface_info_collection_types::const_view&
          surface_infos,
      const strip_pairing_rule_collection_types::const_view& pairing_rules,
      const point3& beam_spot) const override;

 protected:
  /// @name Function(s) to be implemented by derived classes
  /// @{

  /// Payload for the @c count_strip_pairs_kernel function.
  struct count_strip_pairs_kernel_payload {
    /// The number of measurements in the event.
    edm::measurement_collection::const_view::size_type n_measurements;
    /// The input measurements.
    const edm::measurement_collection::const_view& measurements;
    /// Static strip surface information.
    const strip_measurement_surface_info_collection_types::const_view&
        surface_infos;
    const strip_pairing_rule_collection_types::const_view& pairing_rules;
    const point3& beam_spot;
    /// Number of opposite-side strip pairs.
    unsigned int& n_opposite_pairs;
    /// Number of overlap strip pairs.
    unsigned int& n_overlap_pairs;
  };

  /// Launch the strip pair counting kernel.
  virtual void count_strip_pairs_kernel(
      const count_strip_pairs_kernel_payload& payload) const = 0;

  /// Payload for the @c find_strip_pairs_kernel function.
  struct find_strip_pairs_kernel_payload {
    /// The number of measurements in the event.
    edm::measurement_collection::const_view::size_type n_measurements;
    /// The input measurements.
    const edm::measurement_collection::const_view& measurements;
    /// Static strip surface information.
    const strip_measurement_surface_info_collection_types::const_view&
        surface_infos;
    const strip_pairing_rule_collection_types::const_view& pairing_rules;
    /// Beam-spot position.
    const point3& beam_spot;
    /// Next positions in the two output pair buffers.
    unsigned int& opposite_position;
    unsigned int& overlap_position;
    /// Opposite-side and overlap strip pairs.
    strip_pair_collection_types::view& opposite_pairs;
    strip_pair_collection_types::view& overlap_pairs;
  };

  /// Launch the strip pair finding kernel.
  virtual void find_strip_pairs_kernel(
      const find_strip_pairs_kernel_payload& payload) const = 0;

  /// Payload for the @c form_spacepoints_kernel function
  struct form_spacepoints_kernel_payload {
    /// The number of compatible strip pairs in the event.
    strip_pair_collection_types::const_view::size_type n_pairs;
    /// The input measurements.
    const edm::measurement_collection::const_view& measurements;
    /// The compatible strip pairs.
    const strip_pair_collection_types::const_view& pairs;
    /// Static strip surface information.
    const strip_measurement_surface_info_collection_types::const_view&
        surface_infos;
    /// Beam-spot position.
    const point3& beam_spot;
    /// The output spacepoints.
    edm::spacepoint_collection::view& spacepoints;
  };

  /// Launch the spacepoint formation kernel
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void form_spacepoints_kernel(
      const form_spacepoints_kernel_payload& payload) const = 0;

  /// @}

};  // class silicon_strip_spacepoint_formation_algorithm

}  // namespace traccc::device
