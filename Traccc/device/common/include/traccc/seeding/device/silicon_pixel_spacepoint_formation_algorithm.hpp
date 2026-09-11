/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/device/algorithm_base.hpp"

// Project include(s).
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/geometry/detector_buffer.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/memory_resource.hpp"
#include "traccc/utils/messaging.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

/// Algorithm forming space points out of measurements
///
/// This algorithm performs the local-to-global transformation of the 2D (pixel)
/// measurements made on every detector module, into 3D spacepoint coordinates.
///
class silicon_pixel_spacepoint_formation_algorithm
    : public algorithm<edm::spacepoint_collection::buffer(
          const detector_buffer&,
          const edm::measurement_collection::const_view&)>,
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
  silicon_pixel_spacepoint_formation_algorithm(
      const traccc::memory_resource& mr, const vecmem::copy& copy,
      std::unique_ptr<const Logger> logger = getDummyLogger().clone());

  /// Construct spacepoints from 2D silicon pixel measurements
  ///
  /// @param det Detector object
  /// @param measurements A collection of measurements
  /// @return A spacepoint buffer, with one spacepoint for every
  ///         silicon pixel measurement
  ///
  output_type operator()(const detector_buffer& det,
                         const edm::measurement_collection::const_view&
                             measurements) const override;

 protected:
  /// @name Function(s) to be implemented by derived classes
  /// @{

  /// Payload for the @c count_spacepoints_kernel function
  struct count_spacepoints_kernel_payload {
    /// The number of measurements in the event
    edm::measurement_collection::const_view::size_type n_measurements;
    /// The input measurements
    const edm::measurement_collection::const_view& measurements;
    /// One output flag for every measurement
    vecmem::data::vector_view<unsigned int>& spacepoint_flags;
  };

  /// Launch the spacepoint counting kernel
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void count_spacepoints_kernel(
      const count_spacepoints_kernel_payload& payload) const = 0;

  /// Turn the spacepoint flags into a prefix sum, in place
  ///
  /// The scan must be inclusive, so that the last element of @c
  /// spacepoint_flags holds the total number of spacepoints.
  ///
  /// @param spacepoint_flags The flags to scan
  ///
  virtual void scan_spacepoint_flags(
      vecmem::data::vector_view<unsigned int>& spacepoint_flags) const = 0;

  /// Payload for the @c form_spacepoints_kernel function
  struct form_spacepoints_kernel_payload {
    /// The number of measurements in the event
    edm::measurement_collection::const_view::size_type n_measurements;
    /// The detector object
    const detector_buffer& detector;
    /// The input measurements
    const edm::measurement_collection::const_view& measurements;
    /// The prefix sum giving every measurement its spacepoint index
    const vecmem::data::vector_view<const unsigned int>& spacepoint_index;
    /// The output spacepoints
    edm::spacepoint_collection::view& spacepoints;
  };

  /// Launch the spacepoint formation kernel
  ///
  /// @param payload The payload for the kernel
  ///
  virtual void form_spacepoints_kernel(
      const form_spacepoints_kernel_payload& payload) const = 0;

  /// @}

};  // class silicon_pixel_spacepoint_formation_algorithm

}  // namespace traccc::device
