/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/cuda/utils/algorithm_base.hpp"

// Project include(s).
#include "traccc/seeding/device/silicon_strip_spacepoint_formation_algorithm.hpp"

namespace traccc::cuda {

/// Algorithm forming space points out of measurements
///
/// This algorithm forms 3D spacepoints from compatible pairs of 1D strip
/// measurements, using adapter-provided static geometry and pairing rules.
///
class silicon_strip_spacepoint_formation_algorithm
    : public device::silicon_strip_spacepoint_formation_algorithm,
      public cuda::algorithm_base {
 public:
  /// Constructor for the spacepoint formation algorithm
  ///
  /// @param mr The memory resource(s) to use in the algorithm
  /// @param copy The copy object to use for copying data between device
  ///             and host memory blocks
  /// @param str The CUDA stream to use
  /// @param logger The logger instance to use
  ///
  silicon_strip_spacepoint_formation_algorithm(
      const traccc::memory_resource& mr, const vecmem::copy& copy,
      const stream_wrapper& str,
      std::unique_ptr<const Logger> logger = getDummyLogger().clone());

 private:
  /// @name Function(s) inherited from
  /// @c traccc::device::silicon_strip_spacepoint_formation_algorithm
  /// @{

  /// Launch the strip pair counting kernel.
  void count_strip_pairs_kernel(
      const count_strip_pairs_kernel_payload& payload) const override;

  /// Launch the strip pair finding kernel.
  void find_strip_pairs_kernel(
      const find_strip_pairs_kernel_payload& payload) const override;

  /// Launch the spacepoint formation kernel
  ///
  /// @param payload The payload for the kernel
  ///
  void form_spacepoints_kernel(
      const form_spacepoints_kernel_payload& payload) const override;

  /// @}

};  // class silicon_strip_spacepoint_formation_algorithm

}  // namespace traccc::cuda
