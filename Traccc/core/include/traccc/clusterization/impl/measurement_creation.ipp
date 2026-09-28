/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/utils/detray_conversion.hpp"

// System include(s).
#include <cassert>

namespace traccc::details {

template <typename TCell, typename TDesign>
TRACCC_HOST_DEVICE inline vector2 position_from_cell(
    const edm::silicon_cell<TCell>& cell,
    const traccc::detector_design_description_interface<TDesign>& module_dd) {
  // Calculate / construct the local cell position.
  vector2 cell_lower_position = {(module_dd.bin_edges_x()).at(cell.channel0()),
                                 (module_dd.bin_edges_y()).at(cell.channel1())};

  vector2 cell_upper_position = {
      (module_dd.bin_edges_x()).at(cell.channel0() + 1),
      (module_dd.bin_edges_y()).at(cell.channel1() + 1)};

  vector2 cell_middle_position = {
      scalar{0.5f} * (cell_upper_position[0] + cell_lower_position[0]),
      scalar{0.5f} * (cell_upper_position[1] + cell_lower_position[1])};

  return cell_middle_position;
}

}  // namespace traccc::details
