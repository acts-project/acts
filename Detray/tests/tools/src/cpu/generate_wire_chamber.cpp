// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/shapes/line.hpp"

// Detray IO include(s)
#include "detray/io/frontend/detector_writer.hpp"

// Detray test include(s)
#include "detray/options/detector_io_options.hpp"
#include "detray/options/parse_options.hpp"
#include "detray/options/wire_chamber_options.hpp"
#include "detray/test/common/build_wire_chamber.hpp"
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>

// Boost
#include "detray/options/boost_program_options.hpp"

namespace po = boost::program_options;

namespace detray::detail {

/// Build and write the wire chamber according to command line options
template <concepts::scalar scalar_t, typename wire_shape_t>
void write_wire_chamber(int argc, char **argv,
                        boost::program_options::options_description &desc,
                        vecmem::memory_resource *host_mr,
                        detray::io::detector_writer_config &writer_cfg) {
  detray::wire_chamber_config<scalar_t, wire_shape_t> wire_cfg{};

  po::variables_map vm_wire =
      detray::options::parse_options(desc, argc, argv, wire_cfg);

  auto [wire_chamber, names] =
      build_wire_chamber<test::algebra>(*host_mr, wire_cfg);
  detray::io::write_detector(wire_chamber, names, writer_cfg);
}

}  // namespace detray::detail

using namespace detray;

int main(int argc, char **argv) {
  // Options parsing
  po::options_description desc("\nWire chamber generation options");

  desc.add_options()("write_volume_graph", "writes the volume graph to file");
  desc.add_options()("straw_tubes", "build the detector using straw tubes");
  desc.add_options()("drift_cells",
                     "build the detector using wire chamber cells");

  // Configuration
  detray::io::detector_writer_config writer_cfg{};
  writer_cfg.format(detray::io::format::json).replace_files(false);
  // Default output path
  writer_cfg.path("./wire_chamber/");

  // Parse general options
  po::variables_map vm =
      detray::options::parse_options(desc, argc, argv, writer_cfg);

  // General options
  if (vm.count("write_volume_graph")) {
    throw std::invalid_argument("Writing of volume graph not implemented");
  }

  // Determine the wire shape
  if (vm.count("straw_tubes") && vm.count("drift_cells")) {
    throw std::invalid_argument(
        "Cannot build straw tubes and wire cells together");
  }

  // Build the geometry
  vecmem::host_memory_resource host_mr;
  if (vm.count("straw_tubes")) {
    detail::write_wire_chamber<detray::test::scalar, detray::line_circular>(
        argc, argv, desc, &host_mr, writer_cfg);
  } else {
    detail::write_wire_chamber<detray::test::scalar, detray::line_square>(
        argc, argv, desc, &host_mr, writer_cfg);
  }
}
