/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/io/read_detector.hpp"

#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/io/detector.hpp"
#include "traccc/io/utils.hpp"

// System include(s).
#include <string>

namespace {

/// Common implementation for constructing a detector from a set of input files
template <typename detector_t>
void read_detector(traccc::host_detector& detector, vecmem::memory_resource& mr,
                   detray::io::detector_payload& payload,
                   const detray::io::detector_reader_config& cfg) {
  // Geometry payload is already filled
  assert(!payload.geometry.volumes.empty());

  // Read the detector.
  auto det =
      detray::io::read_detector<typename detector_t::host>(mr, cfg, payload);

  detector.set<detector_t>(std::move(det.first));
}

}  // namespace

namespace traccc::io {

void read_detector(host_detector& detector, vecmem::memory_resource& mr,
                   const std::string_view geometry_file,
                   const std::string_view material_file,
                   const std::string_view grid_file,
                   const bool do_consistency_check) {
  detray::io::detector_payload payload{};
  detray::io::convert_json_to_payload(
      payload, {traccc::io::get_absolute_path(geometry_file),
                traccc::io::get_absolute_path(material_file),
                traccc::io::get_absolute_path(grid_file)});

  // Build the detector (of a specific type) from payload
  detray::io::detector_reader_config cfg;
  cfg.do_check(do_consistency_check);

  // TODO: Update this and remove detector name
  std::string_view det_name{payload.detector_name};
  std::string_view metadata{payload.header.metadata};
  if (metadata == "odd_metadata" ||
      det_name == "Cylindrical detector from DD4hep blueprint") {
    ::read_detector<odd_detector>(detector, mr, payload, cfg);
  } else if (metadata == "itk_metadata" || det_name == "detray_detector") {
    ::read_detector<itk_detector>(detector, mr, payload, cfg);
  } else {
    // TODO: Warning here
    ::read_detector<default_detector>(detector, mr, payload, cfg);
  }
}

}  // namespace traccc::io
