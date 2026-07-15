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
                   const std::string_view& geometry_file,
                   const std::string_view& material_file,
                   const std::string_view& grid_file,
                   const bool do_consistency_check = false) {

    // Set up the detector reader configuration.
    detray::io::detector_reader_config cfg;
    cfg.do_check(do_consistency_check);

    cfg.add_file(traccc::io::get_absolute_path(geometry_file));
    if (material_file.empty() == false) {
        cfg.add_file(traccc::io::get_absolute_path(material_file));
    }
    if (grid_file.empty() == false) {
        cfg.add_file(traccc::io::get_absolute_path(grid_file));
    }

    // Read the detector.
    auto det = detray::io::read_detector<typename detector_t::host>(mr, cfg);
    detector.set<detector_t>(std::move(det.first));
}

}  // namespace

namespace traccc::io {

void read_detector(host_detector& detector, vecmem::memory_resource& mr,
                   const std::string_view& geometry_file,
                   const std::string_view& material_file,
                   const std::string_view& grid_file) {

    // Peek at the header to determine the kind of detector that is needed
    const auto header = detray::io::detail::deserialize_json_header(
        traccc::io::get_absolute_path(geometry_file));

    // TODO: Update this
    if (header.detector == "Cylindrical detector from DD4hep blueprint") {
        ::read_detector<odd_detector>(detector, mr, geometry_file,
                                      material_file, grid_file);
    } else if (header.detector == "detray_detector") {
        ::read_detector<itk_detector>(detector, mr, geometry_file,
                                      material_file, grid_file);
    } else {
        // TODO: Warning here
        ::read_detector<default_detector>(detector, mr, geometry_file,
                                          material_file, grid_file);
    }
}

}  // namespace traccc::io
