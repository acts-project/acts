/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Library include(s).
#include "traccc/io/read_detector_description.hpp"

#include "traccc/geometry/host_detector.hpp"
#include "traccc/io/detector.hpp"
#include "traccc/io/read_conditions_config.hpp"
#include "traccc/io/read_detector.hpp"
#include "traccc/io/read_digitization_config.hpp"
#include "traccc/io/utils.hpp"

// Detray include(s)
#include <detray/geometry/tracking_surface.hpp>
#include <detray/utils/type_registry.hpp>

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s).
#include <iterator>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {

void fill_digi_info(traccc::detector_design_description::host& det_desc,
                    const traccc::module_digitization_config& data) {

    det_desc.dimensions().back() = data.dimensions;

    det_desc.bin_edges_x().back().assign(data.bin_edges[0].begin(),
                                         data.bin_edges[0].end());

    det_desc.bin_edges_y().back().assign(data.bin_edges[1].begin(),
                                         data.bin_edges[1].end());
}

template <typename detector_traits_t>
void read_json_dd_impl(traccc::detector_design_description::host& det_desc,
                       traccc::detector_conditions_description::host& det_cond,
                       const traccc::host_detector& detector,
                       const traccc::digitization_config& digi,
                       const traccc::conditions_config& cond)
    requires(traccc::is_detector_traits<detector_traits_t>)
{
    const typename detector_traits_t::host& detector_host =
        detector.as<detector_traits_t>();

    det_desc.reserve(digi.size());
    det_cond.reserve(detector_host.surfaces().size());

    int n_design = 0;
    for (const auto& digi_it : digi) {

        det_desc.resize(det_desc.size() + 1);
        det_desc.design_id().back() = n_design++;
        fill_digi_info(det_desc, digi_it);
    }

    std::vector<int> module_to_design;

    for (const auto& surface_desc : detector_host.surfaces()) {
        const traccc::geometry_id geom_id{surface_desc.source};
        const Acts::GeometryIdentifier acts_geom_id{geom_id};

        if (acts_geom_id.sensitive() == 0) {
            continue;
        }
        // New module — add it to detector conditions description
        det_cond.resize(det_cond.size() + 1);
        det_cond.geometry_id().back() = surface_desc.identifier();
        det_cond.acts_geometry_id().back() = geom_id;

        std::array<detray::dindex_type<traccc::default_algebra>, 2u> subspace =
            {0, 1};
        using annulus_t =
            detray::mask<detray::annulus2D, traccc::default_algebra>;
        using mask_registry_t = typename detector_traits_t::host::masks;
        if constexpr (detray::types::contains<mask_registry_t, annulus_t>) {
            if (surface_desc.mask().id() ==
                detray::types::id<mask_registry_t, annulus_t>) {
                subspace = {1, 0};
            }
        }
        if (!digi.contains(acts_geom_id)) {
            std::ostringstream msg;
            msg << "Could not find digitization config for geometry ID: "
                << acts_geom_id;
            throw std::runtime_error(msg.str());
        } else {
            auto digi_it = digi.find(acts_geom_id);
            if (digi_it != digi.end()) {
                int idx =
                    static_cast<int>(std::distance(digi.begin(), digi_it));
                det_desc.subspace()[static_cast<unsigned long>(idx)] = subspace;
                module_to_design.push_back(idx);
            }
        }

        // Find the module's conditions configuration.
        const traccc::conditions_config::Iterator cond_it =
            cond.find(acts_geom_id);
        if (cond_it == cond.end()) {
            std::ostringstream msg;
            msg << "Could not find conditions config for geometry ID: "
                << acts_geom_id;
            throw std::runtime_error(msg.str());
        }

        det_cond.measurement_translation().back() = cond_it->shift;
    }

    det_cond.module_to_design_id().assign(module_to_design.begin(),
                                          module_to_design.end());
}

void read_json_dd(traccc::detector_design_description::host& det_desc,
                  traccc::detector_conditions_description::host& det_cond,
                  std::string_view geometry_file,
                  const traccc::digitization_config& digi,
                  const traccc::conditions_config& cond) {

    // Construct a (temporary) Detray detector object from the geometry
    // configuration file.
    vecmem::host_memory_resource mr;
    traccc::host_detector detector;
    traccc::io::read_detector(detector, mr, geometry_file);

    // TODO: Implement detector visitor!
    // Peek at the header to determine the kind of detector that is needed
    const auto header = detray::io::detail::deserialize_json_header(
        traccc::io::get_absolute_path(geometry_file));

    if (header.detector == "Cylindrical detector from DD4hep blueprint") {
        read_json_dd_impl<traccc::odd_detector>(det_desc, det_cond, detector,
                                                digi, cond);
    } else if (header.detector == "detray_detector") {
        read_json_dd_impl<traccc::itk_detector>(det_desc, det_cond, detector,
                                                digi, cond);
    } else {
        // TODO: Warning here
        read_json_dd_impl<traccc::default_detector>(det_desc, det_cond,
                                                    detector, digi, cond);
    }
}

}  // namespace

namespace traccc::io {

void read_detector_description(detector_design_description::host& det_desc,
                               detector_conditions_description::host& det_cond,
                               std::string_view geometry_file,
                               std::string_view digitization_file,
                               std::string_view conditions_file,
                               const data_format geometry_format,
                               const data_format digitization_format,
                               const data_format conditions_format) {

    // Read the digitization configuration.
    const digitization_config digi =
        read_digitization_config(digitization_file, digitization_format);
    const conditions_config cond =
        read_conditions_config(conditions_file, conditions_format);
    // Fill the detector description with the correct type of geometry file.
    switch (geometry_format) {
        case data_format::json:
            ::read_json_dd(det_desc, det_cond, geometry_file, digi, cond);
            break;
        default:
            throw std::invalid_argument("Unsupported geometry format.");
    }
}

}  // namespace traccc::io
