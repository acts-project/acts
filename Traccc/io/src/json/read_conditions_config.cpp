/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "read_conditions_config.hpp"

// Acts include(s).
#if __has_include(<ActsPlugins/Json/ActsJson.hpp>)
#include <ActsPlugins/Json/ActsJson.hpp>
#include <ActsPlugins/Json/GeometryHierarchyMapJsonConverter.hpp>
#include <ActsPlugins/Json/UtilitiesJsonConverter.hpp>
#else
#include <Acts/Plugins/Json/ActsJson.hpp>
#include <Acts/Plugins/Json/GeometryHierarchyMapJsonConverter.hpp>
#include <Acts/Plugins/Json/UtilitiesJsonConverter.hpp>
#endif

// System include(s).
#include <fstream>

namespace traccc {

/// Function allowing the read of @c traccc::module_conditions_config objects
///
/// Note that this function must be declared in the same namespace as
/// @c traccc::module_conditions_config for nlohmann_json to work correctly.
///
traccc::conditions_config read_conditions_config(const nlohmann::json& json) {
    traccc::conditions_config result;

    static const char* entries_key = "entries";
    static const char* binningdata_key = "binningdata";
    static const char* geometric = "geometric";
    static const char* segmentation = "segmentation";

    std::vector<traccc::conditions_config::InputElement> elements;

    for (const auto& entry : json[entries_key]) {

        Acts::GeometryIdentifier geoId;
        Acts::GeometryIdentifier::Value null(0u);
        geoId = geoId.withVolume(entry.value("volume", null))
                    .withLayer(entry.value("layer", null))
                    .withSensitive(entry.value("sensitive", null));

        const auto& json_val = entry["value"];
        const auto& json_geom = json_val[geometric];
        const auto& json_segm = json_geom[segmentation];
        const auto& json_binning = json_segm[binningdata_key];
        vector2 shift = {0.f, 0.f};

        if (json_binning.contains("shift")) {
            shift[0] = json_binning["shift"][0].get<float>();
            shift[1] = json_binning["shift"][1].get<float>();
        }

        elements.push_back({geoId, {shift}});
    }

    return traccc::conditions_config(std::move(elements));
}

namespace io::json {

conditions_config read_conditions_config(std::string_view filename) {
    conditions_config result;

    // Open the input file. Relying on exceptions for the error handling.
    std::ifstream infile(filename.data(), std::ifstream::binary);
    infile.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    // Read the contents of the file into a JSON object.
    nlohmann::json json;
    infile >> json;

    // Construct the object from the JSON configuration.
    return traccc::read_conditions_config(json);
}

}  // namespace io::json
}  // namespace traccc
