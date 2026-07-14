/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "write_digitization_config.hpp"

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

nlohmann::json module_digi_config_to_json(
    const Acts::GeometryIdentifier& geoId,
    const module_digitization_config& cfg) {

    static const char* geometric = "geometric";
    static const char* segmentation = "segmentation";
    static const char* binningdata_key = "binningdata";

    static const char* bin_value_labels[] = {"binX", "binY"};

    nlohmann::json entry;

    // Write geometry identifier fields — only write non-zero ones
    // ie the top most value in geometry sructure
    if (geoId.volume() != 0)
        entry["volume"] = geoId.volume();
    if (geoId.layer() != 0)
        entry["layer"] = geoId.layer();
    if (geoId.sensitive() != 0)
        entry["sensitive"] = geoId.sensitive();

    nlohmann::json binning_array = nlohmann::json::array();

    for (std::size_t dim = 0; dim < cfg.bin_edges.size(); ++dim) {
        const auto& edges = cfg.bin_edges[dim];
        int nbins = static_cast<int>(edges.size()) - 1;

        nlohmann::json bindata;
        bindata["bins"] = nbins;
        bindata["option"] = "open";
        bindata["value"] = (dim < 2) ? bin_value_labels[dim] : "binX";

        if (!edges.empty()) {
            bindata["min"] = static_cast<float>(edges.front());
            bindata["max"] = static_cast<float>(edges.back());
        }

        // Correctly write equidistant vs irregular binning
        bool equidistant = true;
        if (nbins > 1) {
            float expected_pitch =
                (edges.back() - edges.front()) / static_cast<float>(nbins);
            for (std::size_t i = 1; i < edges.size(); ++i) {
                if (std::abs((edges[i] - edges[i - 1]) - expected_pitch) >
                    1e-5f) {
                    equidistant = false;
                    break;
                }
            }
        }

        if (equidistant) {
            bindata["type"] = "equidistant";
        } else {
            bindata["type"] = "variable";
            bindata["edges"] = edges;
        }

        binning_array.push_back(bindata);
    }

    entry["value"][geometric][segmentation][binningdata_key] = binning_array;

    return entry;
}

nlohmann::json to_json(const digitization_config& config) {
    nlohmann::json json;

    // Top-level header — must match what Acts GeometryHierarchyMap expects
    json["acts-geometry-hierarchy-map"] = {
        {"format-version", 0},
        {"value-identifier", "digitization-configuration"}};

    nlohmann::json entries = nlohmann::json::array();

    for (unsigned int i = 0; i < config.size(); i++) {
        entries.push_back(
            module_digi_config_to_json(config.idAt(i), config.valueAt(i)));
    }

    json["entries"] = entries;
    return json;
}

namespace io::json {

void write_digitization_config(std::string_view filename,
                               const digitization_config& config) {
    // Construct the JSON object.
    nlohmann::json json = traccc::to_json(config);

    // Open the output file. Relying on exceptions for error handling.
    std::ofstream outfile(filename.data(), std::ofstream::binary);
    outfile.exceptions(std::ofstream::failbit | std::ofstream::badbit);

    // Write the JSON object to the file.
    outfile << json.dump(4);
}

}  // namespace io::json
}  // namespace traccc
