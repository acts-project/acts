/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "read_digitization_config.hpp"

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

traccc::digitization_config read_digitization_config(
    const nlohmann::json& json) {
    traccc::digitization_config result;

    static const char* entries_key = "entries";
    static const char* binningdata_key = "binningdata";
    static const char* geometric = "geometric";
    static const char* segmentation = "segmentation";

    std::vector<traccc::digitization_config::InputElement> elements;

    for (const auto& entry : json[entries_key]) {

        Acts::GeometryIdentifier geoId;
        Acts::GeometryIdentifier::Value null(0u);
        geoId = geoId.withVolume(entry.value("volume", null))
                    .withLayer(entry.value("layer", null))
                    .withSensitive(entry.value("sensitive", null));

        const auto& json_val = entry["value"];
        const auto& json_geom = json_val[geometric];
        const auto& json_segm = json_geom[segmentation];
        std::vector<std::vector<float>> bin_edges;
        unsigned char dimensions = 2;
        for (const auto& bindata : json_segm[binningdata_key]) {
            std::vector<float> bins;

            if (bindata["bins"].get<int>() == 1) {
                dimensions = 1;
                if (bindata["type"].get<std::string>() == "equidistant") {
                    float pitch = (bindata["max"].get<float>() -
                                   bindata["min"].get<float>()) /
                                  bindata["bins"].get<float>();
                    for (int i = 0; i <= bindata["bins"].get<int>(); ++i) {
                        bins.push_back(bindata["min"].get<float>() +
                                       static_cast<float>(i) * pitch);
                    }
                } else {
                    for (const auto& edge : bindata["edges"]) {
                        bins.push_back(edge.get<float>());
                    }
                }
            } else {
                if (bindata["type"].get<std::string>() == "equidistant") {
                    float pitch = (bindata["max"].get<float>() -
                                   bindata["min"].get<float>()) /
                                  bindata["bins"].get<float>();
                    for (int i = 0; i <= bindata["bins"].get<int>(); ++i) {
                        bins.push_back(bindata["min"].get<float>() +
                                       static_cast<float>(i) * pitch);
                    }
                } else {
                    for (const auto& edge : bindata["edges"]) {
                        bins.push_back(edge.get<float>());
                    }
                }
            }

            bin_edges.push_back(bins);
        }

        elements.push_back({geoId, {bin_edges, dimensions}});
    }

    return traccc::digitization_config(std::move(elements));
}

namespace io::json {

digitization_config read_digitization_config(std::string_view filename) {
    // Open the input file. Relying on exceptions for the error handling.
    std::ifstream infile(filename.data(), std::ifstream::binary);
    infile.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    // Read the contents of the file into a JSON object.
    nlohmann::json json;
    infile >> json;

    // Construct the object from the JSON configuration.
    return traccc::read_digitization_config(json);
}

}  // namespace io::json
}  // namespace traccc
