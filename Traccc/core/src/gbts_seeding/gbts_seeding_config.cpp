/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "traccc/gbts_seeding/gbts_seeding_config.hpp"

#include <algorithm>
#include <ranges>

namespace traccc {

// binTables contains pairs of linked layer-eta bins
// the layerInfo should really be calculated from the geoIDBinning
// GeoIDBinning pair is detray geo ID and bin index (corrisponding to the
// layers in layerInfo) minPt in MeV
bool gbts_seedfinder_config::setLinkingScheme(
    const std::vector<std::pair<int, std::vector<int>>>& input_binTables,
    const device::gbts_layerInfo input_layerInfo,
    std::vector<std::pair<uint64_t, short>>& detrayGeoIDBinning,
    float minPt = 900.0f,
    std::unique_ptr<const traccc::Logger> callers_logger =
        getDummyLogger().clone()) {

    TRACCC_LOCAL_LOGGER(std::move(callers_logger));
    // copy layer-eta binning infomation
    layerInfo = input_layerInfo;
    // unroll binTables
    for (std::pair<int, std::vector<int>> binPairs : input_binTables) {
        for (int bin2 : binPairs.second) {
            binTables.push_back(
                std::make_pair(static_cast<unsigned int>(binPairs.first),
                               static_cast<unsigned int>(bin2)));
        }
    }

    for (std::pair<int, int> lI : layerInfo.info)
        n_eta_bins = std::max(n_eta_bins,
                              static_cast<unsigned int>(lI.first + lI.second));

    // bin by volume
    std::ranges::sort(
        detrayGeoIDBinning,
        [](const std::pair<uint64_t, short> a,
           const std::pair<uint64_t, short> b) { return a.first > b.first; });

    unsigned int largest_volume_index =
        detray::geometry::identifier(detrayGeoIDBinning[0].first).volume();
    auto current_volume = static_cast<short>(largest_volume_index);
    if (largest_volume_index >= SHRT_MAX) {
        TRACCC_ERROR(
            "volume index to large to fit in type-short GBTS volume to layer "
            "map");
        return false;
    }

    bool layerChange = false;
    short current_layer = detrayGeoIDBinning[0].second;

    int split_volumes = 0;
    std::vector<std::pair<short, unsigned int>> volumeToLayerMap_unordered;
    detrayGeoIDBinning.push_back(
        std::make_pair(UINT_MAX, -1));  // end-of-vector element
    std::vector<std::array<unsigned int, 2>> surfacesInVolume;
    for (std::pair<uint64_t, short> geoIDLayerPair : detrayGeoIDBinning) {
        detray::geometry::identifier geo_id(geoIDLayerPair.first);
        if (current_volume != static_cast<short>(geo_id.volume())) {
            // reached the end of this volume so add it to the maps
            short bin = current_layer;
            if (layerChange) {
                split_volumes++;
                bin = -1 *
                      static_cast<short>(
                          surfaceToLayerMap.size() +
                          1);  // start of this volume's surfaces in the map + 1
                for (std::array<unsigned int, 2> pair : surfacesInVolume)
                    surfaceToLayerMap.push_back(pair);
            }
            volumeToLayerMap_unordered.push_back(std::make_pair(
                bin, current_volume));  // layerIdx if not split, begin-index in
                                        // the surface map otherwise

            current_volume = static_cast<short>(geo_id.volume());
            current_layer = geoIDLayerPair.second;
            layerChange = false;
            surfacesInVolume.clear();
        }
        // is volume encompassed by a layer
        layerChange |= (current_layer != geoIDLayerPair.second);

        // save surfaces incase volume is not encommpassed by a layer
        surfacesInVolume.push_back(std::array<unsigned int, 2>{
            static_cast<unsigned int>(geo_id.index()),
            static_cast<unsigned int>(geoIDLayerPair.second)});
    }
    // make volume by layer map
    volumeToLayerMap.resize(largest_volume_index + 1);
    for (unsigned int i = 0; i < largest_volume_index + 1; ++i)
        volumeToLayerMap.push_back(SHRT_MAX);
    for (std::pair<short, unsigned int> vLpair : volumeToLayerMap_unordered)
        volumeToLayerMap[vLpair.second] = vLpair.first;
    // scale cuts
    float ptScale = 900.0f / minPt;
    graph_building_params.min_delta_phi *= ptScale;
    graph_building_params.dphi_coeff *= ptScale;
    graph_building_params.min_delta_phi_low_dr *= ptScale;
    graph_building_params.dphi_coeff_low_dr *= ptScale;
    graph_building_params.max_Kappa *= ptScale;

    // contianers sizes
    nLayers = static_cast<unsigned int>(layerInfo.type.size());

    TRACCC_INFO("volume layer map has " << volumeToLayerMap_unordered.size()
                                        << " volumes");
    TRACCC_INFO("The maxium volume index in the layer map is "
                << volumeToLayerMap.size());
    TRACCC_INFO("surface to layer map has " << surfaceToLayerMap.size()
                                            << " geo IDs from " << split_volumes
                                            << " multi-layer volumes");
    TRACCC_INFO("layer info found for " << nLayers << " layers");
    TRACCC_INFO(binTables.size() << " linked layer-eta bins for GBTS");

    if (nLayers == 0) {
        TRACCC_ERROR("no layers input");
        return false;
    } else if (volumeToLayerMap.size() == 0) {
        TRACCC_ERROR("empty volume to layer map");
        return false;
    } else if (surfaceToLayerMap.size() > SHRT_MAX) {
        TRACCC_ERROR("surface to layer map is to large");
        return false;
    }
    return true;
}

}  // namespace traccc
