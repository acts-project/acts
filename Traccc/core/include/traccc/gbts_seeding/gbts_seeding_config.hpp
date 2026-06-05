/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s)
#include <memory>

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/utils/messaging.hpp"

// Detray include(s).
#include <detray/geometry/identifier.hpp>

namespace traccc::device {

struct gbts_layerInfo {
    std::vector<char> type;
    // etaBin0 and numBins
    std::vector<std::pair<int, int>> info;
    // minEta and deltaEta
    std::vector<std::pair<float, float>> geo;

    void reserve(unsigned int n) {
        type.reserve(n);
        info.reserve(n);
        geo.reserve(n);
    }

    void addLayer(char layerType, int firstBin, int nBins, float minEta,
                  float etaBinWidth) {
        type.push_back(layerType);
        info.push_back(std::make_pair(firstBin, nBins));
        geo.push_back(std::make_pair(minEta, etaBinWidth));
    }
};

struct gbts_consts {

    // CCA max iterations -> maxium seed length
    static constexpr unsigned short max_cca_iter = 15;
    // shared memory allocation sizes
    static constexpr unsigned short node_buffer_length = 128;
    static constexpr unsigned short live_path_buffer = 1024;

    // access into output graph
    static constexpr char node1 = 0;
    static constexpr char node2 = 1;
    static constexpr char nNei = 2;
    static constexpr char nei_start = 3;
};

}  // namespace traccc::device

namespace traccc {

struct gbts_graph_building_params {

    // edge making cuts
    float min_delta_phi = 0.015f;
    float dphi_coeff = 2.2e-4f;
    float min_delta_phi_low_dr = 0.002f;
    float dphi_coeff_low_dr = 4.33e-4f;

    float minDeltaRadius = 2.0f;

    float min_z0 = -160.0f;
    float max_z0 = 160.0f;
    float maxOuterRadius = 350.0f;
    // how to get ROI dzdr
    float cut_zMinU = min_z0 - maxOuterRadius * 45.0f;
    float cut_zMaxU = max_z0 + maxOuterRadius * 45.0f;

    float max_Kappa = 3.75e-4f;
    float low_Kappa_d0 = 0.00f;
    float high_Kappa_d0 = 0.0f;

    // tau prediction cut
    float tMin_slope = 6.7f;
    float offset = 0.2f;
    float tMax_min = 1.6f;
    float tMax_correction = 0.15f;
    float tMax_slope = 6.1f;

    float type1_max_width = 0.2f;

    // edge matching cuts
    float cut_dphi_max = 0.012f;
    float cut_dcurv_max = 0.001f;
    float cut_tau_ratio_max = 0.01f;
};

struct gbts_seed_extraction_params {
    // for 900 MeV track at eta=0
    float sigmaMS = 0.016f;
    // 2.5% per layer
    float radLen = 0.025f;

    float sigma_x = 0.08f;
    float sigma_y = 0.25f;

    float weight_x = 0.5f;
    float weight_y = 0.5f;

    float maxDChi2_x = 5.0f;
    float maxDChi2_y = 6.0f;
    // controls if seeds of shorter lengths
    // can win bidding against longer seeds
    float add_hit = 14.0f;
    // seed quality is an int scaled up from a float
    // max qual = add_hit*max_length*qual_scale
    float qual_scale =
        0.01f * static_cast<float>(INT_MAX) /
        (add_hit *
         static_cast<float>(traccc::device::gbts_consts::max_cca_iter));

    float inv_max_curvature = 900.0f;
    float max_z0 = 160.0f;
};

struct gbts_seed_ambi_params {
    // sample multiple triplets when forming seeds to hedge against outliers
    bool use_dropout = true;
    // these curvatures are in 1/m
    float dropout_dcurv_m = 0.007f;
    float force_dropout_max_curv_m = 0.03f;
    float best_hit_frac = 0.49f;
    float tight_bid_cot_threshold = 1.0f;
};

struct gbts_seedfinder_config {
    bool setLinkingScheme(
        const std::vector<std::pair<int, std::vector<int>>>& binTables,
        const device::gbts_layerInfo layerInfo,
        std::vector<std::pair<uint64_t, short>>& detrayGeoIDBinning,
        float minPt, std::unique_ptr<const traccc::Logger> logger);

    // layer linking and geometry
    std::vector<std::pair<unsigned int, unsigned int>> binTables{};
    traccc::device::gbts_layerInfo layerInfo{};
    unsigned int nLayers = 0;

    std::vector<short> volumeToLayerMap{};
    std::vector<std::array<unsigned int, 2>> surfaceToLayerMap{};

    // tuned for 900 MeV pT cut and scaled by input minPt
    gbts_graph_building_params graph_building_params{};
    gbts_seed_extraction_params seed_extraction_params{};
    gbts_seed_ambi_params seed_ambi_params{};

    // node making bin counts
    unsigned int n_eta_bins = 0;  // calculated from input layerInfo
    unsigned int n_phi_bins = 128;
    // graph making maxiums
    unsigned int max_num_neighbours = 10;
    // graph extraction cuts
    int minLevel = 3;  // equivlent to a cut of #seed edges or #spacepoints-1

    // maxium number of edges to be created per node(spacepoint)
    unsigned int max_edges_factor = 10;
};

}  // namespace traccc
