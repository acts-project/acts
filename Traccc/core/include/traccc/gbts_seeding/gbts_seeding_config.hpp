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
  std::vector<std::pair<unsigned int, unsigned int>> info;
  // minEta and deltaEta
  std::vector<std::pair<float, float>> geo;

  void reserve(unsigned int n) {
    type.reserve(n);
    info.reserve(n);
    geo.reserve(n);
  }

  void addLayer(char layerType, unsigned int firstBin, unsigned int nBins,
                float minEta, float etaBinWidth) {
    type.push_back(layerType);
    info.push_back(std::make_pair(firstBin, nBins));
    geo.push_back(std::make_pair(minEta, etaBinWidth));
  }
};

// Named indices into the flat device counter buffer, mirroring the layout in
// traccc/gbts_changes. One memset zeros all of them.
enum gbts_counter : unsigned int {
  nEdges,           // edges created by gbts_make_graph_edges
  nConnections,     // edge-to-edge connections from gbts_match_graph_edges
  nConnectedEdges,  // edges kept after gbts_reindex_edges
  nEdgesLeft,       // edges remaining for CCA (kept for reference parity)
  nPaths,           // total paths reachable from any terminus edge
  nTerminusEdges,   // #terminus edges; then reused as path-store write cursor
  nProps,           // seed proposals from gbts_fit_segments
  nRejected,        // rejected seed proposals
  nCounters         // total number of counters
};

struct gbts_consts {
  // CCA max iterations -> maximum seed length (in edges).
  static constexpr unsigned short max_cca_iter = 15;
  // shared memory allocation sizes (element counts per block).
  // Which is used in the fill_path_store, and store 2 unsigned int.
  static constexpr unsigned short live_path_buffer = 1024;
  static constexpr unsigned short node_buffer_length = 128;

  // Per-edge offsets into the row-major output graph
  // (each edge occupies edge_size = 2 + 1 + max_num_neighbours ints).
  static constexpr unsigned char node1 = 0;
  static constexpr unsigned char node2 = 1;
  static constexpr unsigned char nNei = 2;
  static constexpr unsigned char nei_start = 3;
};

}  // namespace traccc::device

namespace traccc {

// Tau-prediction cuts for device::gbts_sort_nodes.
struct gbts_sort_nodes_params {
  // Slope of the lower-bound tau line: min_tau = tMin_slope * (w - offset).
  float tMin_slope = 6.7f;
  // Minimum cluster width: offset = w_min.
  float offset = 0.2f;
  // Asymptotic lower bound on the upper-tau line (w -> large).
  float tMax_min = 1.6f;
  // Inverse-width correction term on the upper-tau line:
  // tMax_correction/(w+offset).
  float tMax_correction = 0.15f;
  // Slope of the upper-tau line:
  //   max_tau = tMax_min + tMax_correction / (w + offset) + tMax_slope * (w -
  //   offset).
  float tMax_slope = 6.1f;
  // |tau| acceptance for nodes without a usable cluster width.
  //   sinh(|eta_max|) = maxTau   ->   sinh(4.3) ~ 36
  float maxTau = 36.0f;
  // Opt-in: use a tau lookup table instead of the GPU-friendly linear
  // formula. The LUT is laid out as [w_bin_edge, min_tau_0, max_tau_0,
  // min_tau_1, max_tau_1].
  // TODO: Test this with the LUT on CPU
  bool useTauLUT = false;
  // Inverse cluster-width bin size used to index the LUT.
  float tau_lut_inv_bin = 0.0f;
  // Number of float entries in the LUT (bounds the index).
  unsigned int tauLutSize = 0;
};

// Geometric / kinematic edge-making cuts for device::gbts_make_graph_edges.
struct gbts_make_graph_edges_params {
  // Two nodes must be radially separated to form an edge:
  // dr >= minDeltaRadius (mm)
  float minDeltaRadius = 2.0f;
  // z estimate at the beamline: z0 = z1 - r1*tau.
  // required min_z0 <= z0 <= max_z0.
  // +/-160 mm = chosen luminous-region.
  float min_z0 = -160.0f;
  float max_z0 = 160.0f;
  // Outer radius (mm) to which the edge is extrapolated for the ROI z cut.
  // Geometry input (~ outer pixel radius).
  float maxOuterRadius = 350.0f;
  // ROI band for zouter = z0 + maxOuterRadius*tau, required in [cut_zMinU,
  // cut_zMaxU]:  cut_zMin/MaxU = -/+ (|z0|_max + maxOuterRadius * tau_roi).
  float tau_roi = 45.0f;
  float cut_zMinU = min_z0 - maxOuterRadius * tau_roi;
  float cut_zMaxU = max_z0 + maxOuterRadius * tau_roi;
  // Maximum edge curvature = minimum pT:
  //   max_Kappa = curv_max = kappa/2 = c*q*B/(2*pT)
  //   = 0.299792458*2/(2*0.9) m^-1 = 0.333 m^-1 = 3.33e-4 mm^-1 (2 T, 0.9
  //   GeV).
  float max_Kappa = 3.75e-4f;
  // Max transverse impact parameter, per curvature regime:
  //   d0 = r1*r2*(|curv| - max_Kappa) <= {low,high}_Kappa_d0  (mm).
  // 0 = prompt tracks only (no d0 allowance beyond the pT cut).
  float low_Kappa_d0 = 0.00f;
  float high_Kappa_d0 = 0.0f;
};

// Pair-matching cuts for device::gbts_match_graph_edges.
struct gbts_match_graph_edges_params {
  // |dphi_1 - dphi_2| <= cut_dphi_max rad.
  float cut_dphi_max = 0.012f;
  // |curv_1 - curv_2| <= cut_dcurv_max  (curv = dphi/dr).
  float cut_dcurv_max = 0.001f;
  // |tau_2/tau_1 - 1| <= cut_tau_ratio_max (~1%).
  float cut_tau_ratio_max = 0.01f;
};

// Host-side dphi window used to compute bin_pair_dphi before launching
// device::gbts_make_graph_edges.
struct gbts_dphi_window_params {
  // deltaPhi = min_delta_phi + dphi_coeff * maxDeltaR, where maxDeltaR is the
  // maximum radial separation of the pair of nodes.
  float min_delta_phi = 0.015f;
  float dphi_coeff = 2.2e-4f;
  // deltaPhi window used if deltaR < low_dr_threshold (mm)
  float min_delta_phi_low_dr = 0.002f;
  float dphi_coeff_low_dr = 4.33e-4f;
  // delta-R (mm) below which the "low dr" window is used.
  float low_dr_threshold = 60.0f;
};

// Kalman-filter cuts for device::gbts_fit_segments.
struct gbts_fit_segments_params {
  // Per-layer multiple-scattering angle:
  // sigmaMS = E_s / pT = 14.1/900 = 0.0156  (900 MeV, eta=0).
  // 14.1 MeV is the Highland constant
  float sigmaMS = 0.016f;

  // Material per layer in radiation lengths: radLen = x/X0 = 2.5%.
  // Enters the MS covariance as radLen*sigmaMS^2.
  float radLen = 0.025f;

  // Measurement resolutions in the x-y (X) and r-z (Y) coordinates:
  float sigma_x = 0.08f;
  float sigma_y = 0.25f;

  // Relative weights of the x-y (X) / r-z (Y) chi2 terms in the seed quality.
  float weight_x = 0.5f;
  float weight_y = 0.5f;

  // Per-hit chi2 acceptance in x-y (X) / r-z (Y).
  float maxDChi2_x = 5.0f;
  float maxDChi2_y = 6.0f;

  // controls if seeds of shorter lengths
  // can win bidding against longer seeds.
  float add_hit = 14.0f;

  // seed quality is an int scaled up from a float
  // Exact int-scaling so the longest seed maps to ~1% of INT_MAX:
  float qual_scale =
      0.01f * static_cast<float>(INT_MAX) /
      (add_hit * static_cast<float>(traccc::device::gbts_consts::max_cca_iter));

  // Minimum-pT gate in the fit: reject if |X2| * inv_max_curvature > 1
  // inv_max_curvature = 1/curv_max = ~pT[MeV].
  float inv_max_curvature = 900.0f;

  // max_z0 is used from the graph_making to insure concistency
};

// Seed ambiguity / dropout parameters for device::gbts_convert_seeds
struct gbts_convert_seeds_params {
  // sample multiple triplets when forming seeds to hedge against outliers.
  bool use_dropout = true;
  // Curvature thresholds (1/m) for the dropout logic.
  // dcurv = kappa = c*q*B/pT, so 0.007 = ~86GeV
  float dropout_dcurv_m = 0.007f;
  // 0.03 = ~20 GeV,
  float force_dropout_max_curv_m = 0.03f;
  // Fraction of shared hits above which a seed loses a bid (~1/2). Tuning.
  float best_hit_frac = 0.49f;
  // Region switch for "tight" bidding:
  // cot(theta) = sinh(eta) = tight_bid_cot_threshold
  //   -> sinh(0.88) ~ 1.0, so tracks with |eta| < 0.88 get tighter
  //   hit-sharing cuts.
  float tight_bid_cot_threshold = 1.0f;
};

// SP counting cuts for device::gbts_count_spacepoints_by_layer
struct gbts_count_spacepoints_by_layer_params {
  // Maximum cluster width allowed on "type 1" layers.
  float type1_max_width = 0.2f;
  // If true, apply the cluster-width / tau cut at SP-counting time.
  bool doTauCut = true;
};

// Main configuration struct for the GBTS seeding algorithm.
struct gbts_seedfinder_config {
  bool setLinkingScheme(
      const std::vector<std::pair<unsigned int, std::vector<unsigned int>>>&
          binTables,
      const device::gbts_layerInfo layerInfo,
      std::vector<std::pair<uint64_t, int16_t>>& detrayGeoIDBinning,
      const float minPt, std::unique_ptr<const traccc::Logger> logger);

  // layer linking and geometry
  std::vector<std::pair<unsigned int, unsigned int>> binTables{};
  traccc::device::gbts_layerInfo layerInfo{};
  unsigned int nLayers = 0;

  std::vector<int16_t> volumeToLayerMap{};
  std::vector<std::pair<unsigned int, unsigned int>> surfaceToLayerMap{};

  // Per kernel structs
  traccc::gbts_sort_nodes_params gbts_sort_nodes_params{};
  traccc::gbts_make_graph_edges_params gbts_make_graph_edges_params{};
  traccc::gbts_match_graph_edges_params gbts_match_graph_edges_params{};
  traccc::gbts_dphi_window_params gbts_dphi_window_params{};
  traccc::gbts_count_spacepoints_by_layer_params
      gbts_count_spacepoints_by_layer_params{};
  traccc::gbts_fit_segments_params gbts_fit_segments_params{};
  traccc::gbts_convert_seeds_params gbts_convert_seeds_params{};

  // Optional tau lookup table consumed by device::gbts_sort_nodes when
  // gbts_sort_nodes_params.useTauLUT is set.
  // Layout: [w_bin_edge, min_tau_0, max_tau_0, min_tau_1, max_tau_1]
  std::vector<float> tau_lut{};

  // node making bin counts
  // calculated from input layerInfo (geometry)
  unsigned int n_eta_bins = 0;

  // Phi bin width
  unsigned int n_phi_bins = 128;

  // graph making maxiums: max neighbours kept per edge.
  unsigned int max_num_neighbours = 10;

  // graph extraction cuts: minimum number of edges as a CCA level,
  //   nSP_seed = minLevel + 1
  unsigned char minLevel = 3;

  // nMaxEdges = max_edges_factor * nNodes  (fixed buffer size).
  unsigned int max_edges_factor = 10;

  // number of seed-vs-edge bidding rounds during disambiguation.
  unsigned int edge_bidding_rounds = 5;
};

}  // namespace traccc
