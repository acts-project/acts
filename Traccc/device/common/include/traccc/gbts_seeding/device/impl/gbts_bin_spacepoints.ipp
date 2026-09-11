/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/math.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/device/concepts/thread_id.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/device_atomic_ref.hpp>

// Detray include(s).
#include <detray/geometry/identifier.hpp>

// System include(s).
#include <array>
#include <climits>
#include <utility>

namespace traccc::device {

/// Order-preserving float to uint32 transform.
TRACCC_HOST_DEVICE inline unsigned int float_ordered_bits(const float phi) {
  const unsigned int u = std::bit_cast<unsigned int>(phi);
  return ((u & 0x80000000u) != 0u) ? ~u : (u | 0x80000000u);
}

/// Node sort key layout, from the least significant bit up: the low bits
/// of the spacepoint index (a tie-break between nodes of identical eta bin
/// and phi only, not a lookup).
inline constexpr unsigned int gbts_sort_key_index_bits = 12u;
inline constexpr unsigned int gbts_sort_key_phi_bits = 32u;
inline constexpr unsigned int gbts_sort_key_eta_bits =
    64u - gbts_sort_key_phi_bits - gbts_sort_key_index_bits;
inline constexpr unsigned int gbts_sort_key_phi_shift =
    gbts_sort_key_index_bits;
inline constexpr unsigned int gbts_sort_key_eta_shift =
    gbts_sort_key_index_bits + gbts_sort_key_phi_bits;
/// Mask selecting the spacepoint-index bits of a node sort key
inline constexpr unsigned long long int gbts_sort_key_index_mask =
    (1ull << gbts_sort_key_index_bits) - 1ull;
/// Largest number of eta bins the key's eta field can hold
inline constexpr unsigned int gbts_sort_key_max_eta_bins =
    1u << gbts_sort_key_eta_bits;

template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_bin_spacepoints(
    const thread_id_t& thread_id, const gbts_bin_spacepoints_payload& payload) {
  const traccc::edm::spacepoint_collection::const_device spacepoints(
      payload.spacepoints);
  const edm::measurement_collection::const_device measurements(
      payload.measurements);
  const vecmem::device_vector<const short> volumeToLayerMap(
      payload.volumeToLayerMap);
  const vecmem::device_vector<const std::pair<unsigned int, unsigned int>>
      surfaceToLayerMap(payload.surfaceToLayerMap);
  const vecmem::device_vector<const char> layerType(payload.layerType);
  const vecmem::device_vector<const std::pair<unsigned int, unsigned int>>
      d_layer_info(payload.layer_info);
  const vecmem::device_vector<const std::pair<float, float>> d_layer_geo(
      payload.layer_geo);

  vecmem::device_vector<float4> reducedSP(payload.reducedSP);
  vecmem::device_vector<unsigned int> d_eta_node_counter(
      payload.eta_node_counter);
  vecmem::device_vector<unsigned long long int> d_sort_keys(payload.sort_keys);
  vecmem::device_vector<unsigned int> d_sort_values(payload.sort_values);

  const unsigned int globalIdx = thread_id.getGlobalThreadIdX();
  const unsigned int blockDimX = thread_id.getBlockDimX();
  const unsigned int gridDimX = thread_id.getGridDimX();

  for (unsigned int globalIndex = globalIdx; globalIndex < payload.nSp;
       globalIndex += blockDimX * gridDimX) {
    // --- Stage 1: layer assignment -----------
    const auto spacepoint = spacepoints.at(globalIndex);
    const auto measurement = measurements.at(spacepoint.measurement_index_1());

    const detray::geometry::identifier geo_id = measurement.surface_link();
    const unsigned int volume = geo_id.volume();
    const short begin_or_bin =
        (volume < payload.volumeMapSize) ? volumeToLayerMap[volume] : SHRT_MAX;

    if (begin_or_bin == SHRT_MAX) {
      reducedSP[globalIndex].w = -CHAR_MAX - 1;
      continue;
    }
    unsigned int layerIdx = 0u;
    if (begin_or_bin < 0) {
      const unsigned int surface_index =
          static_cast<unsigned int>(geo_id.index());

      for (unsigned int surface =
               static_cast<unsigned int>(-1 * (begin_or_bin + 1));
           surface < payload.surfaceMapSize; surface++) {
        const std::pair<unsigned int, unsigned int> surfaceBinPair =
            surfaceToLayerMap[surface];
        if (surfaceBinPair.first == surface_index) {
          layerIdx = surfaceBinPair.second;
          break;
        }
      }
    } else {
      layerIdx = static_cast<unsigned int>(begin_or_bin);
    }
    float cluster_diameter = measurement.diameter();
    const int type = static_cast<int>(layerType[layerIdx]);
    if (type == 1 &&
        cluster_diameter >
            payload.gbts_count_spacepoints_by_layer_params.type1_max_width) {
      reducedSP[globalIndex].w = -CHAR_MAX - 1;
      continue;
    }
    cluster_diameter =
        (payload.gbts_count_spacepoints_by_layer_params.doTauCut && type != 0)
            ? static_cast<float>(-1 * type)
            : cluster_diameter;

    const std::array<float, 3u> pos = spacepoint.global();
    reducedSP[globalIndex] = float4{pos[0], pos[1], pos[2], cluster_diameter};
    // global x, y, z, and cluster diameter

    // --- Stage 2: node_eta_binning -----------
    const std::pair<unsigned int, unsigned int> layerInfo =
        d_layer_info[layerIdx];
    const unsigned int bin0 = layerInfo.first;
    const unsigned int num_eta_bins = layerInfo.second;
    unsigned int eta_index;
    if (num_eta_bins == 1u) {
      eta_index = bin0;
    } else {
      const std::pair<float, float> layerGeo = d_layer_geo[layerIdx];
      const float min_eta = layerGeo.first;
      const float eta_bin_width = layerGeo.second;
      const float r = math::sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
      const float t1 = pos[2] / r;
      const float eta = -math::log(math::sqrt(1.0f + t1 * t1) - t1);
      const unsigned int binIdx = static_cast<unsigned int>(
          math::max(0.0f, math::min((eta - min_eta) / eta_bin_width,
                                    static_cast<float>(num_eta_bins - 1u))));
      eta_index = bin0 + binIdx;
    }
    vecmem::device_atomic_ref<unsigned int>(d_eta_node_counter[eta_index])
        .fetch_add(1u);

    // --- Stage 3: node_sort_key_write -----------
    // Here we concatenate the eta bin, phi bits, and the spacepoint index
    // into a single 64-bit integer. Which is then used to sort the nodes in the
    // next stage.
    const float Phi = math::atan2(pos[1], pos[0]);
    const unsigned int slot = vecmem::device_atomic_ref<unsigned int>(
                                  d_eta_node_counter[payload.nEtaBins])
                                  .fetch_add(1u);
    d_sort_keys[slot] =
        (static_cast<unsigned long long int>(eta_index)
         << gbts_sort_key_eta_shift) |
        (static_cast<unsigned long long int>(float_ordered_bits(Phi))
         << gbts_sort_key_phi_shift) |
        (static_cast<unsigned long long int>(globalIndex) &
         gbts_sort_key_index_mask);
    d_sort_values[slot] = globalIndex;
  }
}

}  // namespace traccc::device
