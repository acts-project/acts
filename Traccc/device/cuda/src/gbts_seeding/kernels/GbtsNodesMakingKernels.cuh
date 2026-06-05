/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// CUDA include(s)
#include <cuda.h>
#include <cuda_runtime.h>
#include <math_constants.h>
#include <vector_functions.h>

// Project include(s)
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"

// Detray include(s).
#include <detray/geometry/identifier.hpp>

namespace traccc::cuda::kernels {

__global__ void count_sp_by_layer(
    const traccc::edm::spacepoint_collection::const_view spacepoints_view,
    const edm::measurement_collection::const_view measurements_view,
    const short* volumeToLayerMap, const uint2* surfaceToLayerMap,
    const char* d_layerType, float4* reducedSP, int* d_layerCounts,
    short* spacepointsLayer, const float type1_max_width,
    const unsigned int nSp, const long unsigned int volumeMapSize,
    const long unsigned int surfaceMapSize, bool doTauCut = true) {

    const edm::measurement_collection::const_device measurements(
        measurements_view);
    const traccc::edm::spacepoint_collection::const_device spacepoints(
        spacepoints_view);

    for (int spIdx = threadIdx.x + blockDim.x * blockIdx.x; spIdx < nSp;
         spIdx += blockDim.x * gridDim.x) {
        // get the layer of the spacepoint
        const traccc::edm::spacepoint_collection::const_device::const_proxy_type
            spacepoint = spacepoints.at(spIdx);
        const auto measurement =
            measurements.at(spacepoint.measurement_index_1());

        detray::geometry::identifier geo_id = measurement.surface_link();

        // some volume_ids map one to one with layer others need searching
        if (geo_id.volume() > volumeMapSize) {
            reducedSP[spIdx].w = -CHAR_MAX - 1;
            continue;  // unconfigured volume
        }
        short begin_or_bin = volumeToLayerMap[geo_id.volume()];
        if (begin_or_bin == SHRT_MAX) {
            reducedSP[spIdx].w = -CHAR_MAX - 1;
            continue;  // unconfigured volume
        }
        unsigned int layerIdx;
        if (begin_or_bin < 0) {
            unsigned int surface_index =
                static_cast<unsigned int>(geo_id.index());

            for (unsigned int surface = -1 * (begin_or_bin + 1);
                 surface < surfaceMapSize; surface++) {

                uint2 surfaceBinPair = surfaceToLayerMap[surface];
                if (surfaceBinPair.x == surface_index) {
                    layerIdx = surfaceBinPair.y;
                    break;
                }
            }
        } else {
            layerIdx = static_cast<unsigned int>(begin_or_bin);
        }
        float cluster_diameter = measurement.diameter();
        int type = static_cast<int>(d_layerType[layerIdx]);
        if (type == 1 && cluster_diameter > type1_max_width) {
            //-ve cluster_diameter to skip cot(theta) prediction
            // large -ve to skip spacepoint entirely
            reducedSP[spIdx].w = -CHAR_MAX - 1;
            continue;
        }
        cluster_diameter =
            (doTauCut && type != 0) ? -1 * type : cluster_diameter;

        // count and store x,y,z,cw info
        atomicAdd(&d_layerCounts[layerIdx], 1);
        spacepointsLayer[spIdx] = layerIdx;
        const std::array<float, 3u> pos = spacepoint.global();
        reducedSP[spIdx] =
            make_float4(pos[0], pos[1], pos[2], cluster_diameter);
    }
}

// layerCounts is prefix sumed on CPU inbetween count_sp_by_layer and this
// kerenel
__global__ void bin_sp_by_layer(float4* sp_params, float4* reducedSP,
                                int* layerCounts, short* spacepointsLayer,
                                int* original_sp_idx, const unsigned int nSp) {

    for (int spIdx = threadIdx.x + blockDim.x * blockIdx.x; spIdx < nSp;
         spIdx += blockDim.x * gridDim.x) {

        float4 sp = reducedSP[spIdx];
        if (sp.w < -CHAR_MAX) {
            continue;
        }
        short layerIdx = spacepointsLayer[spIdx];
        unsigned int binedIdx = atomicSub(&layerCounts[layerIdx], 1) - 1;
        original_sp_idx[binedIdx] = spIdx;
        sp_params[binedIdx] = reducedSP[spIdx];
    }
}

__global__ void node_phi_binning_kernel(const float4* d_sp_params,
                                        int* d_node_phi_index,
                                        const unsigned int nNodesPerBlock,
                                        const unsigned int nNodes,
                                        const unsigned int nPhiBins) {

    int begin_node = blockIdx.x * nNodesPerBlock;

    float inv_phiSliceWidth = 1 / (2.0f * CUDART_PI_F / nPhiBins);

    for (int idx = threadIdx.x + begin_node; idx < begin_node + nNodesPerBlock;
         idx += blockDim.x) {

        if (idx >= nNodes) {
            continue;
        }
        float4 sp = d_sp_params[idx];

        float Phi = atan2f(sp.y, sp.x);

        int phiIdx = (Phi + CUDART_PI_F) * inv_phiSliceWidth;

        if (phiIdx >= nPhiBins) {
            phiIdx %= nPhiBins;
        } else if (phiIdx < 0) {
            phiIdx += nPhiBins;
            phiIdx %= nPhiBins;
        }
        d_node_phi_index[idx] = phiIdx;
    }
}

__global__ void node_eta_binning_kernel(const float4* d_sp_params,
                                        const int2* d_layer_info,
                                        const float2* d_layer_geo,
                                        int* d_node_eta_index,
                                        int* d_layerCounts,
                                        const unsigned int nLayers) {

    __shared__ int layer_begin;
    __shared__ int layer_end;
    __shared__ int num_eta_bins;
    __shared__ int bin0;
    __shared__ float min_eta;
    __shared__ float eta_bin_width;

    int layerIdx = blockIdx.x;

    if (threadIdx.x == 0) {
        int2 layerInfo = d_layer_info[layerIdx];
        bin0 = layerInfo.x;
        num_eta_bins = layerInfo.y;
        layer_begin = d_layerCounts[layerIdx];
        layer_end = d_layerCounts[layerIdx + 1];
        float2 layerGeo = d_layer_geo[layerIdx];
        min_eta = layerGeo.x;
        eta_bin_width = layerGeo.y;
    }

    __syncthreads();

    if (num_eta_bins == 1) {  // all nodes are in the same bin
        for (int idx = threadIdx.x + layer_begin; idx < layer_end;
             idx += blockDim.x) {

            d_node_eta_index[idx] = bin0;
        }
    } else {  // binIndex needs to be calculated first
        for (int idx = threadIdx.x + layer_begin; idx < layer_end;
             idx += blockDim.x) {

            float4 sp = d_sp_params[idx];

            float r = sqrtf(sp.x * sp.x + sp.y * sp.y);

            float t1 = sp.z / r;

            float eta = -logf(sqrtf(1 + t1 * t1) - t1);

            int binIdx = (int)((eta - min_eta) / eta_bin_width);

            if (binIdx < 0) {
                binIdx = 0;
            }
            if (binIdx >= num_eta_bins) {
                binIdx = num_eta_bins - 1;
            }
            d_node_eta_index[idx] = bin0 + binIdx;
        }
    }
}
// TO-DO fuse kernels?
__global__ void eta_phi_histo_kernel(const int* d_node_phi_index,
                                     const int* d_node_eta_index,
                                     int* d_eta_phi_histo,
                                     const unsigned int nNodesPerBlock,
                                     const unsigned int nNodes,
                                     const unsigned int nPhiBins) {

    int begin_node = blockIdx.x * nNodesPerBlock;

    for (int idx = threadIdx.x + begin_node; idx < begin_node + nNodesPerBlock;
         idx += blockDim.x) {

        if (idx >= nNodes) {
            continue;
        }
        int eta_index = d_node_eta_index[idx];

        int histo_bin = d_node_phi_index[idx] + nPhiBins * eta_index;
        atomicAdd(&d_eta_phi_histo[histo_bin], 1);
    }
}

__global__ void eta_phi_counting_kernel(const int* d_histo,
                                        int* d_eta_node_counter,
                                        int* d_phi_cusums,
                                        const unsigned nBinsPerBlock,
                                        const unsigned maxEtaBin,
                                        const unsigned nPhiBins) {

    int eta_bin_start = nBinsPerBlock * blockIdx.x;

    int eta_bin_idx = eta_bin_start + threadIdx.x;

    if (eta_bin_idx >= maxEtaBin) {
        return;
    }
    int offset = nPhiBins * eta_bin_idx;

    int sum = 0;

    for (int phiIdx = 0; phiIdx < nPhiBins; phiIdx++) {

        d_phi_cusums[offset + phiIdx] = sum;

        sum += d_histo[offset + phiIdx];
    }
    d_eta_node_counter[eta_bin_idx] = sum;
}

__global__ void eta_phi_prefix_sum_kernel(const int* d_eta_node_counter,
                                          int* d_phi_cusums,
                                          const unsigned int nBinsPerBlock,
                                          const unsigned int maxEtaBin,
                                          const unsigned int nPhiBins) {

    int eta_bin_start = nBinsPerBlock * blockIdx.x;

    int eta_bin_idx = eta_bin_start + threadIdx.x;

    if (eta_bin_idx >= maxEtaBin) {
        return;
    }
    if (eta_bin_idx == 0) {
        return;
    }
    int offset = nPhiBins * eta_bin_idx;

    int val0 = d_eta_node_counter[eta_bin_idx - 1];

    for (int phiIdx = 0; phiIdx < nPhiBins; phiIdx++) {
        d_phi_cusums[offset + phiIdx] += val0;
    }
}

__global__ void node_sorting_kernel(
    const float4* d_sp_params, const int* d_node_eta_index,
    const int* d_node_phi_index, int* d_phi_cusums, float* d_node_params,
    int* d_node_index, int* d_original_sp_idx,
    const gbts_graph_building_params* ap, const unsigned int nNodesPerBlock,
    const unsigned int nNodes, const unsigned int nPhiBins) {

    int begin_node = blockIdx.x * nNodesPerBlock;

    for (int idx = threadIdx.x + begin_node; idx < begin_node + nNodesPerBlock;
         idx += blockDim.x) {

        if (idx >= nNodes)
            continue;

        float4 sp = d_sp_params[idx];

        float Phi = atan2f(sp.y, sp.x);
        float r = sqrtf(sp.x * sp.x + sp.y * sp.y);
        float z = sp.z;

        float min_tau = -100.0;
        float max_tau = 100.0;

        if (sp.w > 0) {  // type 0 only
            // linear fit
            min_tau = ap->tMin_slope * (sp.w - ap->offset);
            // linear fit + correction for short clusters
            max_tau = ap->tMax_min + ap->tMax_correction / (sp.w + ap->offset) +
                      ap->tMax_slope * (sp.w - ap->offset);
        }

        int eta_index = d_node_eta_index[idx];
        int histo_bin = d_node_phi_index[idx] + nPhiBins * eta_index;

        int pos = atomicAdd(&d_phi_cusums[histo_bin], 1);

        int o = 5 * pos;

        d_node_params[o] = min_tau;
        d_node_params[o + 1] = max_tau;
        d_node_params[o + 2] = Phi;
        d_node_params[o + 3] = r;
        d_node_params[o + 4] = z;
        // keep the original index of the input spacepoint
        d_node_index[pos] = d_original_sp_idx[idx];
    }
}

__global__ void minmax_rad_kernel(const int2* d_eta_bin_views,
                                  const float* d_node_params,
                                  float2* d_bin_rads,
                                  const unsigned int nBinsPerBlock,
                                  const unsigned int maxEtaBin) {

    int eta_bin_start = nBinsPerBlock * blockIdx.x;

    int eta_bin_idx = eta_bin_start + threadIdx.x;

    if (eta_bin_idx >= maxEtaBin) {
        return;
    }

    int2 view = d_eta_bin_views[eta_bin_idx];
    int node_start = view.x;
    int node_end = view.y;
    if (node_start == node_end) {
        return;
    }
    float min_r = 1e8;
    float max_r = -1e8;

    for (int idx = node_start; idx < node_end; idx++) {
        float r = d_node_params[5 * idx + 3];
        if (r > max_r) {
            max_r = r;
        }
        if (r < min_r) {
            min_r = r;
        }
    }

    d_bin_rads[eta_bin_idx] = make_float2(min_r, max_r);
}

}  // namespace traccc::cuda::kernels
