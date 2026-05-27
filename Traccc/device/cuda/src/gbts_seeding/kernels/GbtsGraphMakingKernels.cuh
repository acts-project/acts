/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// CUDA include(s)
#include <cuda.h>
#include <cuda_fp16.h>
#include <cuda_runtime.h>
#include <math_constants.h>
#include <vector_functions.h>

// Project include(s)
#include "traccc/cuda/gbts_seeding/gbts_seeding_algorithm.hpp"

namespace traccc::cuda::kernels {

struct __align__(8) half4 {
    __half x, y, z, w;
};

inline __device__ __host__ half4 make_half4(const __half x, const __half y,
                                            const __half z, const __half w) {
    half4 t;
    t.x = x;
    t.y = y;
    t.z = z;
    t.w = w;
    return t;
}

__global__ static void graphEdgeMakingKernel(
    const uint4* d_bin_pair_views, const float* d_bin_pair_dphi,
    const float* d_node_params,
    const gbts_graph_building_params* d_graph_building_params,
    unsigned int* d_counters, int2* d_edge_nodes, half4* d_edge_params,
    int* d_num_outgoing_edges, const unsigned int nMaxEdges,
    const unsigned int nPhiBins) {

    __shared__ unsigned int begin_bin1;
    __shared__ unsigned int begin_bin2;
    __shared__ unsigned int num_nodes1;
    __shared__ unsigned int num_nodes2;
    __shared__ float deltaPhi;

    __shared__ float minDeltaRad;
    __shared__ float min_z0;
    __shared__ float max_z0;
    __shared__ float maxOuterRad;
    __shared__ float min_zU;
    __shared__ float max_zU;
    __shared__ float max_kappa;
    __shared__ float low_Kappa_d0;
    __shared__ float high_Kappa_d0;

    __shared__ float tau_min[traccc::device::gbts_consts::node_buffer_length];
    __shared__ float tau_max[traccc::device::gbts_consts::node_buffer_length];
    __shared__ float phi[traccc::device::gbts_consts::node_buffer_length];
    __shared__ float r[traccc::device::gbts_consts::node_buffer_length];
    __shared__ float z[traccc::device::gbts_consts::node_buffer_length];

    if (threadIdx.x == 0) {
        uint4 views = d_bin_pair_views[blockIdx.x];
        deltaPhi = d_bin_pair_dphi[blockIdx.x];

        begin_bin1 = views.x;
        begin_bin2 = views.z;
        num_nodes1 = views.y - begin_bin1;
        num_nodes2 = views.w - begin_bin2;

        minDeltaRad = d_graph_building_params->minDeltaRadius;
        min_z0 = d_graph_building_params->min_z0;
        max_z0 = d_graph_building_params->max_z0;
        maxOuterRad = d_graph_building_params->maxOuterRadius;
        min_zU = d_graph_building_params->cut_zMinU;
        max_zU = d_graph_building_params->cut_zMaxU;
        max_kappa = d_graph_building_params->max_Kappa;
        low_Kappa_d0 = d_graph_building_params->low_Kappa_d0;
        high_Kappa_d0 = d_graph_building_params->high_Kappa_d0;
    }

    __syncthreads();
    for (int idx = threadIdx.x; idx < num_nodes1; idx += blockDim.x) {
        // loading a chunk of nodes1 into shared mem buffers
        int offset = 5 * (idx + begin_bin1);
        tau_min[idx] = d_node_params[offset];
        tau_max[idx] = d_node_params[offset + 1];
        phi[idx] = d_node_params[offset + 2];
        r[idx] = d_node_params[offset + 3];
        z[idx] = d_node_params[offset + 4];
    }

    __syncthreads();

    int last_n1 = 0;  // initial value for the sliding window

    float phi_bin_width = 2.0f * CUDART_PI_F / nPhiBins;

    for (int n2Idx = threadIdx.x; n2Idx < num_nodes2; n2Idx += blockDim.x) {

        int n1Idx = last_n1;

        int globalIdx2 = begin_bin2 + n2Idx;
        int o2 = 5 * globalIdx2;

        float phi2 = d_node_params[2 + o2];

        float min_phi1 = phi2 - deltaPhi;
        float max_phi1 = phi2 + deltaPhi;

        if (min_phi1 < -CUDART_PI_F) {
            min_phi1 += 2.0f * CUDART_PI_F;
        }
        if (max_phi1 > CUDART_PI_F) {
            max_phi1 -= 2.0f * CUDART_PI_F;
        }
        bool boundary = max_phi1 < min_phi1;  // +/- pi wraparound

        // expand over nearest bin boundary
        max_phi1 += phi_bin_width;
        min_phi1 -= phi_bin_width;

        if (!boundary) {
            if (phi[0] > max_phi1) {
                continue;
            }
            if (phi[num_nodes1 - 1] < min_phi1) {
                // if bin1 can't be part of a wraparound
                // from a high-phi node skip it
                if (phi[0] > deltaPhi + phi_bin_width - CUDART_PI_F) {
                    break;
                }
                continue;
            }
        } else {
            if (phi[0] < max_phi1) {
                // if not to large for lower wraparound don't skip it
                n1Idx = 0;
            } else if (phi[num_nodes1 - 1] < min_phi1) {
                continue;
            }
        }

        float tau_min2 = d_node_params[o2];
        float tau_max2 = d_node_params[1 + o2];
        float r2 = d_node_params[3 + o2];
        float z2 = d_node_params[4 + o2];

        for (; n1Idx < num_nodes1; n1Idx++) {
            float phi1 = phi[n1Idx];

            if (!boundary) {
                if (phi1 > max_phi1) {
                    break;
                }
                if (phi1 < min_phi1) {
                    continue;
                }
                last_n1 = n1Idx;
            } else {
                if (phi1 > max_phi1 && phi1 < min_phi1) {
                    // skip to high wraparound after the lower part is done
                    if (n1Idx < last_n1) {
                        n1Idx = last_n1 - 1;
                    }
                    continue;
                }
            }

            float r1 = r[n1Idx];
            float dr = r2 - r1;

            if (dr < minDeltaRad) {
                continue;
            }
            float z1 = z[n1Idx];
            float dz = z2 - z1;
            float tau = dz / dr;
            float ftau = fabsf(tau);

            if (ftau > 36.0f) {
                continue;  // detector acceptance
            }
            if ((ftau < tau_min2) || (ftau > tau_max2)) {
                continue;
            }
            if ((ftau < tau_min[n1Idx]) || (ftau > tau_max[n1Idx])) {
                continue;
            }
            // RZ doublet filter cuts
            float z0 = z1 - r1 * tau;
            if ((z0 < min_z0) || (z0 > max_z0)) {
                continue;
            }
            float zouter = z0 + maxOuterRad * tau;

            if (zouter < min_zU || zouter > max_zU) {
                continue;
            }
            float dphi = phi2 - phi1;
            if (boundary) {
                if (dphi < -CUDART_PI_F)
                    dphi += 2.0f * CUDART_PI_F;
                else if (dphi > CUDART_PI_F)
                    dphi -= 2.0f * CUDART_PI_F;
            }

            // needed for sliding phi window consistancy
            if (fabsf(dphi) > deltaPhi) {
                continue;
            }
            float curv = dphi / dr;
            float d0_for_max_curv = r1 * r2 * (fabsf(curv) - max_kappa);
            float d0_max = (ftau < 4.0f) ? low_Kappa_d0 : high_Kappa_d0;
            if (d0_for_max_curv > d0_max) {
                continue;
            }
            unsigned int nEdges = atomicAdd(&d_counters[0], 1);
            if (nEdges < nMaxEdges) {
                __half exp_eta = __float2half(sqrtf(1 + tau * tau) - tau);
                // edge linking order is inside->out
                atomicAdd(&d_num_outgoing_edges[begin_bin1 + n1Idx], 1);

                d_edge_nodes[nEdges] =
                    make_int2(globalIdx2, begin_bin1 + n1Idx);

                d_edge_params[nEdges] = make_half4(
                    exp_eta, __float2half(curv), __float2half(phi2 + curv * r2),
                    __float2half(phi1 + curv * r1));
            }
        }
    }
}

__global__ static void graphEdgeLinkingKernel(const int2* d_edge_nodes,
                                              int* d_edge_links,
                                              int* d_num_outgoing_edges,
                                              const unsigned int nEdges) {

    int edge_idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (edge_idx >= nEdges) {
        return;
    }
    int sharedNode = d_edge_nodes[edge_idx].y;

    // this converts num_outgoing_edges to the start postion for each node in
    // d_edge_links
    int pos = atomicSub(&d_num_outgoing_edges[sharedNode], 1);
    // provides views of edges leaving the sharedNode for linking
    d_edge_links[pos - 1] = edge_idx;
}

__global__ static void graphEdgeMatchingKernel(
    const gbts_graph_building_params* d_graph_building_params,
    const half4* d_edge_params, const int2* d_edge_nodes,
    const int* d_num_outgoing_edges, const int* d_edge_links,
    unsigned char* d_num_neighbours, int* d_neighbours, int* d_reIndexer,
    unsigned int* d_counters, const unsigned int nEdges,
    const unsigned int nMaxNei) {
    __shared__ __half cut_dphi_max;
    __shared__ __half cut_dcurv_max;
    __shared__ __half cut_tau_ratio_max;
    __shared__ __half PI_h;
    __shared__ __half PI_2_h;
    __shared__ __half ONE_h;
    if (threadIdx.x == 0) {
        cut_dphi_max = __float2half(d_graph_building_params->cut_dphi_max);
        cut_dcurv_max = __float2half(d_graph_building_params->cut_dcurv_max);
        cut_tau_ratio_max =
            __float2half(d_graph_building_params->cut_tau_ratio_max);

        PI_h = __float2half(CUDART_PI_F);
        PI_2_h = __float2half(2 * CUDART_PI_F);
        ONE_h = __float2half(1.0f);
    }
    __syncthreads();

    int edge1_idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (edge1_idx >= nEdges) {
        return;
    }

    int sharedNode = d_edge_nodes[edge1_idx].x;

    int link_begin = d_num_outgoing_edges[sharedNode];
    // the number of edges leaving the sharedNode
    int nLinks = d_num_outgoing_edges[sharedNode + 1] - link_begin;
    if (nLinks == 0) {
        return;
    }
    half4 params1 = d_edge_params[edge1_idx];  // [exp_eta, curv, Phi1, Phi2]

    __half uat_2 = ONE_h / params1.x;
    __half Phi2 = params1.z;
    __half curv2 = params1.y;

    int nei_pos = nMaxNei * edge1_idx;

    unsigned char num_nei = 0;

    for (int k = 0; k < nLinks; k++) {  // loop over potential neighbours

        if (num_nei >= nMaxNei) {
            break;
        }
        int edge2_idx = d_edge_links[link_begin + k];

        half4 params2 = d_edge_params[edge2_idx];

        __half tau_ratio = params2.x * uat_2 - ONE_h;

        if (__habs(tau_ratio) > cut_tau_ratio_max) {  // bad match
            continue;
        }

        __half dPhi = Phi2 - params2.w;  // Phi2

        if (dPhi < -PI_h) {
            dPhi += PI_2_h;
        } else if (dPhi > PI_h) {
            dPhi -= PI_2_h;
        }
        if (__habs(dPhi) > cut_dphi_max) {
            continue;
        }

        __half dcurv = curv2 - params2.y;

        if (__habs(dcurv) > cut_dcurv_max) {
            continue;
        }

        d_neighbours[nei_pos + num_nei] = edge2_idx;
        d_reIndexer[edge2_idx] = 1;

        ++num_nei;
    }

    d_num_neighbours[edge1_idx] = num_nei;

    if (num_nei != 0) {
        d_reIndexer[edge1_idx] = 1;
        atomicAdd(&d_counters[1], num_nei);
    }
}

__global__ void edgeReIndexingKernel(int* d_reIndexer, unsigned int* d_counters,
                                     const unsigned int nEdges) {

    // each thread gets an edge

    int edge_idx = threadIdx.x + blockIdx.x * blockDim.x;

    if (edge_idx >= nEdges) {
        return;
    }
    if (d_reIndexer[edge_idx] == -1) {
        return;
    }
    d_reIndexer[edge_idx] = atomicAdd(&d_counters[2], 1);
}

__global__ static void graphCompressionKernel(
    const int* d_orig_node_index, const int2* d_edge_nodes,
    const unsigned char* d_num_neighbours, const int* d_neighbours,
    const int* d_reIndexer, int* d_output_graph,
    const unsigned int nEdgesPerBlock, const unsigned int nEdges,
    const unsigned int nMaxNei) {

    int begin_edge = blockIdx.x * nEdgesPerBlock;
    int edge_size = 2 + 1 + nMaxNei;

    for (int idx = threadIdx.x + begin_edge; idx < begin_edge + nEdgesPerBlock;
         idx += blockDim.x) {

        if (idx >= nEdges) {
            continue;
        }
        int newIdx = d_reIndexer[idx];
        if (newIdx == -1) {
            continue;
        }
        int pos = edge_size * newIdx;
        int2 edge_nodes = d_edge_nodes[idx];
        int node1_idx = d_orig_node_index[edge_nodes.x];
        d_output_graph[pos + traccc::device::gbts_consts::node1] = node1_idx;
        int node2_idx = d_orig_node_index[edge_nodes.y];
        d_output_graph[pos + traccc::device::gbts_consts::node2] = node2_idx;

        unsigned char nNei = d_num_neighbours[idx];
        d_output_graph[pos + traccc::device::gbts_consts::nNei] = nNei;
        int nei_pos = nMaxNei * idx;
        for (int k = 0; k < nNei; k++) {
            d_output_graph[pos + traccc::device::gbts_consts::nei_start + k] =
                d_reIndexer[d_neighbours[nei_pos + k]];
        }
    }
}

}  // namespace traccc::cuda::kernels
