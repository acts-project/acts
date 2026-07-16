/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/math.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/device/concepts/barrier.hpp"
#include "traccc/device/concepts/thread_id.hpp"
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"
#include "traccc/utils/trigonometric_helpers.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/device_atomic_ref.hpp>

namespace traccc::device {

namespace detail {

// Apply the RZ-doublet + curvature/d0 cuts to a candidate edge (node1->node2)
// and, if it passes, append it to the edge list. Single-use helper for
// gbts_make_graph_edges, kept inline here. Node params are float4 (tau_min,
// tau_max, r, z).
TRACCC_HOST_DEVICE inline void gbts_check_edge_candidate(
    const float4 node_params_1, const float4 node_params_2,
    vecmem::device_vector<uint2>& d_edge_nodes,
    vecmem::device_vector<float4>& d_edge_params,
    vecmem::device_vector<unsigned int>& d_num_outgoing_edges,
    unsigned int& nEdgesCounter, const unsigned int globalIdx2,
    const unsigned int begin_bin1, const unsigned int n1Idx, const float phi1,
    const float phi2, const float deltaPhi,
    const gbts_make_graph_edges_params& ap, const unsigned int nMaxEdges) {

    const float tau_min1 = node_params_1.x;
    const float tau_max1 = node_params_1.y;
    const float r1 = node_params_1.z;
    const float z1 = node_params_1.w;
    const float tau_min2 = node_params_2.x;
    const float tau_max2 = node_params_2.y;
    const float r2 = node_params_2.z;
    const float z2 = node_params_2.w;
    const float dr = r2 - r1;

    if (dr < ap.minDeltaRadius) {
        return;
    }
    const float dz = z2 - z1;
    const float tau = dz / dr;
    const float ftau = math::fabs(tau);

    if ((ftau < tau_min2) || (ftau > tau_max2)) {
        return;
    }
    if ((ftau < tau_min1) || (ftau > tau_max1)) {
        return;
    }
    // RZ doublet filter cuts
    const float z0 = z1 - r1 * tau;
    if ((z0 < ap.min_z0) || (z0 > ap.max_z0)) {
        return;
    }
    const float zouter = z0 + ap.maxOuterRadius * tau;
    if (zouter < ap.cut_zMinU || zouter > ap.cut_zMaxU) {
        return;
    }

    const float dphi = traccc::detail::wrap_phi(phi2 - phi1);
    if (math::fabs(dphi) > deltaPhi) {
        return;
    }
    const float curv = dphi / dr;
    const float d0_for_max_curv = r1 * r2 * (math::fabs(curv) - ap.max_Kappa);
    const float d0_max = (ftau < 4.0f) ? ap.low_Kappa_d0 : ap.high_Kappa_d0;
    if (d0_for_max_curv > d0_max) {
        return;
    }

    const unsigned int nEdges =
        vecmem::device_atomic_ref<unsigned int>(nEdgesCounter).fetch_add(1u);
    if (nEdges < nMaxEdges) {
        const float exp_eta = math::sqrt(1.0f + tau * tau) - tau;
        // edge linking order is inside->out
        vecmem::device_atomic_ref<unsigned int>(
            d_num_outgoing_edges[begin_bin1 + n1Idx])
            .fetch_add(1u);
        d_edge_nodes[nEdges] = uint2{globalIdx2, begin_bin1 + n1Idx};
        d_edge_params[nEdges] =
            float4{exp_eta, curv, phi2 + curv * r2, phi1 + curv * r1};
        // edge params: (exp(-eta), curvature, extrapolated phi at node1,
        //               extrapolated phi at node2)
    }
}

}  // namespace detail

template <concepts::thread_id1 thread_id_t, concepts::barrier barrier_t>
TRACCC_HOST_DEVICE inline void gbts_make_graph_edges(
    const thread_id_t& thread_id, const barrier_t& barrier,
    const gbts_make_graph_edges_payload& payload,
    const gbts_make_graph_edges_shared_payload& shared) {

    const unsigned int blockIndex = thread_id.getBlockIdX();
    const unsigned int threadIndex = thread_id.getLocalThreadIdX();
    const unsigned int blockSize = thread_id.getBlockDimX();

    const vecmem::device_vector<const unsigned int> d_bin_pair_views(
        payload.bin_pair_views);
    const vecmem::device_vector<const float> d_bin_pair_dphi(
        payload.bin_pair_dphi);
    const vecmem::device_vector<const float4> d_node_params(
        payload.node_params);
    const vecmem::device_vector<const float> d_node_phi(payload.node_phi);
    vecmem::device_vector<uint2> d_edge_nodes(payload.edge_nodes);
    vecmem::device_vector<float4> d_edge_params(payload.edge_params);
    vecmem::device_vector<unsigned int> d_num_outgoing_edges(
        payload.num_outgoing_edges);

    vecmem::device_vector<float> shared_phi(shared.phi);
    vecmem::device_vector<float4> shared_node_pack(shared.node_pack);

    unsigned int& nEdgesCounter = *payload.nEdgesCounter;

    const float deltaPhi = d_bin_pair_dphi[blockIndex];

    const unsigned int begin_bin1 = d_bin_pair_views[4u * blockIndex];
    const unsigned int begin_bin2 = d_bin_pair_views[4u * blockIndex + 2u];
    const unsigned int num_nodes1 =
        d_bin_pair_views[4u * blockIndex + 1u] - begin_bin1;
    const unsigned int num_nodes2 =
        d_bin_pair_views[4u * blockIndex + 3u] - begin_bin2;

    for (unsigned int node1_idx = threadIndex; node1_idx < num_nodes1;
         node1_idx += blockSize) {
        // loading a chunk of nodes1 into shared mem buffers
        const unsigned int gidx = node1_idx + begin_bin1;
        shared_node_pack[node1_idx] = d_node_params[gidx];
        shared_phi[node1_idx] = d_node_phi[gidx];
    }

    barrier.blockBarrier();

    const float phi0 = shared_phi[0];
    const float phiN = shared_phi[num_nodes1 - 1u];
    const float phi_bin_width =
        traccc::device::TWO_PI_F / static_cast<float>(payload.nPhiBins);
    const float break_threshold =
        deltaPhi + phi_bin_width - traccc::device::PI_F;

    unsigned int last_n1 = 0u;  // initial value for the sliding window

    for (unsigned int n2Idx = threadIndex; n2Idx < num_nodes2;
         n2Idx += blockSize) {

        const unsigned int globalIdx2 = begin_bin2 + n2Idx;
        const float phi2 = d_node_phi[globalIdx2];

        unsigned int n1Idx = last_n1;

        float min_phi1 = phi2 - deltaPhi;
        float max_phi1 = phi2 + deltaPhi;

        if (min_phi1 < -traccc::device::PI_F) {
            min_phi1 += traccc::device::TWO_PI_F;
        }
        if (max_phi1 > traccc::device::PI_F) {
            max_phi1 -= traccc::device::TWO_PI_F;
        }
        const bool boundary = max_phi1 < min_phi1;  // +/- pi wraparound

        // expand over nearest bin boundary
        max_phi1 += phi_bin_width;
        min_phi1 -= phi_bin_width;

        if (!boundary) {
            if (phi0 > max_phi1) {
                continue;
            }
            if (phiN < min_phi1) {
                // if bin1 can't be part of a wraparound
                // from a high-phi node skip it
                if (phi0 > break_threshold) {
                    break;
                }
                continue;
            }
        } else {
            if (phi0 < max_phi1) {
                // if not too large for lower wraparound don't skip it
                n1Idx = 0u;
            } else if (phiN < min_phi1) {
                continue;
            }
        }

        const float4 np2 = d_node_params[globalIdx2];
        if (!boundary) {
            for (; n1Idx < num_nodes1; n1Idx++) {
                const float phi1 = shared_phi[n1Idx];

                if (phi1 > max_phi1) {
                    break;
                }
                if (phi1 < min_phi1) {
                    continue;
                }
                last_n1 = n1Idx;

                const float4 np1 = shared_node_pack[n1Idx];
                detail::gbts_check_edge_candidate(
                    np1, np2, d_edge_nodes, d_edge_params, d_num_outgoing_edges,
                    nEdgesCounter, globalIdx2, begin_bin1, n1Idx, phi1, phi2,
                    deltaPhi, payload.gbts_make_graph_edges_params,
                    payload.nMaxEdges);
            }
        } else {
            for (; n1Idx < num_nodes1; n1Idx++) {
                const float phi1 = shared_phi[n1Idx];
                if (phi1 > max_phi1 && phi1 < min_phi1) {
                    continue;
                }
                last_n1 = n1Idx;

                const float4 np1 = shared_node_pack[n1Idx];
                detail::gbts_check_edge_candidate(
                    np1, np2, d_edge_nodes, d_edge_params, d_num_outgoing_edges,
                    nEdgesCounter, globalIdx2, begin_bin1, n1Idx, phi1, phi2,
                    deltaPhi, payload.gbts_make_graph_edges_params,
                    payload.nMaxEdges);
            }
        }
    }
}

}  // namespace traccc::device
