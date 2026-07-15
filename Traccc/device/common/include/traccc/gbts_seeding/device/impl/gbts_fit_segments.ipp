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
#include "traccc/device/concepts/thread_id.hpp"
#include "traccc/gbts_seeding/device/details/gbts_create_seed_candidate.hpp"
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"
#include "traccc/gbts_seeding/gbts_types.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>
#include <vecmem/memory/device_atomic_ref.hpp>

// System include(s).
#include <cstring>

namespace traccc::device {

namespace details {

// ===========================================================================
// edgeState -- Kalman-filter state for a track-segment fit.
// ===========================================================================
// Two decoupled Kalman fits per seed, in a frame along the first doublet
//   x-y bending plane: parabola fit, m_X = [eta, deta/dA, curvature (~q/pT)]
//   r-z plane:         line fit,     m_Y = [z, tau]  (tau = dz/dr = cot theta)
// m_Cx / m_Cy are the packed symmetric covariances of m_X / m_Y.
//
struct edgeState {

    TRACCC_HOST_DEVICE inline float& m_Cx(const int i, const int j) {
        return Cx[i + j + 1 * (i != 0) * (j != 0)];
    }
    TRACCC_HOST_DEVICE inline float& m_Cy(const int i, const int j) {
        return Cy[i + j];
    }
    TRACCC_HOST_DEVICE inline const float& m_Cx(const int i,
                                                const int j) const {
        return Cx[i + j + 1 * (i != 0) * (j != 0)];
    }
    TRACCC_HOST_DEVICE inline const float& m_Cy(const int i,
                                                const int j) const {
        return Cy[i + j];
    }

    // Initialize the edgeState from the two seed nodes.
    // The nodes are float4 (x, y, z, width)
    TRACCC_HOST_DEVICE inline void initialize(
        const traccc::float4& node1_params,
        const traccc::float4& node2_params) {
        m_J = 0.0f;
        m_head_node_type = (node1_params.w < 0);

        // differenc in x-y and its length L.
        const float dx = node1_params.x - node2_params.x;
        const float dy = node1_params.y - node2_params.y;
        const float L = math::sqrt(dx * dx + dy * dy);

        // Radius of the two seed nodes.
        const float r1 = math::sqrt(node1_params.x * node1_params.x +
                                    node1_params.y * node1_params.y);
        const float r2 = math::sqrt(node2_params.x * node2_params.x +
                                    node2_params.y * node2_params.y);

        m_s = dy / L;
        m_c = dx / L;

        m_refX = node2_params.x * m_c + node2_params.y * m_s;
        m_refY = r2;

        // Bending parabola, m_X = [eta, deta/dA, curvature]:
        //   m_X[0] = eta = -node2.x*m_s + node2.y*m_c
        //   m_X[1] = deta/dA = bending slope     -> 0
        //   m_X[2] = curvature (~q/pT)           -> 0
        m_X[0] = -node2_params.x * m_s + node2_params.y * m_c;
        m_X[1] = 0.0f;
        m_X[2] = 0.0f;

        // r-z line, m_Y = [z, tau]:
        //   m_Y[0] = z
        //   m_Y[1] = tau
        m_Y[0] = node2_params.z;
        m_Y[1] = (node1_params.z - node2_params.z) / (r1 - r2);

        std::memset(&m_Cx(0, 0), 0, sizeof(Cx));
        std::memset(&m_Cy(0, 0), 0, sizeof(Cy));

        m_Cx(0, 0) = 0.25f;
        m_Cx(1, 1) = 0.001f;
        m_Cx(2, 2) = 0.001f;

        m_Cy(0, 0) = 1.5f;
        m_Cy(1, 1) = 0.001f;
    }

    float m_X[3], m_Y[2];
    float m_c, m_s, m_refX, m_refY;
    float m_J;
    bool m_head_node_type;
    float Cx[6];
    float Cy[3];
};

TRACCC_HOST_DEVICE inline bool gbts_kalman_update(
    edgeState* new_ts, const edgeState* ts, const traccc::float4& node1_params,
    const gbts_fit_segments_params& KF_params, const float max_z0) {

    const float tau2 = ts->m_Y[1] * ts->m_Y[1];
    const float invSin2 = 1 + tau2;

    const float lenCorr = (node1_params.w != -1) ? invSin2 : invSin2 / tau2;
    const float minPtFrac = fabsf(ts->m_X[2]) * KF_params.inv_max_curvature;

    const float corrMS = KF_params.sigmaMS * minPtFrac;  // MS angle at this pT
    const float sigma2 =
        KF_params.radLen * lenCorr * corrMS * corrMS;  // MS var

    // Inflate the slope variance of each filter by the MS variance.
    const float m_Cx11 = ts->m_Cx(1, 1) + sigma2;
    const float m_Cy11 = ts->m_Cy(1, 1) + sigma2;

    float mx, my;
    const float r = math::sqrt(node1_params.x * node1_params.x +
                               node1_params.y * node1_params.y);

    new_ts->m_refX = node1_params.x * ts->m_c + node1_params.y * ts->m_s;
    mx = -node1_params.x * ts->m_s + node1_params.y * ts->m_c;
    new_ts->m_refY = r;
    my = node1_params.z;

    // Step variables from the previous state: A along the doublet
    const float A = new_ts->m_refX - ts->m_refX;
    const float B = (0.5f) * A * A;
    const float dr = new_ts->m_refY - ts->m_refY;

    // Extrapolate the parabola to the hit:
    //   eta += slope*A + curvature*A^2/2,  slope += curvature*A.
    new_ts->m_X[0] = ts->m_X[0] + ts->m_X[1] * A + ts->m_X[2] * B;
    new_ts->m_X[1] = ts->m_X[1] + ts->m_X[2] * A;
    new_ts->m_X[2] = ts->m_X[2];

    // Propagate its covariance, Cx' = F Cx F^T (MS folded into m_Cx11).
    //   F = [ 1  A  B ]   (B = A^2/2)
    //       [ 0  1  A ]
    //       [ 0  0  1 ]
    new_ts->m_Cx(0, 0) = ts->m_Cx(0, 0) + 2 * ts->m_Cx(0, 1) * A +
                         2 * ts->m_Cx(0, 2) * B + A * m_Cx11 * A +
                         2 * A * ts->m_Cx(1, 2) * B + B * ts->m_Cx(2, 2) * B;
    new_ts->m_Cx(0, 1) = ts->m_Cx(0, 1) + m_Cx11 * A + ts->m_Cx(1, 2) * B +
                         ts->m_Cx(0, 2) * A + A * ts->m_Cx(1, 2) * A +
                         A * ts->m_Cx(2, 2) * B;
    new_ts->m_Cx(0, 2) =
        ts->m_Cx(0, 2) + ts->m_Cx(1, 2) * A + ts->m_Cx(2, 2) * B;
    new_ts->m_Cx(1, 1) =
        m_Cx11 + 2 * A * ts->m_Cx(1, 2) + A * ts->m_Cx(2, 2) * A;
    new_ts->m_Cx(1, 2) = ts->m_Cx(1, 2) + ts->m_Cx(2, 2) * A;
    new_ts->m_Cx(2, 2) = ts->m_Cx(2, 2);

    // Extrapolate the r-z line to the hit:  z += tau*dr  (tau unchanged).
    new_ts->m_Y[0] = ts->m_Y[0] + ts->m_Y[1] * dr;
    new_ts->m_Y[1] = ts->m_Y[1];
    // Propagate its covariance (MS folded into m_Cy11).
    new_ts->m_Cy(0, 0) =
        ts->m_Cy(0, 0) + 2 * ts->m_Cy(0, 1) * dr + dr * m_Cy11 * dr;
    new_ts->m_Cy(0, 1) = ts->m_Cy(0, 1) + dr * m_Cy11;
    new_ts->m_Cy(1, 1) = m_Cy11;

    // Residuals = measured - predicted (eta in the bending plane, z in r-z).
    const float resid_x = mx - new_ts->m_X[0];
    const float resid_y = my - new_ts->m_Y[0];

    // r-z measurement error: sigma_y (*tau for strips, w < 0).
    float sigma_rz = 0;
    if (!ts->m_head_node_type) {
        sigma_rz = KF_params.sigma_y;
    } else {
        sigma_rz = KF_params.sigma_y * ts->m_Y[1];
    }

    const float inv_Dx =
        new_ts->m_Cx(0, 0) + KF_params.sigma_x * KF_params.sigma_x;
    const float Dx = 1 / inv_Dx;
    const float Dy = 1 / (new_ts->m_Cy(0, 0) + sigma_rz * sigma_rz);

    // Per-hit chi2 = resid^2 * D; reject the hit if either plane is too large.
    const float dchi2_x = resid_x * resid_x * Dx;
    const float dchi2_y = resid_y * resid_y * Dy;

    if (dchi2_x > KF_params.maxDChi2_x || dchi2_y > KF_params.maxDChi2_y) {
        return false;
    }

    // Accumulate seed quality: +add_hit per accepted hit, -weighted chi2.
    new_ts->m_J = ts->m_J + (KF_params.add_hit - dchi2_x * KF_params.weight_x -
                             dchi2_y * KF_params.weight_y);

    for (unsigned int i = 0u; i < 3u; i++) {
        new_ts->m_X[i] += Dx * new_ts->m_Cx(0, static_cast<int>(i)) * resid_x;
    }

    // pT cut: reject once the updated curvature exceeds curv_max.
    if (math::fabs(new_ts->m_X[2]) * KF_params.inv_max_curvature > 1.0f) {
        return false;
    }

    // Measurement update of the r-z line.
    for (unsigned int i = 0u; i < 2u; i++) {
        new_ts->m_Y[i] += Dx * new_ts->m_Cy(0, static_cast<int>(i)) * resid_y;
    }

    // z0 cut: extrapolate the line back to r = 0 and reject large |z0|.
    const float z0 = new_ts->m_Y[0] - new_ts->m_refY * ts->m_Y[1];
    if (math::fabs(z0) > max_z0) {
        return false;
    }

    // Covariance update Cx <- (I - K H) Cx, in place from high index to low.
    //   I = identity, H = [1 0 0],
    //   K = Dx * [Cx(0,0), Cx(0,1), Cx(0,2)]^T  (Kalman gain).
    new_ts->m_Cx(2, 2) = Dx * (new_ts->m_Cx(2, 2) * inv_Dx -
                               new_ts->m_Cx(0, 2) * new_ts->m_Cx(0, 2));
    new_ts->m_Cx(1, 2) = Dx * (new_ts->m_Cx(1, 2) * inv_Dx -
                               new_ts->m_Cx(0, 1) * new_ts->m_Cx(0, 2));
    new_ts->m_Cx(1, 1) = Dx * (new_ts->m_Cx(1, 1) * inv_Dx -
                               new_ts->m_Cx(0, 1) * new_ts->m_Cx(0, 1));
    new_ts->m_Cx(0, 2) = Dx * (new_ts->m_Cx(0, 2) * inv_Dx -
                               new_ts->m_Cx(0, 0) * new_ts->m_Cx(0, 2));
    new_ts->m_Cx(0, 1) = Dx * (new_ts->m_Cx(0, 1) * inv_Dx -
                               new_ts->m_Cx(0, 0) * new_ts->m_Cx(0, 1));
    new_ts->m_Cx(0, 0) *= Dx * (KF_params.sigma_x * KF_params.sigma_x);

    // Covariance update Cy <- (I - K H) Cy (uses Dy).
    new_ts->m_Cy(1, 1) -= Dy * new_ts->m_Cy(0, 1) * new_ts->m_Cy(0, 1);
    new_ts->m_Cy(0, 1) -= Dy * new_ts->m_Cy(0, 0) * new_ts->m_Cy(0, 1);
    new_ts->m_Cy(0, 0) -= Dy * new_ts->m_Cy(0, 0) * new_ts->m_Cy(0, 0);

    // Carry the frozen frame forward
    new_ts->m_c = ts->m_c;
    new_ts->m_s = ts->m_s;
    new_ts->m_head_node_type = (node1_params.w < 0);

    return true;
}

}  // namespace details

template <concepts::thread_id1 thread_id_t>
TRACCC_HOST_DEVICE inline void gbts_fit_segments(
    const thread_id_t& thread_id, const gbts_fit_segments_payload& payload) {

    const vecmem::device_vector<const float4> d_sp_reduced(payload.reducedSP);
    const vecmem::device_vector<const unsigned int> d_output_graph(
        payload.output_graph);
    const vecmem::device_vector<const int2> d_path_store(payload.path_store);

    const gbts_fit_segments_params& fit_params =
        payload.gbts_fit_segments_params;
    const float max_z0 = payload.max_z0;

    // Row-major output graph: each edge owns a contiguous block of
    // edge_size = 2 + 1 + max_num_neighbours ints.
    const unsigned int edge_size = 2u + 1u + payload.max_num_neighbours;

    const unsigned int store_size = *payload.nPathStoreSize;

    const unsigned int globalIdx = thread_id.getGlobalThreadIdX();
    const unsigned int blockDimX = thread_id.getBlockDimX();
    const unsigned int gridDimX = thread_id.getGridDimX();

    for (unsigned int globalIndex = globalIdx;
         globalIndex + payload.nTerminusEdges < store_size;
         globalIndex += blockDimX * gridDimX) {

        const unsigned int path_idx = globalIndex + payload.nTerminusEdges;

        unsigned char length = 1;
        bool toggle = false;
        details::edgeState state1;
        details::edgeState state2;

        int2 path = d_path_store[path_idx];

        const unsigned int nodeidx1 =
            d_output_graph[edge_size * static_cast<unsigned int>(path.x) +
                           gbts_consts::node1];
        traccc::float4 node1 = d_sp_reduced[nodeidx1];
        const unsigned int nodeidx2 =
            d_output_graph[edge_size * static_cast<unsigned int>(path.x) +
                           gbts_consts::node2];
        traccc::float4 node2 = d_sp_reduced[nodeidx2];

        state1.initialize(node2, node1);
        while (path.y >= 0) {
            path = d_path_store[static_cast<unsigned int>(path.y)];
            node2 = d_sp_reduced
                [d_output_graph[edge_size * static_cast<unsigned int>(path.x) +
                                gbts_consts::node2]];
            if (toggle) {
                if (!details::gbts_kalman_update(&state1, &state2, node2,
                                                 fit_params, max_z0)) {
                    state1 = state2;
                    break;
                }
            } else if (!details::gbts_kalman_update(&state2, &state1, node2,
                                                    fit_params, max_z0)) {
                break;
            }
            toggle = !toggle;
            length++;
        }
        if (length < payload.minLevel) {
            continue;
        }
        int qual = 0;
        if (toggle) {
            qual = static_cast<int>(fit_params.qual_scale * state2.m_J);
        } else {
            qual = static_cast<int>(fit_params.qual_scale * state1.m_J);
        }
        const unsigned int prop_idx =
            vecmem::device_atomic_ref<unsigned int>(*payload.nPropsCounter)
                .fetch_add(1u);
        details::gbts_create_seed_candidate(
            qual, static_cast<int>(path_idx), prop_idx, payload.seed_ambiguity,
            payload.seed_proposals, payload.edge_bids, payload.path_store, 1);
    }
}

}  // namespace traccc::device
