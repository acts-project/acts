/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// cuda includes
#include <cuda.h>
#include <cuda_fp16.h>
#include <cuda_runtime.h>
#include <math_constants.h>
#include <vector_functions.h>

namespace traccc::cuda::kernels {

struct edgeState {

    __device__ inline void initialize(const float4& node1_params,
                                      const float4& node2_params);

    __device__ inline float& m_Cx(const int i, const int j) {
        return Cx[i + j + 1 * (i != 0) * (j != 0)];
    }
    __device__ inline float& m_Cy(const int i, const int j) {
        return Cy[i + j];
    }
    __device__ inline const float& m_Cx(const int i, const int j) const {
        return Cx[i + j + 1 * (i != 0) * (j != 0)];
    }
    __device__ inline const float& m_Cy(const int i, const int j) const {
        return Cy[i + j];
    }

    float m_X[3], m_Y[2];
    float m_c, m_s, m_refX, m_refY;

    float m_J;

    bool m_head_node_type;

    // upper triangle of the Cov matrix for the parabola in the x,y plane since
    // symetry gives the rest
    float Cx[6];  //(0,0), (0,1), (0,2), (1,1), (1,2), (2,2)
    // Cov matrix for the linear fit of eta and z
    float Cy[3];  //(0,0), (0,1), (1,1)
};

struct Tracklet {
    unsigned int nodes[traccc::device::gbts_consts::max_cca_iter + 1];
    int size;
};

/** @brief Performs one iteration of the CCA over the graph to calculate
 * potential seed length
 *
 * also counts the size of d_path_store to describe all the paths back
 * to the inner-most (terminus) edges
 *
 *  @param[in] d_output_graph see comments in device_context.h
 *  @param[in] d_levels is the maximum seed length originating at each edge
 *  @param[in/out] d_active_edges stores flags for edges that need more CCA
 * iterations
 *  @param[out] d_outgoing_paths [#paths needed, is-terminus]
 */
__global__ static void CCA_IterationKernel(const int* d_output_graph,
                                           char* d_levels, char* d_active_edges,
                                           short2* d_outgoing_paths, int iter,
                                           unsigned int nEdges,
                                           unsigned int max_num_neighbours,
                                           int minLevel) {

    unsigned int edge_size = 2 + 1 + max_num_neighbours;

    int toggle = iter % 2;
    int levelLoad = toggle * nEdges;
    int levelStore = (1 - toggle) * nEdges;

    for (int edgeIdx = threadIdx.x + blockIdx.x * blockDim.x; edgeIdx < nEdges;
         edgeIdx += blockDim.x * gridDim.x) {

        if (iter != 0) {
            if (d_active_edges[edgeIdx] != iter) {
                continue;
            }
        }
        int edge_pos = edge_size * edgeIdx;

        int nNei = d_output_graph[edge_pos + traccc::device::gbts_consts::nNei];

        char next_level = d_levels[levelLoad + edgeIdx];

        bool localChange = false;
        for (int nIdx = 0; nIdx < nNei;
             nIdx++) {  // loop over neighbouring edges

            int nextEdgeIdx =
                d_output_graph[edge_pos +
                               traccc::device::gbts_consts::nei_start + nIdx];
            char forward_level = d_levels[levelLoad + nextEdgeIdx];

            if (next_level == forward_level) {
                next_level = forward_level + 1;
                localChange = true;
                break;
            }
        }
        // add all remianing edges to level_views on the last iteration
        if (localChange) {
            if (iter == traccc::device::gbts_consts::max_cca_iter - 1) {
                // shorten paths longer than max_cca_iter
                d_outgoing_paths[edgeIdx].y = -1;
                d_active_edges[edgeIdx] = -1;
            } else {
                // flag edge for the next iteration
                d_active_edges[edgeIdx] = iter + 1;
            }
        } else {
            d_active_edges[edgeIdx] = -1;
            short out_paths = 0;
            for (int nIdx = 0; nIdx < nNei; ++nIdx) {
                int nextEdgeIdx =
                    d_output_graph[edge_pos +
                                   traccc::device::gbts_consts::nei_start +
                                   nIdx];
                if (next_level == 1 + d_levels[nextEdgeIdx]) {
                    // calculate the #d_state_store nodes for segment extraction
                    // starting at this edge
                    out_paths += 1 + d_outgoing_paths[nextEdgeIdx].x;
                }
                // flag as not terminus edge
                d_outgoing_paths[nextEdgeIdx].y = -1;
            }
            // flag as long enough segement to become a seed
            d_outgoing_paths[edgeIdx] =
                make_short2(out_paths, (next_level >= minLevel) - 1);
        }
        // store new level and ensure all final
        // levels are on both sides of the array
        d_levels[levelStore + edgeIdx] = next_level;
    }
}

void __global__ count_terminus_edges(int2* d_path_store,
                                     short2* d_outgoing_paths,
                                     unsigned int* d_counters,
                                     unsigned int nEdges) {

    // count in shared first to reduce global atomics
    __shared__ int outgoingCount;

    if (threadIdx.x == 0) {
        outgoingCount = 0;
    }
    __syncthreads();

    int edge_idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (edge_idx < nEdges) {

        short2 out_paths = d_outgoing_paths[edge_idx];
        // only count terminus edges that could lead to a seed
        // fill the first part of path_store so fitting can skip go-nowhere
        // paths
        if (out_paths.y != -1) {
            d_outgoing_paths[edge_idx].y = atomicAdd(&d_counters[7], 1);
            atomicAdd(&outgoingCount, out_paths.x);
        }
    }
    __syncthreads();
    if (threadIdx.x == 0) {
        atomicAdd(&d_counters[6], outgoingCount);
    }
}

void __global__ add_terminus_to_path_store(int2* d_path_store,
                                           short2* d_outgoing_paths,
                                           unsigned int* d_counters,
                                           unsigned int nEdges) {

    int edge_idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (edge_idx >= nEdges) {
        return;
    }
    short2 out_paths = d_outgoing_paths[edge_idx];
    if (out_paths.y == -1) {
        return;
    }
    // -1 flags as the terminus of a path
    d_path_store[out_paths.y] = make_int2(edge_idx, -1);
}

// each node in the path_store defines a unique path through the graph
// but includes subsets like 0->1->2 and 1->2
// The paths order the edges in reverse order compared to graph linking
// extracting all paths here allows for fitting to occor down know paths in
// registers
void __global__ fill_path_store(int2* d_path_store, int* d_output_graph,
                                char* d_levels, unsigned int* d_counters,
                                unsigned int nTerminus,
                                unsigned int nTerminusPerBlock,
                                unsigned int max_num_neighbours,
                                unsigned int nPaths) {

    __shared__ int2 live_paths[traccc::device::gbts_consts::live_path_buffer];
    __shared__ int n_live_paths;

    if (threadIdx.x == 0) {
        n_live_paths = 0;
    }
    __syncthreads();

    int edge_size = 2 + 1 + max_num_neighbours;
    unsigned int path_idx = threadIdx.x + blockIdx.x * nTerminusPerBlock;
    // populate live_paths with terminus to start exploration from
    if (threadIdx.x < nTerminusPerBlock && path_idx < nTerminus) {
        int2 path = d_path_store[path_idx];
        int nNei = d_output_graph[traccc::device::gbts_consts::nNei +
                                  edge_size * path.x];
        char level = d_levels[path.x];
        for (int nei = 0; nei < nNei; ++nei) {
            int edge_idx =
                d_output_graph[traccc::device::gbts_consts::nei_start + nei +
                               edge_size * path.x];
            // only search down longest path
            if (level != d_levels[edge_idx] + 1) {
                continue;
            }
            int live_idx = atomicAdd(&n_live_paths, 1);
            if (live_idx >= traccc::device::gbts_consts::live_path_buffer) {
                break;
            }
            int new_path_idx = atomicAdd(&d_counters[7], 1);
            // head edge idx, link back
            d_path_store[new_path_idx] = make_int2(edge_idx, path_idx);
            live_paths[live_idx] = make_int2(edge_idx, new_path_idx);
        }
    }
    __syncthreads();

    int2 path = make_int2(0, 0);
    bool has_path = false;

    while (n_live_paths > 0) {
        has_path = false;
        if (threadIdx.x == 0) {
            n_live_paths = min(n_live_paths,
                               traccc::device::gbts_consts::live_path_buffer);
        }
        __syncthreads();
        // get path
        if (threadIdx.x < n_live_paths) {
            path = live_paths[n_live_paths - threadIdx.x - 1];
            has_path = true;
        }
        __syncthreads();
        if (threadIdx.x == 0) {
            n_live_paths =
                n_live_paths < blockDim.x ? 0 : n_live_paths - blockDim.x;
        }
        __syncthreads();
        if (has_path) {
            int nNei = d_output_graph[traccc::device::gbts_consts::nNei +
                                      edge_size * path.x];
            char level = d_levels[path.x];
            for (int nei = 0; nei < nNei; ++nei) {
                int edge_idx =
                    d_output_graph[traccc::device::gbts_consts::nei_start +
                                   nei + edge_size * path.x];
                // only search down longest segments
                if (level != d_levels[edge_idx] + 1) {
                    continue;
                }
                path_idx = atomicAdd(&d_counters[7], 1);
                if (path_idx >= nPaths) {
                    break;
                }
                int live_idx = atomicAdd(&n_live_paths, 1);
                if (live_idx >= traccc::device::gbts_consts::live_path_buffer) {
                    break;
                }
                // head edge idx, link back
                d_path_store[path_idx] = make_int2(edge_idx, path.y);
                live_paths[live_idx] = make_int2(edge_idx, path_idx);
            }
        }
        // wait for live_paths to repopulate
        __syncthreads();
    }
}

/** @brief initialize the Kalman filter for this new edgeState from the starting
 * edge (2 nodes)
 *
 *  @param[in] node1_params / node2_params the nodes of the starting edge. Node
 *  1 is the inner node and we filter outside in
 */
__device__ inline void edgeState::initialize(
    const float4& node1_params, const float4& node2_params) {  // x, y, z,type

    m_J = 0.0f;
    m_head_node_type = (node1_params.w < 0);
    // n2->n1

    float dx = node1_params.x - node2_params.x;
    float dy = node1_params.y - node2_params.y;
    float L = sqrtf(dx * dx + dy * dy);

    float r1 = sqrtf(node1_params.x * node1_params.x +
                     node1_params.y * node1_params.y);
    float r2 = sqrtf(node2_params.x * node2_params.x +
                     node2_params.y * node2_params.y);

    m_s = dy / L;
    m_c = dx / L;

    // transform for extrapolation and update
    //  x' =  x*m_c + y*m_s
    //  y' = -x*m_s + y*m_c

    m_refY = r2;
    m_refX = node2_params.x * m_c + node2_params.y * m_s;

    // X-state: y, dy/dx, d2y/dx2

    m_X[0] = -node2_params.x * m_s + node2_params.y * m_c;
    m_X[1] = 0.0f;
    m_X[2] = 0.0f;

    // Y-state: z, dz/dr

    m_Y[0] = node2_params.z;
    m_Y[1] = (node1_params.z - node2_params.z) / (r1 - r2);

    memset(&m_Cx(0, 0), 0, sizeof(Cx));
    memset(&m_Cy(0, 0), 0, sizeof(Cy));

    m_Cx(0, 0) = 0.25f;
    m_Cx(1, 1) = 0.001f;
    m_Cx(2, 2) = 0.001f;

    m_Cy(0, 0) = 1.5f;
    m_Cy(1, 1) = 0.001f;
}

/** Attempts to update the edgeState to include node1
 *
 *  This is a Kalamn filter update fitting to strait line in z,r and a parabola
 * in x',y' with the transformation defined by the first edge when the state is
 * initialized It's main output is the m_J seed quality used for disambiguation
 *  Seed extraction goes outside in
 *
 *  @param[out] new_ts output edgeState for the updated seed including node1
 *  this is also used for register space
 *  @param[in] ts input edgeState
 *  @param[in] node1_params params of the inner node of the new edge to be added
 *  to the seed
 */
inline __device__ bool update(edgeState* new_ts, const edgeState* ts,
                              const float4& node1_params,
                              // node params are x, y, z, type
                              const gbts_seed_extraction_params& KF_params) {

    float tau2 = ts->m_Y[1] * ts->m_Y[1];
    float invSin2 = 1 + tau2;

    float lenCorr = (node1_params.w != -1) ? invSin2 : invSin2 / tau2;
    float minPtFrac = fabsf(ts->m_X[2]) * KF_params.inv_max_curvature;

    float corrMS = KF_params.sigmaMS * minPtFrac;
    float sigma2 = KF_params.radLen * lenCorr * corrMS * corrMS;  // /invSin2;

    // add ms.
    float m_Cx11 = ts->m_Cx(1, 1) + sigma2;
    float m_Cy11 = ts->m_Cy(1, 1) + sigma2;

    // extrapolation

    // float refX, refY;
    float mx, my;

    float r = sqrtf(node1_params.x * node1_params.x +
                    node1_params.y * node1_params.y);

    // using new_ts as register storage where possible
    new_ts->m_refX = node1_params.x * ts->m_c + node1_params.y * ts->m_s;
    mx = -node1_params.x * ts->m_s + node1_params.y * ts->m_c;  // measured X[0]
    new_ts->m_refY = r;
    my = node1_params.z;  // measured Y[0]

    float A = new_ts->m_refX - ts->m_refX;
    float B = (0.5f) * A * A;
    float dr = new_ts->m_refY - ts->m_refY;

    new_ts->m_X[0] = ts->m_X[0] + ts->m_X[1] * A + ts->m_X[2] * B;
    new_ts->m_X[1] = ts->m_X[1] + ts->m_X[2] * A;
    new_ts->m_X[2] = ts->m_X[2];

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

    new_ts->m_Y[0] = ts->m_Y[0] + ts->m_Y[1] * dr;
    new_ts->m_Y[1] = ts->m_Y[1];

    new_ts->m_Cy(0, 0) =
        ts->m_Cy(0, 0) + 2 * ts->m_Cy(0, 1) * dr + dr * m_Cy11 * dr;

    new_ts->m_Cy(0, 1) = ts->m_Cy(0, 1) + dr * m_Cy11;
    new_ts->m_Cy(1, 1) = m_Cy11;

    // chi2 test
    float resid_x = mx - new_ts->m_X[0];
    float resid_y = my - new_ts->m_Y[0];

    float sigma_rz = 0;

    if (!ts->m_head_node_type) {
        // barrel TO-DO: split into barrel Pixel and barrel SCT
        sigma_rz = KF_params.sigma_y;
    } else {
        sigma_rz = KF_params.sigma_y * ts->m_Y[1];
    }

    float inv_Dx = new_ts->m_Cx(0, 0) + KF_params.sigma_x * KF_params.sigma_x;
    float Dx = 1 / inv_Dx;

    float Dy = 1 / (new_ts->m_Cy(0, 0) + sigma_rz * sigma_rz);

    float dchi2_x = resid_x * resid_x * Dx;
    float dchi2_y = resid_y * resid_y * Dy;

    if (dchi2_x > KF_params.maxDChi2_x || dchi2_y > KF_params.maxDChi2_y) {
        return false;
    }

    // state update
    new_ts->m_J = ts->m_J + (KF_params.add_hit - dchi2_x * KF_params.weight_x -
                             dchi2_y * KF_params.weight_y);

    for (int i = 0; i < 3; i++) {
        new_ts->m_X[i] += Dx * new_ts->m_Cx(0, i) * resid_x;
    }

    if (fabsf(new_ts->m_X[2]) * KF_params.inv_max_curvature > 1.0f) {
        return false;
    }

    for (int i = 0; i < 2; i++) {
        new_ts->m_Y[i] += Dx * new_ts->m_Cy(0, i) * resid_y;
    }

    float z0 = new_ts->m_Y[0] - new_ts->m_refY * ts->m_Y[1];
    if (fabsf(z0) > KF_params.max_z0) {
        return false;
    }

    // less loss from float precsion this way (helps prevent sign change)
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

    new_ts->m_Cy(1, 1) -= Dy * new_ts->m_Cy(0, 1) * new_ts->m_Cy(0, 1);
    new_ts->m_Cy(0, 1) -= Dy * new_ts->m_Cy(0, 0) * new_ts->m_Cy(0, 1);
    new_ts->m_Cy(0, 0) -= Dy * new_ts->m_Cy(0, 0) * new_ts->m_Cy(0, 0);

    new_ts->m_c = ts->m_c;
    new_ts->m_s = ts->m_s;
    new_ts->m_head_node_type = (node1_params.w < 0);

    return true;
}

/** @brief Performs seed disambiguation through seeds biding to use edges with
 * seed quality
 *
 *  @param[in] qual is the quality metric output by the Kalman filter
 *  @param[in] path_idx the index of inital path for this seed
 *  @param[in] d_path_store stores the path each seed took through the graph in
 *  reverse order
 *  @param[in] prop_idx the index of this new seeds proposition in
 *  d_seed_proposals
 *  @param[in/out] d_edge_bids is [int_m_J, prop_idx] so that atomicMax will
 *  swap it out with higher quality bids. The index is then used to flag the
 *  replaced seed as maybe fake
 *  @param[out] d_seed_proposals stores the information needed to construct an
 *  output Tracklet for this seed
 *  @param[out] d_seed_ambiguity here is 0 if the seed is the highest quality
 *  seed using all of its edges and -1 otherwise
 */
inline __device__ void add_seed_proposal(const int qual, const int path_idx,
                                         const unsigned int prop_idx,
                                         char* d_seed_ambiguity,
                                         int2* d_seed_proposals,
                                         unsigned long long int* d_edge_bids,
                                         const int2* d_path_store,
                                         char depth = -1) {

    // new seed bids for its edges
    d_seed_proposals[prop_idx] = make_int2(qual, path_idx);
    d_seed_ambiguity[prop_idx] = 0;
    __threadfence();  // ensure above proposal info is written before biding

    unsigned long long int seed_bid =
        (static_cast<unsigned long long int>(qual) << 32) |
        (static_cast<unsigned long long int>(prop_idx));

    // dummy path to start the loop
    int2 path = make_int2(0, d_seed_proposals[prop_idx].y);
    while (path.y >= 0 && depth != 0) {
        path = d_path_store[path.y];
        depth--;

        unsigned long long int competing_offer =
            atomicMax(&d_edge_bids[path.x], seed_bid);
        if (competing_offer > seed_bid) {
            d_seed_ambiguity[prop_idx] = -1;
        } else if (competing_offer != 0) {
            d_seed_ambiguity[competing_offer & 0xFFFFFFFFLL] = -1;
        }  // default bids are 0 so no need to replace
    }
}

void __global__ fit_segments(
    float4* d_sp_params, int* d_output_graph, int2* d_path_store,
    int2* d_seed_proposals, unsigned long long int* d_edge_bids,
    char* d_seed_ambiguity, char* d_levels, unsigned int* d_counters,
    unsigned int nTerminusEdges, int minLevel, unsigned int max_num_neighbours,
    gbts_seed_extraction_params seed_extraction_params) {

    // take an extracted path and fit it to produce a quality score
    unsigned int path_idx =
        threadIdx.x + blockIdx.x * blockDim.x + nTerminusEdges;
    if (path_idx >= d_counters[7]) {
        return;
    }
    int edge_size = 2 + 1 + max_num_neighbours;

    char length = 1;

    bool toggle = false;
    edgeState state1;
    edgeState state2;

    int2 path = d_path_store[path_idx];

    int nodeidx =
        d_output_graph[traccc::device::gbts_consts::node1 + edge_size * path.x];
    float4 node1 = d_sp_params[nodeidx];
    nodeidx =
        d_output_graph[traccc::device::gbts_consts::node2 + edge_size * path.x];
    float4 node2 = d_sp_params[nodeidx];

    state1.initialize(node2, node1);
    while (path.y >= 0) {
        path = d_path_store[path.y];

        node2 = d_sp_params[d_output_graph[traccc::device::gbts_consts::node2 +
                                           edge_size * path.x]];
        if (toggle) {
            if (!update(&state1, &state2, node2, seed_extraction_params)) {
                state1 = state2;
                break;
            }
        } else if (!update(&state2, &state1, node2, seed_extraction_params)) {
            break;
        }
        toggle = !toggle;
        length++;
    }
    // only keep long enough seeds
    if (length < minLevel) {
        return;
    }
    int qual = 0;
    if (toggle) {
        qual = static_cast<int>(seed_extraction_params.qual_scale * state2.m_J);
    } else {
        qual = static_cast<int>(seed_extraction_params.qual_scale * state1.m_J);
    }
    int prop_idx = atomicAdd(&d_counters[8], 1);
    // perform first round of bidding for disambiguation
    // only on the outermost edge
    add_seed_proposal(qual, path_idx, prop_idx, d_seed_ambiguity,
                      d_seed_proposals, d_edge_bids, d_path_store, 1);
}

void __global__ reset_edge_bids(int2* d_path_store, int2* d_seed_proposals,
                                unsigned long long int* d_edge_bids,
                                char* d_seed_ambiguity,
                                unsigned int* d_counters, int round) {

    int nProps = d_counters[8];
    // first round find best seed starting at each edge
    for (int prop_idx = threadIdx.x + blockIdx.x * blockDim.x;
         prop_idx < nProps; prop_idx += blockDim.x * gridDim.x) {

        char ambi = d_seed_ambiguity[prop_idx];
        if (round == -1) {
            if (ambi == 0) {
                // rebid 'best seed from edge' in later rounds
                d_seed_ambiguity[prop_idx] = 1;
                continue;
            } else {
                d_seed_ambiguity[prop_idx] = -2;
                // count rejected props to calculate nSeeds
                atomicAdd(&d_counters[9], 1);
                continue;
            }
        } else if (ambi == -2 | ambi == 0) {
            // only reset maybes
            continue;
        }
        int2 prop = d_seed_proposals[prop_idx];

        bool isgood = true;

        // dummy path to start the loop
        int2 path = make_int2(0, prop.y);
        while (path.y >= 0) {
            path = d_path_store[path.y];

            unsigned long long int best_bid = d_edge_bids[path.x];
            if (d_seed_ambiguity[best_bid & 0xFFFFFFFFLL] == 0) {
                isgood = false;
                break;
            }
        }
        if (isgood) {
            d_seed_ambiguity[prop_idx] = 1;
        }  // flag as maybe seed, shares with a loser
        else {
            d_seed_ambiguity[prop_idx] = -2;
            atomicAdd(&d_counters[9], 1);
            // definate fake, shares with a winner
        }
    }
}

// TO-DO?: reset prop count each iter and make new props like CCA active_edges
void __global__ seeds_rebid_for_edges(int2* d_path_store,
                                      int2* d_seed_proposals,
                                      unsigned long long int* d_edge_bids,
                                      char* d_seed_ambiguity,
                                      unsigned int nProps) {

    for (int prop_idx = threadIdx.x + blockIdx.x * blockDim.x;
         prop_idx < nProps; prop_idx += blockDim.x * gridDim.x) {

        char ambi = d_seed_ambiguity[prop_idx];
        if (ambi == -2 | ambi == 0) {
            // only rebid for maybes
            continue;
        }
        int2 prop = d_seed_proposals[prop_idx];

        add_seed_proposal(prop.x, prop.y, prop_idx, d_seed_ambiguity,
                          d_seed_proposals, d_edge_bids, d_path_store, -1);
    }
}

void __global__ seeds_bid_for_hits(int* d_output_graph, int2* d_seed_proposals,
                                   int2* d_path_store, char* d_seed_ambiguity,
                                   unsigned long long int* d_hit_bids,
                                   const unsigned int nProps, int edge_size) {

    for (unsigned int prop_idx = threadIdx.x + blockDim.x * blockIdx.x;
         prop_idx < nProps; prop_idx += gridDim.x * blockDim.x) {
        if (d_seed_ambiguity[prop_idx] == -2) {
            continue;
        }
        int2 prop = d_seed_proposals[prop_idx];
        unsigned long long int seed_bid =
            (static_cast<unsigned long long int>(prop.x) << 32) |
            (static_cast<unsigned long long int>(prop_idx));

        int2 path = make_int2(0, prop.y);
        while (path.y >= 0) {
            path = d_path_store[path.y];
            int sp_idx = d_output_graph[traccc::device::gbts_consts::node1 +
                                        edge_size * path.x];
            atomicMax(&d_hit_bids[sp_idx], seed_bid);
        }
        int sp_idx = d_output_graph[traccc::device::gbts_consts::node2 +
                                    edge_size * path.x];
        atomicMax(&d_hit_bids[sp_idx], seed_bid);
    }
}

inline __device__ float2 estimate_params(float4 sps[3]) {

    // conformal mapping with the center at the middle spacepoint

    float u[2], v[2];

    const float x0 = sps[1].x;
    const float y0 = sps[1].y;

    const float r0 = sqrtf(x0 * x0 + y0 * y0);

    const float cosA = x0 / r0;

    const float sinA = y0 / r0;

    for (unsigned int k = 0; k < 2; k++) {

        int sp_idx = (k == 1) ? 2 : k;

        const float dx = sps[sp_idx].x - x0;

        const float dy = sps[sp_idx].y - y0;

        const float r2_inv = 1.0f / (dx * dx + dy * dy);

        const float xn = dx * cosA + dy * sinA;

        const float yn = -dx * sinA + dy * cosA;

        u[k] = xn * r2_inv;
        v[k] = yn * r2_inv;
    }

    const float du = u[0] - u[1];
    if (du == 0.0f) {
        return make_float2(0.0f, 0.0f);
    }

    const float A = (v[0] - v[1]) / du;

    const float B = v[1] - A * u[1];

    // signed curvature in 1/m
    const float curv = 1000.0f * B / sqrtf(1 + A * A);
    const float cot_t = (sps[2].z - sps[1].z) /
                        (sqrtf(sps[2].x * sps[2].x + sps[2].y * sps[2].y) - r0);
    return make_float2(curv, cot_t);
}

void __global__ gbts_seed_conversion_kernel(
    int2* d_seed_proposals, char* d_seed_ambiguity, int2* d_path_store,
    int* d_output_graph, float4* d_sp_params,
    edm::seed_collection::view output_seeds, unsigned long long int* d_hit_bids,
    const unsigned int nProps, const unsigned int max_num_neighbours,
    const float dcurv_cut_m, const float force_dropout_max_curv_m,
    const float best_hit_frac, const float tight_bid_cot_threshold,
    const bool use_dropout) {

    int edge_size = 2 + 1 + max_num_neighbours;
    edm::seed_collection::device seeds_device(output_seeds);
    for (int prop_idx = threadIdx.x + blockIdx.x * blockDim.x;
         prop_idx < nProps; prop_idx += blockDim.x * gridDim.x) {

        if (d_seed_ambiguity[prop_idx] == -2) {
            // drop seeds that lost the bidding
            continue;
        }
        // collect seed hits and reject those that lost the hit bidding
        char best_for_hit = 0;
        Tracklet seed;
        seed.size = 0;
        // dummy path to start the loop
        int2 prop = d_seed_proposals[prop_idx];
        int2 path = make_int2(0, prop.y);
        while (path.y >= 0) {
            path = d_path_store[path.y];
            seed.nodes[seed.size++] =
                d_output_graph[traccc::device::gbts_consts::node1 +
                               edge_size * path.x];
            best_for_hit +=
                (prop_idx ==
                 (d_hit_bids[seed.nodes[seed.size - 1]] & 0xFFFFFFFFLL));
        }
        seed.nodes[seed.size++] =
            d_output_graph[traccc::device::gbts_consts::node2 +
                           edge_size * path.x];
        best_for_hit += (prop_idx == (d_hit_bids[seed.nodes[seed.size - 1]] &
                                      0xFFFFFFFFLL));

        if ((best_for_hit < best_hit_frac * seed.size)) {
            continue;
        }
        char diff_code = 0;
        bool force_dropout = false;
        if (use_dropout) {
            float4 sps[3];
            // seed 1
            sps[0] = d_sp_params[seed.nodes[seed.size - 1]];
            sps[1] = d_sp_params[seed.nodes[(seed.size - 1) / 2 + 1]];
            sps[2] = d_sp_params[seed.nodes[0]];
            float2 curv_cot_1 = estimate_params(sps);
            // seed 2
            sps[1] = d_sp_params[seed.nodes[(seed.size - 1) / 2]];
            float2 curv_cot_2 = estimate_params(sps);
            sps[0] = d_sp_params[seed.nodes[seed.size - 2]];
            // seed 3
            float2 curv_cot_3 = estimate_params(sps);
            // for low eta (higher fake rate) seeds perform a stronger cut
            if ((best_for_hit < seed.size - 1) &
                (abs(curv_cot_1.y + curv_cot_2.y + curv_cot_3.y) <
                 3.0f * tight_bid_cot_threshold) &
                (seed.size < 5)) {
                continue;
            }
            float diff[3] = {abs(curv_cot_1.x - curv_cot_2.x),
                             abs(curv_cot_2.x - curv_cot_3.x),
                             abs(curv_cot_1.x - curv_cot_3.x)};
            diff_code = 4 * (diff[0] < dcurv_cut_m) +
                        2 * (diff[1] < dcurv_cut_m) + (diff[2] < dcurv_cut_m);
            // for high pt the diff may pass dispite bad estimates
            force_dropout = abs(curv_cot_1.x + curv_cot_2.x + curv_cot_3.x) <
                            3.0f * force_dropout_max_curv_m;
            force_dropout |= (abs(curv_cot_1.y + curv_cot_2.y + curv_cot_3.y) <
                              3.0f * tight_bid_cot_threshold) &
                             diff_code == 0;
        }
        float quality = static_cast<float>(prop.x);
        // use one seed from a consistant pair/set + the inconsistant one
        // sample spacepoints from tracklet to create seeds
        // include 1st order unless either 2 or 3 are consitant with the other
        // and 1
        if (diff_code != 3 & diff_code != 6 | force_dropout) {
            seeds_device.push_back({seed.nodes[seed.size - 1],
                                    seed.nodes[(seed.size - 1) / 2 + 1],
                                    seed.nodes[0], quality});
        }
        // include 2nd order if it consistant with 1 and 3 or only 1 and 3 are
        // consistant
        if (diff_code == 1 | diff_code == 6) {
            seeds_device.push_back({seed.nodes[seed.size - 1],
                                    seed.nodes[(seed.size - 1) / 2],
                                    seed.nodes[0], quality});
        }
        // include 3rd order if it is consistant with 1 and 2 or only 1 and 2
        // are consistant or if only 2 and 3 are consistant
        if (diff_code == 2 | diff_code == 3 | diff_code == 4 | force_dropout) {
            seeds_device.push_back({seed.nodes[seed.size - 2],
                                    seed.nodes[(seed.size - 1) / 2],
                                    seed.nodes[0], quality});
        }
    }
}

}  // namespace traccc::cuda::kernels
