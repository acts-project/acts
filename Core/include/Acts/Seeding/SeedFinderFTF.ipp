// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// SeedFinderFTF.ipp
// TODO: update to C++17 style

#include "Acts/Definitions/Algebra.hpp"  //for M_PI
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedFinder.hpp"
#include "Acts/Seeding/SeedFinderFTFConfig.hpp"
#include "Acts/Seeding/SeedFinderUtils.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <list>
#include <numeric>
#include <type_traits>
#include <vector>

// core so in ACTS namespace

namespace Acts {

template <typename external_spacepoint_t>
SeedFinderFTF<external_spacepoint_t>::SeedFinderFTF(
    const SeedFinderFTFConfig<external_spacepoint_t>& config,
    const TrigFTF_GNN_Geometry<external_spacepoint_t>& GNNgeo)
    : m_config(config) {
  m_storage = new TrigFTF_GNN_DataStorage(GNNgeo);
}

template <typename external_spacepoint_t>
SeedFinderFTF<external_spacepoint_t>::~SeedFinderFTF() {
  delete m_storage;

  m_storage = nullptr;
}

// define loadspace points function
template <typename external_spacepoint_t>
void SeedFinderFTF<external_spacepoint_t>::loadSpacePoints(
    const std::vector<FTF_SP<external_spacepoint_t>>& FTF_SP_vect) {
  for (const auto& FTF_sp : FTF_SP_vect) {
    bool is_Pixel = FTF_sp.isPixel();
    if (!is_Pixel) {
      continue;
    }
    m_storage->addSpacePoint(FTF_sp, (m_config.m_useClusterWidth > 0));
  }

  m_config.m_nMaxPhiSlice = 1;
  m_config.m_phiSliceWidth = 2 * M_PI / m_config.m_nMaxPhiSlice;

  m_storage->sortByPhi();

  m_storage->generatePhiIndexing(1.5 * m_config.m_phiSliceWidth);
}

template <typename external_spacepoint_t>
void SeedFinderFTF<external_spacepoint_t>::runGNN_TrackFinder(
    std::vector<GNN_TrigTracklet<external_spacepoint_t>>& vTracks,
    const Acts::RoiDescriptor& roi,
    const Acts::TrigFTF_GNN_Geometry<external_spacepoint_t>& gnngeo) {
  // long term move these to ftf finder config, then m_config. to access them
  const int MaxEdges = 2000000;

  const float cut_dphi_max = 0.012;
  const float cut_dcurv_max = 0.001;
  const float cut_tau_ratio_max = 0.007;
  const float min_z0 = -2800;  // roiDescriptor->zedMinus(); //blank for now,
                               // get from config eventually
  const float max_z0 = 2800;   // roiDescriptor->zedPlus();

  const float maxOuterRadius = 550.0;
  const float cut_zMinU =
      min_z0 +
      maxOuterRadius * roi.dzdrMinus();  // dzdr can only find =0 in athena
  const float cut_zMaxU = max_z0 + maxOuterRadius * roi.dzdrPlus();

  float m_minR_squ = 1;  // set earlier
  float m_maxCurv = 1;

  const float maxKappa_high_eta = 0.8 / m_minR_squ;
  const float maxKappa_low_eta = 0.6 / m_minR_squ;

  // 1. loop over stages

  int currentStage = 0;

  const Acts::FasTrackConnector& fastrack = *(gnngeo.fastrack());

  std::vector<Acts::TrigFTF_GNN_Edge<external_spacepoint_t>> edgeStorage;

  edgeStorage.reserve(MaxEdges);

  int nEdges = 0;

  for (std::map<int, std::vector<FasTrackConnector::LayerGroup>>::const_iterator
           it = fastrack.m_layerGroups.begin();
       it != fastrack.m_layerGroups.end(); ++it, currentStage++) {
    // loop over L1 layers for the current stage

    for (const auto& layerGroup : (*it).second) {
      unsigned int dst = layerGroup.m_dst;  // n1 : inner nodes

      const TrigFTF_GNN_Layer<external_spacepoint_t>* pL1 =
          gnngeo.getTrigFTF_GNN_LayerByKey(dst);

      if (pL1 == nullptr) {
        continue;
      }

      for (const auto& conn :
           layerGroup.m_sources) {  // loop over L2(L1) for the current stage

        unsigned int src = conn->m_src;  // n2 : the new connectors

        const TrigFTF_GNN_Layer<external_spacepoint_t>* pL2 =
            gnngeo.getTrigFTF_GNN_LayerByKey(src);

        if (pL2 == nullptr) {
          continue;
        }
        int nDstBins = pL1->m_bins.size();
        int nSrcBins = pL2->m_bins.size();

        for (int b1 = 0; b1 < nDstBins; b1++) {  // loop over bins in Layer 1

          const TrigFTF_GNN_EtaBin<external_spacepoint_t>& B1 =
              m_storage->getEtaBin(pL1->m_bins.at(b1));

          if (B1.empty()) {
            continue;
          }

          float rb1 = pL1->getMinBinRadius(b1);

          // 3. loops over source eta-bins

          for (int b2 = 0; b2 < nSrcBins; b2++) {  // loop over bins in Layer 2

            if (m_config.m_useEtaBinning && (nSrcBins + nDstBins > 2)) {
              if (conn->m_binTable[b1 + b2 * nDstBins] != 1) {
                continue;  // using precomputed LUT
              }
            }

            const TrigFTF_GNN_EtaBin<external_spacepoint_t>& B2 =
                m_storage->getEtaBin(pL2->m_bins.at(b2));

            if (B2.empty()) {
              continue;
            }

            float rb2 = pL2->getMaxBinRadius(b2);

            // calculated delta Phi for rb1 ---> rb2 extrapolation

            float deltaPhi =
                0.5f *
                m_config
                    .m_phiSliceWidth;  // the default sliding window along phi

            if (m_config.m_useEtaBinning) {
              deltaPhi = 0.001f + m_maxCurv * std::fabs(rb2 - rb1);
            }

            unsigned int first_it = 0;
            for (typename std::vector<
                     TrigFTF_GNN_Node<external_spacepoint_t>*>::const_iterator
                     n1It = B1.m_vn.begin();
                 n1It != B1.m_vn.end(); ++n1It) {  // loop over nodes in Layer 1

              TrigFTF_GNN_Node<external_spacepoint_t>* n1 = (*n1It);

              if (n1->m_in.size() >= MAX_SEG_PER_NODE) {
                continue;
              }

              float r1 = n1->m_sp_FTF.SP->r();
              float x1 = n1->m_sp_FTF.SP->x();
              float y1 = n1->m_sp_FTF.SP->y();
              float z1 = n1->m_sp_FTF.SP->z();
              float phi1 = std::atan(x1 / y1);

              float minPhi = phi1 - deltaPhi;
              float maxPhi = phi1 + deltaPhi;

              for (unsigned int n2PhiIdx = first_it;
                   n2PhiIdx < B2.m_vPhiNodes.size();
                   n2PhiIdx++) {  // sliding window over nodes in Layer 2

                float phi2 = B2.m_vPhiNodes.at(n2PhiIdx).first;

                if (phi2 < minPhi) {
                  first_it = n2PhiIdx;
                  continue;
                }
                if (phi2 > maxPhi) {
                  break;
                }

                TrigFTF_GNN_Node<external_spacepoint_t>* n2 =
                    B2.m_vn.at(B2.m_vPhiNodes.at(n2PhiIdx).second);

                if (n2->m_out.size() >= MAX_SEG_PER_NODE) {
                  continue;
                }
                if (n2->isFull()) {
                  continue;
                }

                float r2 = n2->m_sp_FTF.SP->r();

                float dr = r2 - r1;

                if (dr < m_config.m_minDeltaRadius) {
                  continue;
                }

                float z2 = n2->m_sp_FTF.SP->z();

                float dz = z2 - z1;
                float tau = dz / dr;
                float ftau = std::fabs(tau);
                if (ftau > 36.0) {
                  continue;
                }

                if (ftau < n1->m_minCutOnTau) {
                  continue;
                }
                if (ftau < n2->m_minCutOnTau) {
                  continue;
                }
                if (ftau > n1->m_maxCutOnTau) {
                  continue;
                }
                if (ftau > n2->m_maxCutOnTau) {
                  continue;
                }

                if (m_config.m_doubletFilterRZ) {
                  float z0 = z1 - r1 * tau;

                  if (z0 < min_z0 || z0 > max_z0) {
                    continue;
                  }

                  float zouter = z0 + maxOuterRadius * tau;

                  if (zouter < cut_zMinU || zouter > cut_zMaxU) {
                    continue;
                  }
                }

                float dx = n2->m_sp_FTF.SP->x() - x1;
                float dy = n2->m_sp_FTF.SP->y() - y1;

                float L2 = 1 / (dx * dx + dy * dy);

                float D =
                    (n2->m_sp_FTF.SP->y() * x1 - y1 * n2->m_sp_FTF.SP->x()) /
                    (r1 * r2);

                float kappa = D * D * L2;

                if (ftau < 4.0) {  // eta = 2.1
                  if (kappa > maxKappa_low_eta) {
                    continue;
                  }

                } else {
                  if (kappa > maxKappa_high_eta) {
                    continue;
                  }
                }

                // match edge candidate against edges incoming to n2

                float exp_eta = std::sqrt(1 + tau * tau) - tau;

                bool isGood =
                    n2->m_in.size() <=
                    2;  // we must have enough incoming edges to decide

                if (!isGood) {
                  float uat_1 = 1.0f / exp_eta;

                  for (const auto& n2_in_idx : n2->m_in) {
                    float tau2 = edgeStorage.at(n2_in_idx).m_p[0];
                    float tau_ratio = tau2 * uat_1 - 1.0f;

                    if (std::fabs(tau_ratio) > cut_tau_ratio_max) {  // bad
                                                                     // match
                      continue;
                    }
                    isGood = true;  // good match found
                    break;
                  }
                }
                if (!isGood) {
                  continue;  // no moatch found, skip creating [n1 <- n2] edge
                }

                float curv = D * std::sqrt(L2);  // signed curvature
                float dPhi2 = std::asin(curv * r2);
                float dPhi1 = std::asin(curv * r1);

                if (nEdges < MaxEdges) {
                  edgeStorage.emplace_back(n1, n2, exp_eta, curv, phi1 + dPhi1,
                                           phi2 + dPhi2);

                  n1->addIn(nEdges);
                  n2->addOut(nEdges);

                  nEdges++;
                }
              }  // loop over n2 (outer) nodes
            }    // loop over n1 (inner) nodes
          }      // loop over source eta bins
        }        // loop over dst eta bins
      }          // loop over L2(L1) layers
    }            // loop over dst layers
  }              // loop over the stages of doublet making

  std::vector<const TrigFTF_GNN_Node<external_spacepoint_t>*> vNodes;

  m_storage->getConnectingNodes(vNodes);

  if (vNodes.empty()) {
    return;
  }

  int nNodes = vNodes.size();

  for (int nodeIdx = 0; nodeIdx < nNodes; nodeIdx++) {
    const TrigFTF_GNN_Node<external_spacepoint_t>* pN = vNodes.at(nodeIdx);

    std::vector<std::pair<float, int>> in_sort, out_sort;
    in_sort.resize(pN->m_in.size());
    out_sort.resize(pN->m_out.size());

    for (int inIdx = 0; inIdx < static_cast<int>(pN->m_in.size()); inIdx++) {
      int inEdgeIdx = pN->m_in.at(inIdx);
      Acts::TrigFTF_GNN_Edge<external_spacepoint_t>* pS =
          &(edgeStorage.at(inEdgeIdx));
      in_sort[inIdx].second = inEdgeIdx;
      in_sort[inIdx].first = pS->m_p[0];
    }
    for (int outIdx = 0; outIdx < static_cast<int>(pN->m_out.size());
         outIdx++) {
      int outEdgeIdx = pN->m_out.at(outIdx);
      Acts::TrigFTF_GNN_Edge<external_spacepoint_t>* pS =
          &(edgeStorage.at(outEdgeIdx));
      out_sort[outIdx].second = outEdgeIdx;
      out_sort[outIdx].first = pS->m_p[0];
    }

    std::sort(in_sort.begin(), in_sort.end());
    std::sort(out_sort.begin(), out_sort.end());

    unsigned int last_out = 0;

    for (unsigned int in_idx = 0; in_idx < in_sort.size();
         in_idx++) {  // loop over incoming edges

      int inEdgeIdx = in_sort[in_idx].second;

      Acts::TrigFTF_GNN_Edge<external_spacepoint_t>* pS =
          &(edgeStorage.at(inEdgeIdx));

      pS->m_nNei = 0;
      float tau1 = pS->m_p[0];
      float uat_1 = 1.0f / tau1;
      float curv1 = pS->m_p[1];
      float Phi1 = pS->m_p[2];

      for (unsigned int out_idx = last_out; out_idx < out_sort.size();
           out_idx++) {
        int outEdgeIdx = out_sort[out_idx].second;

        Acts::TrigFTF_GNN_Edge<external_spacepoint_t>* pNS =
            &(edgeStorage.at(outEdgeIdx));

        float tau2 = pNS->m_p[0];
        float tau_ratio = tau2 * uat_1 - 1.0f;

        if (tau_ratio < -cut_tau_ratio_max) {
          last_out = out_idx;
          continue;
        }
        if (tau_ratio > cut_tau_ratio_max) {
          break;
        }

        float dPhi = pNS->m_p[3] - Phi1;

        if (dPhi < -M_PI) {
          dPhi += 2 * M_PI;
        } else if (dPhi > M_PI) {
          dPhi -= 2 * M_PI;
        }

        if (dPhi < -cut_dphi_max || dPhi > cut_dphi_max) {
          continue;
        }

        float curv2 = pNS->m_p[1];
        float dcurv = curv2 - curv1;

        if (dcurv < -cut_dcurv_max || dcurv > cut_dcurv_max) {
          continue;
        }

        pS->m_vNei[pS->m_nNei++] = outEdgeIdx;
        if (pS->m_nNei >= N_SEG_CONNS) {
          break;
        }
      }
    }
  }

  const int maxIter = 15;

  int maxLevel = 0;

  int iter = 0;

  std::vector<Acts::TrigFTF_GNN_Edge<external_spacepoint_t>*> v_old;

  for (int edgeIndex = 0; edgeIndex < nEdges; edgeIndex++) {
    Acts::TrigFTF_GNN_Edge<external_spacepoint_t>* pS =
        &(edgeStorage.at(edgeIndex));
    if (pS->m_nNei == 0) {
      continue;
    }
    v_old.push_back(pS);  // TO-DO: increment level for segments as they already
                          // have at least one neighbour
  }

  for (; iter < maxIter; iter++) {
    // generate proposals
    std::vector<Acts::TrigFTF_GNN_Edge<external_spacepoint_t>*> v_new;
    v_new.clear();

    for (auto pS : v_old) {
      int next_level = pS->m_level;

      for (int nIdx = 0; nIdx < pS->m_nNei; nIdx++) {
        unsigned int nextEdgeIdx = pS->m_vNei[nIdx];

        Acts::TrigFTF_GNN_Edge<external_spacepoint_t>* pN =
            &(edgeStorage.at(nextEdgeIdx));

        if (pS->m_level == pN->m_level) {
          next_level = pS->m_level + 1;
          v_new.push_back(pS);
          break;
        }
      }

      pS->m_next = next_level;  // proposal
    }
    // update

    int nChanges = 0;

    for (auto pS : v_new) {
      if (pS->m_next != pS->m_level) {
        nChanges++;
        pS->m_level = pS->m_next;
        if (maxLevel < pS->m_level) {
          maxLevel = pS->m_level;
        }
      }
    }

    if (nChanges == 0) {
      break;
    }

    v_old = std::move(v_new);
    v_new.clear();
  }

  int minLevel = 3;  // a triplet + 2 confirmation

  std::vector<Acts::TrigFTF_GNN_Edge<external_spacepoint_t>*> vSeeds;

  vSeeds.reserve(MaxEdges / 2);

  for (int edgeIndex = 0; edgeIndex < nEdges; edgeIndex++) {
    Acts::TrigFTF_GNN_Edge<external_spacepoint_t>* pS =
        &(edgeStorage.at(edgeIndex));

    if (pS->m_level < minLevel) {
      continue;
    }

    vSeeds.push_back(pS);
  }

  m_triplets.clear();

  std::sort(
      vSeeds.begin(), vSeeds.end(),
      typename Acts::TrigFTF_GNN_Edge<external_spacepoint_t>::CompareLevel());

  if (vSeeds.empty()) {
    return;
  }

  // backtracking

  TrigFTF_GNN_TrackingFilter<external_spacepoint_t> tFilter(
      m_config.m_layerGeometry, edgeStorage);

  for (auto pS : vSeeds) {
    if (pS->m_level == -1) {
      continue;
    }

    TrigFTF_GNN_EdgeState<external_spacepoint_t> rs(false);

    tFilter.followTrack(pS, rs);

    if (!rs.m_initialized) {
      continue;
    }

    if (static_cast<int>(rs.m_vs.size()) < minLevel) {
      continue;
    }

    std::vector<const FTF_SP<external_spacepoint_t>*> vSP;

    for (typename std::vector<Acts::TrigFTF_GNN_Edge<external_spacepoint_t>*>::
             reverse_iterator sIt = rs.m_vs.rbegin();
         sIt != rs.m_vs.rend(); ++sIt) {
      (*sIt)->m_level = -1;  // mark as collected

      if (sIt == rs.m_vs.rbegin()) {
        vSP.push_back(&(*sIt)->m_n1->m_sp_FTF);
      }
      vSP.push_back(&(*sIt)->m_n2->m_sp_FTF);
    }

    if (vSP.size() < 3) {
      continue;
    }

    // making triplets

    unsigned int nTriplets = 0;

    std::vector<TrigInDetTriplet<external_spacepoint_t>> output;

    for (unsigned int idx_m = 1; idx_m < vSP.size() - 1; idx_m++) {
      const FTF_SP<external_spacepoint_t>& spM = *vSP.at(idx_m);
      const double pS_r = spM.SP->r();
      const double pS_x = spM.SP->x();
      const double pS_y = spM.SP->y();
      const double cosA = pS_x / pS_r;
      const double sinA = pS_y / pS_r;

      for (unsigned int idx_o = idx_m + 1; idx_o < vSP.size(); idx_o++) {
        const FTF_SP<external_spacepoint_t>& spO = *vSP.at(idx_o);

        double dx = spO.SP->x() - pS_x;
        double dy = spO.SP->y() - pS_y;
        double R2inv = 1.0 / (dx * dx + dy * dy);
        double xn = dx * cosA + dy * sinA;
        double yn = -dx * sinA + dy * cosA;

        const double uo = xn * R2inv;
        const double vo = yn * R2inv;

        for (unsigned int idx_i = 0; idx_i < idx_m; idx_i++) {
          const FTF_SP<external_spacepoint_t>& spI = *vSP.at(idx_i);

          dx = spI.SP->x() - pS_x;
          dy = spI.SP->y() - pS_y;
          R2inv = 1.0 / (dx * dx + dy * dy);

          xn = dx * cosA + dy * sinA;
          yn = -dx * sinA + dy * cosA;

          const double ui = xn * R2inv;
          const double vi = yn * R2inv;

          // 1. pT estimate

          const double du = uo - ui;
          if (du == 0.0) {
            continue;
          }
          const double A = (vo - vi) / du;
          const double B = vi - A * ui;
          const double R_squ = (1 + A * A) / (B * B);

          if (R_squ < m_minR_squ) {
            continue;
          }

          // 2. d0 cut

          const double fabs_d0 = std::abs(pS_r * (B * pS_r - A));

          if (fabs_d0 > m_config.m_tripletD0Max) {
            continue;
          }

          // 3. phi0 cut

          // if (!roi.isFullscan()) {
          //   const double uc = 2 * B * pS_r - A;
          //   // const double phi0 = std::atan2(sinA - uc * cosA, cosA + uc *
          //   sinA);
          //   //           if ( !RoiUtil::containsPhi( *roiDescriptor, phi0 ) )
          //   {
          //   //             continue;
          //   //           }
          // }

          // 4. add new triplet

          const double Q = fabs_d0 * fabs_d0;

          output.emplace_back(spI, spM, spO, Q);

          nTriplets++;

          if (nTriplets >= m_config.m_maxTripletBufferLength) {
            break;
          }
        }
        if (nTriplets >= m_config.m_maxTripletBufferLength) {
          break;
        }
      }
      if (nTriplets >= m_config.m_maxTripletBufferLength) {
        break;
      }
    }

    if (output.empty()) {
      continue;
    }

    vTracks.emplace_back(vSP, output);
  }
}

template <typename external_spacepoint_t>
void SeedFinderFTF<external_spacepoint_t>::createSeeds(
    const Acts::RoiDescriptor& roi,
    const Acts::TrigFTF_GNN_Geometry<external_spacepoint_t>& gnngeo) {
  std::vector<GNN_TrigTracklet<external_spacepoint_t>>
      vTracks;  // make empty vector

  vTracks.reserve(5000);

  runGNN_TrackFinder(vTracks, roi, gnngeo);  // returns filled vector

  if (vTracks.empty()) {
    return;
  }

  m_triplets.clear();  // member of class , saying not declared, maybe public?

  for (auto& track : vTracks) {
    for (auto& seed : track.m_seeds) {  // access mmeber of GNN_TrigTracklet

      float newQ = seed.Q();  // function of TrigInDetTriplet
      if (m_config.m_LRTmode) {
        // In LRT mode penalize pixels in Triplets
        if (seed.s1().isPixel()) {
          newQ += 1000;  // functions of TrigSiSpacePointBase
        }
        if (seed.s2().isPixel()) {
          newQ += 1000;
        }
        if (seed.s3().isPixel()) {
          newQ += 1000;
        }
      } else {
        // In normal (non LRT) mode penalise SSS by 1000, PSS (if enabled) and
        // PPS by 10000
        if (seed.s3().isSCT()) {
          newQ += seed.s1().isSCT() ? 1000.0 : 10000.0;
        }
      }
      seed.Q(newQ);
      m_triplets.emplace_back(seed);
    }
  }
  vTracks.clear();
}

// // still to be developed
template <typename external_spacepoint_t>
template <typename input_container_t, typename output_container_t,
          typename callable_t>
void SeedFinderFTF<external_spacepoint_t>::createSeeds_old(
    const Acts::SeedFinderOptions& /*options*/,
    const input_container_t& /*spacePoints*/, output_container_t& /*out_cont*/,
    callable_t&& /*extract_coordinates*/) const {}

template <typename external_spacepoint_t>
template <typename input_container_t, typename callable_t>
std::vector<Seed<external_spacepoint_t>>
SeedFinderFTF<external_spacepoint_t>::createSeeds_old(
    const Acts::SeedFinderOptions& options,
    const input_container_t& spacePoints,
    callable_t&& extract_coordinates) const {
  std::vector<seed_t> r;
  createSeeds_old(options, spacePoints, r,
                  std::forward<callable_t>(extract_coordinates));
  return r;
}

}  // namespace Acts
