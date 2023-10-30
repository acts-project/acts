// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// TODO: update to C++17 style

#include "Acts/EventData/SpacePointData.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Seeding/SeedFinderFTFConfig.hpp"
#include "Acts/TrackFinding/RoiDescriptor.hpp"
#include "Acts/Utilities/KDTree.hpp"

#include <array>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace Acts {

template <typename external_spacepoint_t>
struct GNN_TrigTracklet {
 public:
  GNN_TrigTracklet(std::vector<const FTF_SP<external_spacepoint_t> *> &vSP,
                   std::vector<TrigInDetTriplet<external_spacepoint_t>> &tbuf)
      : m_track(vSP), m_seeds(tbuf) {}

  std::vector<const FTF_SP<external_spacepoint_t> *> m_track;
  std::vector<TrigInDetTriplet<external_spacepoint_t>> m_seeds;
};

template <typename external_spacepoint_t>
class SeedFinderFTF {
 public:
  static constexpr std::size_t NDims = 3;

  using seed_t = Seed<external_spacepoint_t>;
  //   using internal_sp_t = InternalSpacePoint<external_spacepoint_t>;
  //   using tree_t = KDTree<NDims, internal_sp_t *, ActsScalar, std::array, 4>;

  // constructors
  SeedFinderFTF(const SeedFinderFTFConfig<external_spacepoint_t> &config,
                const TrigFTF_GNN_Geometry<external_spacepoint_t> &GNNgeo);

  ~SeedFinderFTF();  //!!! is it dangerous not to use default? got def in ipp
  SeedFinderFTF() = default;
  SeedFinderFTF(const SeedFinderFTF<external_spacepoint_t> &) = delete;
  SeedFinderFTF<external_spacepoint_t> &operator=(
      const SeedFinderFTF<external_spacepoint_t> &) = delete;

  void loadSpacePoints(
      const std::vector<FTF_SP<external_spacepoint_t>> &FTF_SP_vect);

  void createSeeds(
      const Acts::RoiDescriptor &roi,
      const Acts::TrigFTF_GNN_Geometry<external_spacepoint_t> &gnngeo);

  // create seeeds function
  template <typename input_container_t, typename output_container_t,
            typename callable_t>
  void createSeeds_old(const Acts::SeedFinderOptions &options,
                       const input_container_t &spacePoints,
                       output_container_t &out_cont,
                       callable_t &&extract_coordinates) const;

  template <typename input_container_t, typename callable_t>
  std::vector<seed_t> createSeeds_old(const Acts::SeedFinderOptions &options,
                                      const input_container_t &spacePoints,
                                      callable_t &&extract_coordinates) const;

 private:
  enum Dim { DimPhi = 0, DimR = 1, DimZ = 2 };

  // config object
  SeedFinderFTFConfig<external_spacepoint_t> m_config;

  void runGNN_TrackFinder(
      std::vector<GNN_TrigTracklet<external_spacepoint_t>> &vTracks,
      const Acts::RoiDescriptor &roi,
      const Acts::TrigFTF_GNN_Geometry<external_spacepoint_t> &gnngeo);

  // needs to be member of class so can accessed by all member functions
  TrigFTF_GNN_DataStorage<external_spacepoint_t> *m_storage;

  // for create seeds:
  std::vector<TrigInDetTriplet<external_spacepoint_t>> m_triplets;
};

}  // namespace Acts

#include "Acts/Seeding/SeedFinderFTF.ipp"
