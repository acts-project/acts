// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointData.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Seeding/SeedFinderOrthogonalConfig.hpp"
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
class SeedFinderOrthogonal {
 public:
  /**
   * @brief Set the number of dimensions in which to embed points. This is just
   * 3 for now (phi, r, and z), but we might want to increase or decrease this
   * number in the future.
   */
  static constexpr std::size_t NDims = 3;

  /**
   * @brief The seed type used by this seeder internally.
   */
  using seed_t = Seed<external_spacepoint_t>;

  /**
   * @brief The spacepoint type used by this seeder internally.
   */
  using internal_sp_t = InternalSpacePoint<external_spacepoint_t>;

  /**
   * @brief The k-d tree type used by this seeder internally, which is
   * three-dimensional, contains internal spacepoint pointers, uses the Acts
   * scalar type for coordinates, stores its coordinates in std::arrays, and
   * has leaf size 4.
   */
  using tree_t = KDTree<NDims, internal_sp_t *, ActsScalar, std::array, 4>;

  /**
   * @brief Construct a new orthogonal seed finder.
   *
   * @param config The configuration parameters for this seed finder.
   */
  SeedFinderOrthogonal(
      const Acts::SeedFinderOrthogonalConfig<external_spacepoint_t> &config);

  /**
   * @brief Destroy the orthogonal seed finder object.
   */
  ~SeedFinderOrthogonal() = default;

  SeedFinderOrthogonal() = default;
  SeedFinderOrthogonal(const SeedFinderOrthogonal<external_spacepoint_t> &) =
      delete;
  SeedFinderOrthogonal<external_spacepoint_t> &operator=(
      const SeedFinderOrthogonal<external_spacepoint_t> &) = default;

  /**
   * @brief Perform seed finding, appending seeds to a container.
   *
   * This method performs seed finding through an orthogonal range search
   * strategy. This strategy differs from binning approaches because it selects
   * seeds _constructively_ rather than _destructively_; instead of trying a
   * large number of possible space point combinations and then rejecting many
   * of them, this algorithm tries to only consider valid seed candidates to
   * reduce combinatorics.
   *
   * In addition, this algorithm replaces the binning step used in other seed
   * finding algorithms with the construction of a k-d tree, which allows us to
   * efficiently search for space points within a given range.
   *
   * The core idea behind this algorithm is to create axis-aligned bounding
   * boxes around the region of validity for a seed candidate (be it a bottom
   * spacepoint for a given middle, a top for a given middle, a middle for a
   * given bottom, or any other combination), and then searching the detector
   * volume for points that lie inside that AABB.
   *
   * @tparam input_container_t The type of the input spacepoint container.
   * @tparam output_container_t The type of the output seed container.
   *
   * @param options frequently changing configuration (like beam position)
   * @param spacePoints The input spacepoints from which to create seeds.
   * @param out_cont The output container to write seeds to.
   * @param extract_coordinates User-defined function for extracting global position and
   * covariance of the external space point
   */
  template <typename input_container_t, typename output_container_t,
            typename callable_t>
  void createSeeds(const Acts::SeedFinderOptions &options,
                   const input_container_t &spacePoints,
                   output_container_t &out_cont,
                   callable_t &&extract_coordinates) const;

  /**
   * @brief Perform seed finding, returning a new container of seeds.
   *
   * This is a filterCandidates method for scenarios where a non-inserter API is
   * more ergonomic. In general, the inserter-based method should be preferred
   * as it is more flexible and usually more performant. For more information
   * about the seeding algorithm, please see that function.
   *
   * @tparam input_container_t The type of the input spacepoint container.
   * @param options frequently changing configuration (like beam position)
   * @param spacePoints The input spacepoints from which to create seeds.
   * @param extract_coordinates User-defined function for extracting global position and
   * covariance of the external space point
   *
   * @return A vector of seeds.
   */
  template <typename input_container_t, typename callable_t>
  std::vector<seed_t> createSeeds(const Acts::SeedFinderOptions &options,
                                  const input_container_t &spacePoints,
                                  callable_t &&extract_coordinates) const;

 private:
  /**
   * @brief Enumeration of the different dimensions in which we can apply cuts.
   */
  enum Dim { DimPhi = 0, DimR = 1, DimZ = 2 };

  /**
   * @brief Return the AABB rearch range for a given spacepoint, searching
   * upwards.
   *
   * This function calculates an axis-aligned bounding box around the volume of
   * validity for the next spacepoint in a pair, given that the lower
   * spacepoint is given. Thus, this method either takes a bottom spacepoint
   * and returns a range for the middle spacepoint, or it takes a middle
   * spacepoint and returns the range for the top spacepoint.
   *
   * @param low The lower spacepoint to find a partner for.
   *
   * @return An N-dimensional axis-aligned search range.
   */
  typename tree_t::range_t validTupleOrthoRangeLH(
      const internal_sp_t &low) const;

  /**
   * @brief Return the AABB rearch range for a given spacepoint, searching
   * downward.
   *
   * This function calculates an axis-aligned bounding box around the volume of
   * validity for the next spacepoint in a pair, given that the upper
   * spacepoint is given. Thus, this method either takes a middle spacepoint
   * and returns a range for the bottom spacepoint, or it takes a top
   * spacepoint and returns the range for the middle spacepoint.
   *
   * @param high The upper spacepoint to find a partner for.
   *
   * @return An N-dimensional axis-aligned search range.
   */
  typename tree_t::range_t validTupleOrthoRangeHL(
      const internal_sp_t &high) const;

  /**
   * @brief Check whether two spacepoints form a valid tuple.
   *
   * This method checks whether the cuts that we have for pairs of space points
   * hold.
   *
   * @warning This method checks ONLY those constraints that cannot be exactly
   * represented as bounding boxes. Thus, this method should not be used on
   * pairs of points that were not generated using a constrained spatial search
   * strategy.
   *
   * @param options frequently changing configuration (like beam position)
   * @param low The lower spacepoint.
   * @param high The upper spacepoint.
   * @param isMiddleInverted If middle spacepoint is in the negative z region
   *
   * @return True if the two points form a valid pair, false otherwise.
   */
  bool validTuple(const SeedFinderOptions &options, const internal_sp_t &low,
                  const internal_sp_t &high, bool isMiddleInverted) const;

  /**
   * @brief Create a k-d tree from a set of spacepoints.
   *
   * @param spacePoints The spacepoints to create a tree from.
   *
   * @return A k-d tree containing the given spacepoints.
   */
  tree_t createTree(const std::vector<internal_sp_t *> &spacePoints) const;

  /**
   * @brief Filter potential candidate pairs, and output seeds into an
   * iterator.
   *
   * @param options frequently changing configuration (like beam position)
   * @param middle The (singular) middle spacepoint.
   * @param bottom The (vector of) candidate bottom spacepoints.
   * @param top The (vector of) candidate top spacepoints.
   * @param seedFilterState  holds quantities used in seed filter
   * @param candidates_collector The container to write the resulting
   * seed candidates to.
   * @param spacePointData Auxiliary variables used by the seeding
   */
  void filterCandidates(
      const SeedFinderOptions &options, internal_sp_t &middle,
      std::vector<internal_sp_t *> &bottom, std::vector<internal_sp_t *> &top,
      SeedFilterState seedFilterState,
      CandidatesForMiddleSp<const InternalSpacePoint<external_spacepoint_t>>
          &candidates_collector,
      Acts::SpacePointData &spacePointData) const;

  /**
   * @brief Search for seeds starting from a given middle space point.
   *
   * @param options frequently changing configuration (like beam position)
   * @tparam NDims Number of dimensions for our spatial embedding (probably 3).
   * @tparam output_container_t Type of the output container.
   *
   * @param tree The k-d tree to use for searching.
   * @param out_cont The container write output seeds to.
   * @param middle_p The middle spacepoint to find seeds for.
   * @param spacePointData Aux data for the spacepoints
   */
  template <typename output_container_t>
  void processFromMiddleSP(const SeedFinderOptions &options, const tree_t &tree,
                           output_container_t &out_cont,
                           const typename tree_t::pair_t &middle_p,
                           Acts::SpacePointData &spacePointData) const;

  /**
   * @brief The configuration for the seeding algorithm.
   */
  Acts::SeedFinderOrthogonalConfig<external_spacepoint_t> m_config;
};
}  // namespace Acts

#include "Acts/Seeding/SeedFinderOrthogonal.ipp"
