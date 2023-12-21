// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "Acts/Utilities/GridBinFinder.hpp"
#include "Acts/Utilities/GridIterator.hpp"
#include "Acts/Utilities/Holders.hpp"
#include "Acts/Seeding/BinnedSPGroupIterator.hpp"

#include <memory>
#include <vector>

namespace Acts {

template <typename grid_t>
class BinnedGroup {
public:
#ifndef DOXYGEN
  friend BinnedGroupIterator<grid_t>;
#endif

  static constexpr std::size_t DIM = grid_t::DIM;

  /// @brief Default constructor
  BinnedGroup() = delete;

  /// brief Constructor
  BinnedGroup(std::unique_ptr<grid_t> grid,
	      std::shared_ptr<const Acts::GridBinFinder<DIM>> bottomFinder,
	      std::shared_ptr<const Acts::GridBinFinder<DIM>> topFinder,
	      std::array<std::vector<std::size_t>, DIM> navigation = std::array<std::vector<std::size_t>, DIM>());
  
  /// @brief Copy constructor
  /// @param [in] other The BinnedGroup to copy
  BinnedGroup(const BinnedGroup<grid_t>& other) = default;
  /// @brief Copy assignment
  /// @param [in] other The BinnedGroup to copy
  /// @return The copied BinnedGroup
  BinnedGroup<grid_t>& operator=(const BinnedGroup<grid_t>& other) = default;

  /// @brief Move Constructor
  /// @param [in] other The BinnedGroup to move
  BinnedGroup(BinnedGroup<grid_t>&& other) noexcept = default;
  /// @brief Move Assignment
  /// @param [in] other The BinnedGroup to move
  /// @return The moved BinnedGroup
  BinnedGroup<grid_t>& operator=(BinnedGroup<grid_t>&& other) noexcept = default;
  
  /// @brief Default destructor
  ~BinnedGroup() = default;

  /// @brief Retrieve const reference to the Grid
  /// @return Const reference to the stored grid
  const grid_t& grid() const;
  /// @brief Retrieve mutable reference to the Grid
  /// @return Mutable reference to the stored grid
  grid_t& grid();

  /// @brief Get the begin iterator
  /// @return The iterator
  Acts::BinnedGroupIterator<grid_t> begin() const;
  /// @brief Get the end iterator
  /// @return The iterator
  Acts::BinnedGroupIterator<grid_t> end() const;

  /// @brief Fill the grid
  template < typename external_spacepoint_t,
	     typename spacepoint_iterator_t,
	     typename callable_t>  
  void fill(const Acts::SeedFinderConfig<external_spacepoint_t>& config,
	    const SeedFinderOptions& options,
	    spacepoint_iterator_t spBegin, spacepoint_iterator_t spEnd,
	    callable_t&& toGlobal,
	    Acts::Extent& rRangeSPExtent);
  
private:
  /// @brief The N-dimentional grid
  std::unique_ptr<grid_t> m_grid{nullptr};
  /// @brief The Grid Bin Finder for bottom candidates
  std::shared_ptr<const Acts::GridBinFinder<DIM>> m_bottomBinFinder{nullptr};
  /// @brief The Grid Bin Finder for top candidates
  std::shared_ptr<const Acts::GridBinFinder<DIM>> m_topBinFinder{nullptr};
  /// @brief Order of bins to loop over when searching for SPs
  std::array<std::vector<std::size_t>, DIM> m_bins{};
};

} // namespace Acts

#include "Acts/Seeding/BinnedSPGroup.ipp"

namespace Acts {
  template <typename external_spacepoint_t>
  using BinnedSPGroup = Acts::BinnedGroup<Acts::SpacePointGrid<external_spacepoint_t>>;
}
