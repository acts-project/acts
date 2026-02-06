// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/EventData/SpacePointProxy2.hpp"
#include "Acts/EventData/Types.hpp"

class TChain;
class TTree;

namespace ActsPlugins {
/// @addtogroup root_plugin
/// @{

/// Root space point reading and writing utility
class RootSpacePointIo {
 public:
  /// @brief sets the branch connection for writing to a file
  ///
  /// @param ttree the TTree to write to
  /// @param spacePoints the space points to write
  void connectForWrite(TTree& ttree,
                       const Acts::SpacePointContainer2& spacePoints);

  /// @brief sets the branch connection for reading from a file
  ///
  /// @param tchain the TChain to read from
  /// @param spacePoints the space points to read into
  void connectForRead(TChain& tchain,
                      const Acts::SpacePointContainer2& spacePoints);

  /// @brief Write a space point to the tree
  /// @note the caller has to do the TTree::Fill() after this call
  ///
  /// @param spacePoint the space point to write
  void write(const Acts::ConstSpacePointProxy2& spacePoint);

  /// @brief Write the space points to the tree
  ///
  /// @param spacePoints the space points to write
  /// @param ttree the TTree to write to
  void write(const Acts::SpacePointContainer2& spacePoints, TTree& ttree);

  /// @brief Read a space point from the tree
  /// @note the caller has to do the TChain::GetEntry() before this call
  ///
  /// @param spacePoint the space point to read into
  /// @param index the original index of the space point in the ROOT file
  void read(Acts::MutableSpacePointProxy2& spacePoint,
            Acts::SpacePointIndex2 index);

  /// @brief Read the space points from the tree
  ///
  /// @param tchain the TChain to read from
  /// @param spacePoints the space points to read into
  void read(TChain& tchain, Acts::SpacePointContainer2& spacePoints);

 private:
  float m_x = 0;
  float m_y = 0;
  float m_z = 0;

  float m_t = 0;

  float m_r = 0;

  float m_varZ = 0;
  float m_varR = 0;
};

/// @}
}  // namespace ActsPlugins
