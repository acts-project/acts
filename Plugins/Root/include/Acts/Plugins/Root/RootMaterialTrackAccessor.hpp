// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <memory>
#include <utility>
#include <vector>

class TChain;
class TTree;

namespace Acts {

class Surface;
class TrackingVolume;

/// Simple payload class that can be wrapped for reading
/// and writing.
class RootMaterialTrackAccessor {
 public:
  struct Config {
    /// Whether to store pre- and post-step information
    bool prePostStepInfo = false;
    /// Whether to store surface information
    bool surfaceInfo = false;
    /// Whether to store volume information
    bool volumeInfo = false;
    /// Whether to recalculate totals from the steps
    bool recalculateTotals = false;
  };

  /// @brief Constructor from config struct
  ///
  /// @param cfg the configuration for the accessor
  explicit RootMaterialTrackAccessor(const Config& cfg) : m_cfg(cfg) {}

  /// @brief Destructor
  ~RootMaterialTrackAccessor() = default;

  /// @brief sets the branch connection for reading from a file
  ///
  /// @param materialChain the TChain to read the material track from
  void connectForRead(TChain& materialChain);

  /// @brief sets the branch connection for writing to a file
  ///
  /// @param materialTree the TTree to write the material track to
  void connectForWrite(TTree& materialTree);

  /// @brief Write the material track to the tree
  ///
  /// @param gctx the geometry context
  /// @param eventNum the event number to write
  /// @param materialTrack the material track to be written out
  ///
  /// @note the caller has to do the TTree::Fill() after this call
  void write(const GeometryContext& gctx, std::uint32_t eventNum,
             const RecordedMaterialTrack& materialTrack);

  /// @brief Read the material track from the tree
  ///
  /// @note the caller has to do the TChain::GetEntry() before this call
  ///
  /// @return the material track
  RecordedMaterialTrack read() const;

 private:
  /// The configuration for the accessor
  Config m_cfg;

  /// Event identifier.
  std::uint32_t m_eventId = 0;

  /// The entry numbers for accessing events in increased order (there could be
  /// multiple entries corresponding to one event number)
  std::vector<long long> m_entryNumbers = {};

  /// start global x
  float m_vX = 0;
  /// start global y
  float m_vY = 0;
  /// start global z
  float m_vZ = 0;
  /// start global momentum x
  float m_vPx = 0;
  /// start global momentum y
  float m_vPy = 0;
  /// start global momentum z
  float m_vPz = 0;
  /// start phi direction
  float m_vPhi = 0;
  /// start eta direction
  float m_vEta = 0;
  /// thickness in X0/L0
  float m_tX0 = 0;
  /// thickness in X0/L0
  float m_tL0 = 0;

  /// (pre-)step x position (optional)
  std::vector<float> m_stepXs;
  /// (pre-)step y position (optional)
  std::vector<float> m_stepYs;
  /// (pre-)step z position (optional)
  std::vector<float> m_stepZs;
  /// (post-)step x position (optional)
  std::vector<float> m_stepXe;
  /// (post-)step y position (optional)
  std::vector<float> m_stepYe;
  /// (post-)step z position (optional)
  std::vector<float> m_stepZe;

  /// step x position
  std::vector<float> m_stepX;
  std::vector<float>* m_stepXPtr = &m_stepX;
  /// step y position
  std::vector<float> m_stepY;
  std::vector<float>* m_stepYPtr = &m_stepY;
  /// step z position
  std::vector<float> m_stepZ;
  std::vector<float>* m_stepZPtr = &m_stepZ;
  /// step radial position
  std::vector<float> m_stepR;
  /// step x direction
  std::vector<float> m_stepDx;
  std::vector<float>* m_stepDxPtr = &m_stepDx;
  /// step y direction
  std::vector<float> m_stepDy;
  std::vector<float>* m_stepDyPtr = &m_stepDy;
  /// step z direction
  std::vector<float> m_stepDz;
  std::vector<float>* m_stepDzPtr = &m_stepDz;
  /// step length
  std::vector<float> m_stepLength;
  std::vector<float>* m_stepLengthPtr = &m_stepLength;
  /// step material x0
  std::vector<float> m_stepMatX0;
  std::vector<float>* m_stepMatX0Ptr = &m_stepMatX0;
  /// step material l0
  std::vector<float> m_stepMatL0;
  std::vector<float>* m_stepMatL0Ptr = &m_stepMatL0;
  /// step material A
  std::vector<float> m_stepMatA;
  std::vector<float>* m_stepMatAPtr = &m_stepMatA;
  /// step material Z
  std::vector<float> m_stepMatZ;
  std::vector<float>* m_stepMatZPtr = &m_stepMatZ;
  /// step material rho
  std::vector<float> m_stepMatRho;
  std::vector<float>* m_stepMatRhoPtr = &m_stepMatRho;

  /// ID of the surface associated with the step
  std::vector<std::uint64_t> m_surfaceId;
  std::vector<std::uint64_t>* m_surfaceIdPtr = &m_surfaceId;
  /// x position of the center of the surface associated with the step
  std::vector<float> m_surfaceX;
  std::vector<float>* m_surfaceXPtr = &m_surfaceX;
  /// y position of the center of the surface associated with the step
  std::vector<float> m_surfaceY;
  std::vector<float>* m_surfaceYPtr = &m_surfaceY;
  /// z position of the center of the surface associated with the step
  std::vector<float> m_surfaceZ;
  std::vector<float>* m_surfaceZPtr = &m_surfaceZ;
  /// r position of the center of the surface associated with the step
  std::vector<float> m_surfaceR;
  /// distance to the surface
  std::vector<float> m_surfaceDistance;
  /// path correction when associating material to the given surface
  std::vector<float> m_surfacePathCorrection;
  std::vector<float>* m_surfacePathCorrectionPtr = &m_surfacePathCorrection;

  /// ID of the volume associated with the step
  std::vector<std::uint64_t> m_volumeId;
};

}  // namespace Acts
