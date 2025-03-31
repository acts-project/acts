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
class RootMaterialTrackPayload {
 public:
  /// @brief Constructor
  explicit RootMaterialTrackPayload(bool prePostStepInfo = false,
                                    bool surfaceInfo = false,
                                    bool volumeInfo = false,
                                    bool recalculateTotals = false)
      : m_prePostStepInfo(prePostStepInfo),
        m_surfaceInfo(surfaceInfo),
        m_volumeInfo(volumeInfo),
        m_recalculateTotals(recalculateTotals) {}

  /// @brief Destructor
  ~RootMaterialTrackPayload();

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
  bool m_prePostStepInfo = false;
  bool m_surfaceInfo = false;
  bool m_volumeInfo = false;
  bool m_recalculateTotals = false;

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
  std::vector<float>* m_stepXs = new std::vector<float>();
  /// (pre-)step y position (optional)
  std::vector<float>* m_stepYs = new std::vector<float>();
  /// (pre-)step z position (optional)
  std::vector<float>* m_stepZs = new std::vector<float>();
  /// (post-)step x position (optional)
  std::vector<float>* m_stepXe = new std::vector<float>();
  /// (post-)step y position (optional)
  std::vector<float>* m_stepYe = new std::vector<float>();
  /// (post-)step z position (optional)
  std::vector<float>* m_stepZe = new std::vector<float>();

  /// step x position
  std::vector<float>* m_stepX = new std::vector<float>();
  /// step y position
  std::vector<float>* m_stepY = new std::vector<float>();
  /// step z position
  std::vector<float>* m_stepZ = new std::vector<float>();
  /// step radial position
  std::vector<float>* m_stepR = new std::vector<float>();
  /// step x direction
  std::vector<float>* m_stepDx = new std::vector<float>();
  /// step y direction
  std::vector<float>* m_stepDy = new std::vector<float>();
  /// step z direction
  std::vector<float>* m_stepDz = new std::vector<float>();
  /// step length
  std::vector<float>* m_stepLength = new std::vector<float>();
  /// step material x0
  std::vector<float>* m_stepMatX0 = new std::vector<float>();
  /// step material l0
  std::vector<float>* m_stepMatL0 = new std::vector<float>();
  /// step material A
  std::vector<float>* m_stepMatA = new std::vector<float>();
  /// step material Z
  std::vector<float>* m_stepMatZ = new std::vector<float>();
  /// step material rho
  std::vector<float>* m_stepMatRho = new std::vector<float>();

  /// ID of the surface associated with the step
  std::vector<std::uint64_t>* m_surfaceId = new std::vector<std::uint64_t>();
  /// x position of the center of the surface associated with the step
  std::vector<float>* m_surfaceX = new std::vector<float>();
  /// y position of the center of the surface associated with the step
  std::vector<float>* m_surfaceY = new std::vector<float>();
  /// z position of the center of the surface associated with the step
  std::vector<float>* m_surfaceZ = new std::vector<float>();
  /// r position of the center of the surface associated with the step
  std::vector<float>* m_surfaceR = new std::vector<float>();
  /// distance to the surface
  std::vector<float>* m_surfaceDistance = new std::vector<float>();
  /// path correction when associating material to the given surface
  std::vector<float>* m_surfacePathCorrection = new std::vector<float>();

  /// ID of the volume associated with the step
  std::vector<std::uint64_t>* m_volumeId = new std::vector<std::uint64_t>();
};

}  // namespace Acts
