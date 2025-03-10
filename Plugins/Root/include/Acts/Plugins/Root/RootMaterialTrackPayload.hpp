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
  RootMaterialTrackPayload(bool prePostStepInfo = false,
                           bool surfaceInfo = false, bool volumeInfo = false,
                           bool collapseInteractions = false,
                           bool recalculateTotals = false)
      : m_prePostStepInfo(prePostStepInfo),
        m_surfaceInfo(surfaceInfo),
        m_volumeInfo(volumeInfo),
        m_collapseInteractions(collapseInteractions),
        m_recalculateTotals(recalculateTotals) {}

  /// @brief Destructor
  ~RootMaterialTrackPayload();

  /// @brief sets the branch connection for reading from a file
  ///
  /// @param materainChain the TChain to read the material track from
  void connectForRead(TChain& materainChain);

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
  bool m_collapseInteractions = false;
  bool m_recalculateTotals = false;

  /// Event identifier.
  std::uint32_t m_eventId = 0;

  /// The entry numbers for accessing events in increased order (there could be
  /// multiple entries corresponding to one event number)
  std::vector<long long> m_entryNumbers = {};

  /// start global x
  float m_v_x = 0;
  /// start global y
  float m_v_y = 0;
  /// start global z
  float m_v_z = 0;
  /// start global momentum x
  float m_v_px = 0;
  /// start global momentum y
  float m_v_py = 0;
  /// start global momentum z
  float m_v_pz = 0;
  /// start phi direction
  float m_v_phi = 0;
  /// start eta direction
  float m_v_eta = 0;
  /// thickness in X0/L0
  float m_tX0 = 0;
  /// thickness in X0/L0
  float m_tL0 = 0;

  /// (pre-)step x position (optional)
  std::vector<float>* m_step_sx = new std::vector<float>();
  /// (pre-)step y position (optional)
  std::vector<float>* m_step_sy = new std::vector<float>();
  /// (pre-)step z position (optional)
  std::vector<float>* m_step_sz = new std::vector<float>();
  /// (post-)step x position (optional)
  std::vector<float>* m_step_ex = new std::vector<float>();
  /// (post-)step y position (optional)
  std::vector<float>* m_step_ey = new std::vector<float>();
  /// (post-)step z position (optional)
  std::vector<float>* m_step_ez = new std::vector<float>();

  /// step x position
  std::vector<float>* m_step_x = new std::vector<float>();
  /// step y position
  std::vector<float>* m_step_y = new std::vector<float>();
  /// step z position
  std::vector<float>* m_step_z = new std::vector<float>();
  /// step radial position
  std::vector<float>* m_step_r = new std::vector<float>();
  /// step x direction
  std::vector<float>* m_step_dx = new std::vector<float>();
  /// step y direction
  std::vector<float>* m_step_dy = new std::vector<float>();
  /// step z direction
  std::vector<float>* m_step_dz = new std::vector<float>();
  /// step length
  std::vector<float>* m_step_length = new std::vector<float>();
  /// step material x0
  std::vector<float>* m_step_X0 = new std::vector<float>();
  /// step material l0
  std::vector<float>* m_step_L0 = new std::vector<float>();
  /// step material A
  std::vector<float>* m_step_A = new std::vector<float>();
  /// step material Z
  std::vector<float>* m_step_Z = new std::vector<float>();
  /// step material rho
  std::vector<float>* m_step_rho = new std::vector<float>();

  /// ID of the surface associated with the step
  std::vector<std::uint64_t>* m_sur_id = new std::vector<std::uint64_t>();
  /// x position of the center of the surface associated with the step
  std::vector<float>* m_sur_x = new std::vector<float>();
  /// y position of the center of the surface associated with the step
  std::vector<float>* m_sur_y = new std::vector<float>();
  /// z position of the center of the surface associated with the step
  std::vector<float>* m_sur_z = new std::vector<float>();
  /// r position of the center of the surface associated with the step
  std::vector<float>* m_sur_r = new std::vector<float>();
  /// distance to the surface
  std::vector<float>* m_sur_distance = new std::vector<float>();
  /// path correction when associating material to the given surface
  std::vector<float>* m_sur_pathCorrection = new std::vector<float>();

  /// ID of the volume associated with the step
  std::vector<std::uint64_t>* m_vol_id = new std::vector<std::uint64_t>();
};

}  // namespace Acts