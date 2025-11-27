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
}  // namespace Acts

namespace ActsPlugins {
/// @addtogroup root_plugin
/// @{

/// Simple payload class that can be wrapped for reading
/// and writing.
class RootMaterialTrackIo {
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
  explicit RootMaterialTrackIo(const Config& cfg) : m_cfg(cfg) {}

  /// @brief Destructor
  ~RootMaterialTrackIo() = default;

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
  void write(const Acts::GeometryContext& gctx, std::uint32_t eventNum,
             const Acts::RecordedMaterialTrack& materialTrack);

  /// @brief Read the material track from the tree
  ///
  /// @note the caller has to do the TChain::GetEntry() before this call
  ///
  /// @return the material track
  Acts::RecordedMaterialTrack read() const;

 private:
  struct MateriaSummaryPayload {
    /// start global x
    float vX = 0;
    /// start global y
    float vY = 0;
    /// start global z
    float vZ = 0;
    /// start global momentum x
    float vPx = 0;
    /// start global momentum y
    float vPy = 0;
    /// start global momentum z
    float vPz = 0;
    /// start phi direction
    float vPhi = 0;
    /// start eta direction
    float vEta = 0;
    /// thickness in X0/L0
    float tX0 = 0;
    /// thickness in X0/L0
    float tL0 = 0;
  };

  struct MaterialStepPayload {
    /// (pre-)step x position (optional)
    std::vector<float> stepXs;
    /// (pre-)step y position (optional)
    std::vector<float> stepYs;
    /// (pre-)step z position (optional)
    std::vector<float> stepZs;
    /// (post-)step x position (optional)
    std::vector<float> stepXe;
    /// (post-)step y position (optional)
    std::vector<float> stepYe;
    /// (post-)step z position (optional)
    std::vector<float> stepZe;

    /// step x position
    std::vector<float> stepX;
    std::vector<float>* stepXPtr = &stepX;
    /// step y position
    std::vector<float> stepY;
    std::vector<float>* stepYPtr = &stepY;
    /// step z position
    std::vector<float> stepZ;
    std::vector<float>* stepZPtr = &stepZ;
    /// step radial position
    std::vector<float> stepR;
    /// step x direction
    std::vector<float> stepDx;
    std::vector<float>* stepDxPtr = &stepDx;
    /// step y direction
    std::vector<float> stepDy;
    std::vector<float>* stepDyPtr = &stepDy;
    /// step z direction
    std::vector<float> stepDz;
    std::vector<float>* stepDzPtr = &stepDz;
    /// step length
    std::vector<float> stepLength;
    std::vector<float>* stepLengthPtr = &stepLength;
    /// step material x0
    std::vector<float> stepMatX0;
    std::vector<float>* stepMatX0Ptr = &stepMatX0;
    /// step material l0
    std::vector<float> stepMatL0;
    std::vector<float>* stepMatL0Ptr = &stepMatL0;
    /// step material A
    std::vector<float> stepMatA;
    std::vector<float>* stepMatAPtr = &stepMatA;
    /// step material Z
    std::vector<float> stepMatZ;
    std::vector<float>* stepMatZPtr = &stepMatZ;
    /// step material rho
    std::vector<float> stepMatRho;
    std::vector<float>* stepMatRhoPtr = &stepMatRho;
  };

  struct MaterialSurfacePayload {
    /// ID of the surface associated with the step
    std::vector<std::uint64_t> surfaceId;
    std::vector<std::uint64_t>* surfaceIdPtr = &surfaceId;
    /// x position of the center of the surface associated with the step
    std::vector<float> surfaceX;
    std::vector<float>* surfaceXPtr = &surfaceX;
    /// y position of the center of the surface associated with the step
    std::vector<float> surfaceY;
    std::vector<float>* surfaceYPtr = &surfaceY;
    /// z position of the center of the surface associated with the step
    std::vector<float> surfaceZ;
    std::vector<float>* surfaceZPtr = &surfaceZ;
    /// r position of the center of the surface associated with the step
    std::vector<float> surfaceR;
    /// distance to the surface
    std::vector<float> surfaceDistance;
    /// path correction when associating material to the given surface
    std::vector<float> surfacePathCorrection;
    std::vector<float>* surfacePathCorrectionPtr = &surfacePathCorrection;
  };

  struct MaterialVolumePayload {
    /// ID of the volume associated with the step
    std::vector<std::uint64_t> volumeId;
  };

  /// The configuration for the accessor
  Config m_cfg;

  /// Event identifier.
  std::uint32_t m_eventId = 0;

  /// The entry numbers for accessing events in increased order (there could be
  /// multiple entries corresponding to one event number)
  std::vector<long long> m_entryNumbers = {};

  /// The material summary payload
  MateriaSummaryPayload m_summaryPayload = {};

  /// The material step payload
  MaterialStepPayload m_stepPayload = {};

  /// The material surface payload
  MaterialSurfacePayload m_surfacePayload = {};

  /// The material volume payload
  MaterialVolumePayload m_volumePayload = {};
};

/// @}
}  // namespace ActsPlugins
