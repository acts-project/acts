// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <string>

class TFile;

namespace Acts {

class Surface;
class GeometryContext;

/// Simple payload class that can be wrapped for reading
/// and writing.
class RootMaterialMapAccessor {
 public:
  struct Config {
    /// The name of the output surface tree
    std::string folderSurfaceNameBase = "SurfaceMaterial";
    /// The name of the output volume tree
    std::string folderVolumeNameBase = "VolumeMaterial";
    /// The volume identification string
    std::string voltag = "_vol";
    /// The boundary identification string
    std::string boutag = "_bou";
    /// The layer identification string
    std::string laytag = "_lay";
    /// The approach identification string
    std::string apptag = "_app";
    /// The sensitive identification string
    std::string sentag = "_sen";
    /// The bin number tag
    std::string ntag = "n";
    /// The value tag -> binning values: AxisZ, AxisR, AxisPhi, etc.
    std::string vtag = "v";
    /// The option tag -> binning options: open, closed
    std::string otag = "o";
    /// The range min tag: min value
    std::string mintag = "min";
    /// The range max tag: max value
    std::string maxtag = "max";
    /// The thickness tag
    std::string ttag = "t";
    /// The x0 tag
    std::string x0tag = "x0";
    /// The l0 tag
    std::string l0tag = "l0";
    /// The A tag
    std::string atag = "A";
    /// The Z tag
    std::string ztag = "Z";
    /// The rho tag
    std::string rhotag = "rho";
  };

  /// @brief Constructor from config struct
  /// /// @param cfg the configuration for the accessor
  explicit RootMaterialMapAccessor(const Config& cfg) : m_cfg(cfg) {}

  /// @brief Destructor
  ~RootMaterialMapAccessor() = default;

  /// Write the material to file
  /// /// @param rFile the file to write to
  /// @param gctx the geometry context
  /// @param surface is the surface associated with the material
  void write(TFile& rFile, const GeometryContext& gctx,
             const Surface& surface) const;

 private:
  /// @brief 
  /// @param rDirectory 
  /// @param binnedMaterial 
  void writeBinnedSurfaceMaterial(TDirectory& rDirectory,
                                  const BinnedSurfaceMaterial& binnedMaterial);

  /// The configuration for the accessor
  Config m_cfg;

  /// Central store for homogeneous material
  std::unique_ptr<TTree*> m_hTree = nullptr;
  /// geometry identifier
  std::vector<uint64_t> m_hGeoId;
  /// thickness
  std::vector<float> m_ht;
  /// X0 
  std::vector<float> m_hX0;
  /// L0
  std::vector<float> m_hL0;
  /// A
  std::vector<float> m_hA;
  /// Z
  std::vector<float> m_hZ
  /// Rho
  std::vector<float> m_hRho;

};

}  // namespace Acts