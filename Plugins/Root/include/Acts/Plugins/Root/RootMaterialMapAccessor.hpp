// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <string>

#include <TTree.h>

class TFile;
class TDirectory;

namespace Acts {

class Surface;
class GeometryContext;
class HomogeneousSurfaceMaterial;
class BinnedSurfaceMaterial;

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

  struct MaterialTreePayload {
    std::vector<uint64_t> hGeoId;
    std::vector<uint64_t>* hGeoIdPtr = &hGeoId;
    /// thickness
    std::vector<float> ht;
    std::vector<float>* htPtr = &ht;
    /// X0
    std::vector<float> hX0;
    std::vector<float>* hX0Ptr = &hX0;
    /// L0
    std::vector<float> hL0;
    std::vector<float>* hL0Ptr = &hL0;
    /// A
    std::vector<float> hA;
    std::vector<float>* hAPtr = &hA;
    /// Z
    std::vector<float> hZ;
    std::vector<float>* hZPtr = &hZ;
    /// Rho
    std::vector<float> hRho;
    std::vector<float>* hRhoPtr = &hRho;
  };

  /// @brief Constructor from config struct
  /// @param cfg the configuration for the accessor
  explicit RootMaterialMapAccessor(const Config& cfg) : cfg(cfg) {}

  /// @brief Destructor
  ~RootMaterialMapAccessor() = default;

  /// Write the material to file
  /// @param rFile the file to write to
  /// @param gctx the geometry context
  /// @param surface is the surface associated with the material
  void write(TFile& rFile, const GeometryContext& gctx, const Surface& surface);

 private:
  /// @brief Write the homogeneous material to the file
  /// @param homogeneousMaterial the homogeneous material to write
  void writeHomogeneousMaterial(
      const HomogeneousSurfaceMaterial& homogeneousMaterial);

  /// @brief Connect the homogeneous material tree for writing
  /// @param rTree the tree to connect to
  /// @param treePayload the payload to connect to the tree
  void connectForWrite(TTree& rTree, MaterialTreePayload& treePayload);


  /// @brief Connect the homogeneous material tree for writing
  /// @param rTree the tree to connect to
  /// @param treePayload the payload to connect to the tree
  void connectForRead(const TTree& rTree, MaterialTreePayload& treePayload);

  /// @brief
  /// @param rDirectory
  /// @param binnedMaterial
  void writeBinnedSurfaceMaterial(TDirectory& rDirectory,
                                  const BinnedSurfaceMaterial& binnedMaterial);

  /// The configuration for the accessor
  Config m_cfg;

  /// The homogeneous material tree
  std::unique_ptr<TTree> m_hTree = nullptr;
  MaterialTreePayload m_homogenousMaterialTreePayload;

  /// The globally indexed material tree
  std::unique_ptr<TTree> m_gTree = nullptr;
  MaterialTreePayload m_globallyIndexedMaterialTreePayload;


};

}  // namespace Acts