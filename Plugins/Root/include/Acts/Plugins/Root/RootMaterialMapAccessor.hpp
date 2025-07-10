// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/TrackingGeometryMaterial.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <map>
#include <memory>
#include <string>
#include <tuple>

class TTree;
class TFile;
class TDirectory;

namespace Acts {

class GeometryIdentifier;
class ISurfaceMaterial;
class IVolumeMaterial;
class HomogeneousSurfaceMaterial;
class MaterialSlab;
class BinnedSurfaceMaterial;

/// Simple payload class that can be wrapped for reading
/// and writing.
class RootMaterialMapAccessor {
 public:
  /// @brief Configuration for the accessor
  /// Contains the tags used for writing and reading, tag names are
  /// configuration, as they are not very likely to change.
  struct Config {
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
    /// The index tag
    std::string itag = "i";
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

  /// @brief Options for writing the material maps
  /// Folder names are optional as it allows to write more maps into one
  /// file, e.g. for the same detector with different configurations.
  struct Options {
    /// The name of the homogeneous material tree
    std::string homogeneousMaterialTreeName = "HomogeneousMaterial";
    /// The name of the indexed material tree
    std::string indexedMaterialTreeName = "IndexedMaterial";
    /// The name of the output surface tree
    std::string folderSurfaceNameBase = "SurfaceMaterial";
    /// The name of the output volume tree
    std::string folderVolumeNameBase = "VolumeMaterial";
    /// Use an indexed material tree
    bool indexedMaterial = false;
  };

  struct MaterialTreePayload {
    std::size_t index = 0;
    /// geometry identifier
    std::int64_t hGeoId = 0;
    /// thickness
    float ht = 0.0f;
    /// X0
    float hX0 = 0.0f;
    /// L0
    float hL0 = 0.0f;
    /// A
    float hA = 0.0f;
    /// Z
    float hZ = 0.0f;
    /// Rho
    float hRho = 0.0f;
  };

  /// @brief Constructor from config struct
  /// @param cfg the configuration for the accessor
  /// @param mLogger the logger to use, default is INFO level
  explicit RootMaterialMapAccessor(
      const Config& cfg,
      std::unique_ptr<const Logger> mLogger =
          getDefaultLogger("RootMaterialMapAccessor", Logging::INFO))
      : m_cfg(cfg), m_logger(std::move(mLogger)) {}

  /// @brief Destructor
  ~RootMaterialMapAccessor() = default;

  /// Write the detector maps
  /// @param rFile the file to write to
  /// @param detectorMaterial the detector material maps
  /// @param options the options for writing
  void write(TFile& rFile, const TrackingGeometryMaterial& detectorMaterial,
             const Options& options);

  /// Write the material to file
  /// @param rFile the file to write to
  /// @param geoID the geometry identifier for the surface
  /// @param surfaceMaterial is the surface associated with the material
  /// @param options the options for writing
  void write(TFile& rFile, const GeometryIdentifier& geoID,
             const ISurfaceMaterial& surfaceMaterial, const Options& options);

  /// Read the detector maps
  /// @param rFile the file to read from
  /// @param options the options for reading
  TrackingGeometryMaterial read(TFile& rFile, const Options& options);

 private:
  /// @brief Connect the homogeneous material tree for writing
  /// @param rTree the tree to connect to
  /// @param treePayload the payload to connect to the tree
  void connectForWrite(TTree& rTree, MaterialTreePayload& treePayload);

  /// @brief Connect the homogeneous material tree for writing
  /// @param rTree the tree to connect to
  /// @param treePayload the payload to connect to the tree
  void connectForRead(TTree& rTree, MaterialTreePayload& treePayload);

  /// Fill the material slab
  /// @param payload the tree payload to fill
  /// @param materialSlab the material slab to fill
  void fillMaterialSlab(MaterialTreePayload& payload,
                        const MaterialSlab& materialSlab);

  /// @brief Fill the Binned Surface material as histograms - legacy mode
  /// @param bsMaterial the binned surface material to write
  void fillBinnedSurfaceMaterial(const BinnedSurfaceMaterial& bsMaterial);

  /// @brief Fill the Binned Surface material as histograms - indexed mode
  /// @param payload the tree payload to fill
  /// @param bsMaterial the binned surface material to write
  void fillBinnedSurfaceMaterial(MaterialTreePayload& payload,
                                 const BinnedSurfaceMaterial& bsMaterial);

  /// Read the a texture Surface material
  /// @param rFile the file to read from
  /// @param tdName the name of the texture directory
  /// @param indexedMaterialTree the indexed material tree, if available
  /// @return a shared pointer to the ISurfaceMaterial
  std::shared_ptr<const ISurfaceMaterial> readTextureSurfaceMaterial(
      TFile& rFile, const std::string& tdName,
      TTree* indexedMaterialTree = nullptr);

  /// Read the a grid Surface material
  const Logger& logger() const { return *m_logger; }

  /// The configuration for the accessor
  Config m_cfg;

  /// The logger for this accessor
  std::unique_ptr<const Logger> m_logger;

  /// The homogeneous material tree
  TTree* m_hTree = nullptr;
  MaterialTreePayload m_homogenousMaterialTreePayload;

  /// The globally indexed material tree
  TTree* m_gTree = nullptr;
  MaterialTreePayload m_indexedMaterialTreePayload;
};

}  // namespace Acts
