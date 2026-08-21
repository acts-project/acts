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

#include <memory>
#include <string>
#include <vector>

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
class GridSurfaceMaterial;
}  // namespace Acts

namespace ActsPlugins {
/// @addtogroup root_plugin
/// @{

/// Simple payload class that can be wrapped for reading
/// and writing.
class RootMaterialMapIo {
 public:
  /// Configuration for the accessor
  /// Contains the tags used for writing and reading, tag names are
  /// configuration, as they are not very likely to change.
  struct Config {
    /// The volume identification string
    std::string volumePrefix = "_vol";
    /// The boundary identification string
    std::string portalPrefix = "_bou";
    /// The layer identification string
    std::string layerPrefix = "_lay";
    /// The approach identification string
    std::string passivePrefix = "_app";
    /// The sensitive identification string
    std::string sensitivePrefix = "_sen";
    /// The bin number tag
    std::string nBinsHistName = "n";
    /// The axis direction histogram name: AxisZ, AxisR, AxisPhi, etc.
    std::string axisDirHistName = "v";
    /// The axis boundary type hist name
    std::string axisBoundaryTypeHistName = "o";
    /// The range histogram name: min value
    std::string minRangeHistName = "min";
    /// The range histogram name: max value
    std::string maxRangeHistName = "max";
    /// The thickness histogram name
    std::string thicknessHistName = "t";
    /// The x0 histogram name
    std::string x0HistName = "x0";
    /// The l0 histogram name
    std::string l0HistName = "l0";
    /// The A histogram name
    std::string aHistName = "A";
    /// The Z thisogram name
    std::string zHistName = "Z";
    /// The rho thisogram name
    std::string rhoHistName = "rho";
    /// The index histogram name
    std::string indexHistName = "i";
    /// The material kind marker name: a 1-bin histogram holding a
    /// @c RootMaterialMapIo::MaterialKind value, present only on directories
    /// written by the current writer. Its presence also disambiguates
    /// @c axisBoundaryTypeHistName's encoding on read (see @c MaterialKind);
    /// its absence marks a directory written by the legacy
    /// @c BinnedSurfaceMaterial writer, which never had this concept.
    std::string materialKindHistName = "k";
  };

  /// Options for writing the material maps
  /// Folder names are optional as it allows to write more maps into one
  /// file, e.g. for the same detector with different configurations.
  struct Options {
    /// The storage backend a @c GridSurfaceMaterial is reconstructed with on
    /// read - only used as a fallback for legacy @c BinnedSurfaceMaterial
    /// data, which carries no persisted @c MaterialKind and is always
    /// upgraded to a @c GridSurfaceMaterial. Directories carrying a
    /// @c MaterialKind (i.e. written by the current writer) always
    /// reconstruct that kind instead, ignoring this option.
    enum class GridStorageKind { Direct, Indexed, GloballyIndexed };

    /// The name of the homogeneous material tree - read-only, kept for
    /// backward compatibility with files written before HomogeneousMaterial
    /// was folded into the per-directory model (see @c MaterialKind); the
    /// current writer never creates this tree.
    std::string homogeneousMaterialTreeName = "HomogeneousMaterial";
    /// The name of the indexed material tree
    std::string indexedMaterialTreeName = "IndexedMaterial";
    /// The name of the output surface tree
    std::string folderSurfaceNameBase = "SurfaceMaterial";
    /// The name of the output volume tree
    std::string folderVolumeNameBase = "VolumeMaterial";
    /// Use an indexed material tree
    bool indexedMaterial = false;
    /// The fallback storage backend for legacy BinnedSurfaceMaterial data.
    /// Defaults to GloballyIndexed so that legacy-indexed files (i.e. ones
    /// already carrying a shared material tree) are reconstructed sharing
    /// one material vector across the whole file, same as a natively written
    /// GloballyIndexed GridSurfaceMaterial would be. Legacy raw-histogram
    /// data has no shared tree to begin with, so it gets one independent
    /// vector per surface regardless of this setting.
    GridStorageKind gridStorageKind = GridStorageKind::GloballyIndexed;
  };

  /// The concrete kind of surface material a directory holds, persisted
  /// alongside its axis/index histograms so a directory is fully
  /// self-describing: a GeometryIdentifier (via the directory name), a
  /// histogram of indices into the single shared material storage, and this
  /// indication of what to reconstruct from them.
  enum class MaterialKind : int {
    Homogeneous = 0,
    Direct = 1,
    Indexed = 2,
    GloballyIndexed = 3
  };

  /// Payload structure for material tree data
  struct MaterialTreePayload {
    /// Material index
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

  /// Constructor from config struct
  /// @param cfg the configuration for the accessor
  /// @param mLogger the logger to use, default is INFO level
  explicit RootMaterialMapIo(const Config& cfg,
                             std::unique_ptr<const Acts::Logger> mLogger =
                                 Acts::getDefaultLogger("RootMaterialMapIo",
                                                        Acts::Logging::INFO))
      : m_cfg(cfg), m_logger(std::move(mLogger)) {}

  /// Write the detector maps
  /// @param rFile the file to write to
  /// @param detectorMaterial the detector material maps
  /// @param options the options for writing
  void write(TFile& rFile,
             const Acts::TrackingGeometryMaterial& detectorMaterial,
             const Options& options);

  /// Write the material to file
  /// @param rFile the file to write to
  /// @param geoID the geometry identifier for the surface
  /// @param surfaceMaterial is the surface associated with the material
  /// @param options the options for writing
  void write(TFile& rFile, const Acts::GeometryIdentifier& geoID,
             const Acts::ISurfaceMaterial& surfaceMaterial,
             const Options& options);

  /// Read the detector maps
  /// @param rFile the file to read from
  /// @param options the options for reading
  /// @return TrackingGeometryMaterial with material read from file
  Acts::TrackingGeometryMaterial read(TFile& rFile, const Options& options);

 private:
  /// Connect the homogeneous material tree for writing
  /// @param rTree the tree to connect to
  /// @param treePayload the payload to connect to the tree
  void connectForWrite(TTree& rTree, MaterialTreePayload& treePayload);

  /// Connect the homogeneous material tree for writing
  /// @param rTree the tree to connect to
  /// @param treePayload the payload to connect to the tree
  void connectForRead(TTree& rTree, MaterialTreePayload& treePayload);

  /// Fill the material slab
  /// @param payload the tree payload to fill
  /// @param materialSlab the material slab to fill
  void fillMaterialSlab(MaterialTreePayload& payload,
                        const Acts::MaterialSlab& materialSlab);

  /// Fill the Binned Surface material as histograms - legacy mode
  /// @param bsMaterial the binned surface material to write
  void fillBinnedSurfaceMaterial(const Acts::BinnedSurfaceMaterial& bsMaterial);

  /// Fill the Binned Surface material as histograms - indexed mode
  /// @param payload the tree payload to fill
  /// @param bsMaterial the binned surface material to write
  void fillBinnedSurfaceMaterial(MaterialTreePayload& payload,
                                 const Acts::BinnedSurfaceMaterial& bsMaterial);

  /// Build the per-surface directory name from a geometry identifier
  /// @param geoID the geometry identifier for the surface
  /// @param options the options for writing/reading
  /// @return the directory name
  std::string surfaceDirectoryName(const Acts::GeometryIdentifier& geoID,
                                   const Options& options) const;

  /// Write a GridSurfaceMaterial's axis description (n/v/o/min/max
  /// histograms) into the current directory
  /// @param gridMaterial the grid surface material whose axes to write
  void writeGridAxisHistograms(const Acts::GridSurfaceMaterial& gridMaterial);

  /// Write the material kind marker into the current directory
  /// @param kind the material kind to record
  void writeMaterialKind(MaterialKind kind);

  /// Fill a GridSurfaceMaterial into the current directory: uniformly for
  /// all three storage backends, an index histogram pointing into the
  /// shared indexed material tree.
  ///
  /// A GridSurfaceMaterial's on-disk index values are always absolute row
  /// numbers into the shared tree - this makes the file directly usable by
  /// tooling that operates on "the" global material vector (e.g. clustering
  /// similar slabs and re-indexing), without that tooling needing to know
  /// about any per-surface offset.
  ///
  /// This requires GloballyIndexed storage to carry the *same* material
  /// vector (checked by shared_ptr identity) for every surface in one write
  /// session - that is the actual meaning of "globally" indexed. The first
  /// GloballyIndexed surface encountered establishes that vector and writes
  /// it once; later surfaces referencing the same vector reuse its row
  /// offset for free.
  ///
  /// @param gridMaterial the grid surface material to write
  /// @throws std::invalid_argument if a GloballyIndexed surface references a
  ///         material vector different from the one already established in
  ///         this write session - use Indexed storage for a per-surface,
  ///         non-shared vector instead
  void fillGridSurfaceMaterial(const Acts::GridSurfaceMaterial& gridMaterial);

  /// Read a surface material from a surface directory
  ///
  /// Handles every on-disk shape transparently: natively written material of
  /// any MaterialKind (including Homogeneous, folded into this same
  /// per-directory model), legacy indexed BinnedSurfaceMaterial (same shape
  /// as Direct/Indexed/GloballyIndexed, different axis boundary type
  /// encoding) and legacy raw-histogram BinnedSurfaceMaterial (no shared
  /// tree involved). A directory carrying a MaterialKind always reconstructs
  /// that kind; one that does not (i.e. legacy BinnedSurfaceMaterial data)
  /// falls back to @p legacyStorageKind.
  ///
  /// @param rFile the file to read from
  /// @param tdName the name of the surface directory
  /// @param globalMaterial the fully read-in shared indexed material tree,
  ///        or nullptr if the file has none
  /// @param legacyStorageKind the storage backend to reconstruct legacy
  ///        BinnedSurfaceMaterial data with
  /// @return a shared pointer to the surface material, or nullptr if the
  ///         directory does not hold a (2D) grid or homogeneous material
  std::shared_ptr<const Acts::ISurfaceMaterial> readSurfaceMaterial(
      TFile& rFile, const std::string& tdName,
      const std::shared_ptr<std::vector<Acts::MaterialSlab>>& globalMaterial,
      Options::GridStorageKind legacyStorageKind);

  const Acts::Logger& logger() const { return *m_logger; }

  /// The configuration for the accessor
  Config m_cfg;

  /// The logger for this accessor
  std::unique_ptr<const Acts::Logger> m_logger;

  /// Payload reused when reading a legacy HomogeneousMaterial tree, see
  /// Options::homogeneousMaterialTreeName
  MaterialTreePayload m_homogenousMaterialTreePayload;

  /// The globally indexed material tree
  TTree* m_gTree = nullptr;
  MaterialTreePayload m_indexedMaterialTreePayload;

  /// Identity (the shared_ptr's stored pointer) of the one GloballyIndexed
  /// material vector written to m_gTree so far in this write session, or
  /// nullptr if none has been written yet. There is only ever one: that is
  /// the meaning of "globally" indexed, see fillGridSurfaceMaterial.
  const void* m_globalMaterialIdentity = nullptr;

  /// Row offset of m_globalMaterialIdentity's entries in m_gTree
  std::size_t m_globalMaterialOffset = 0;
};

/// @}
}  // namespace ActsPlugins
