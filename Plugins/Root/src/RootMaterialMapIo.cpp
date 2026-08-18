// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Root/RootMaterialMapIo.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/GridSurfaceMaterial.hpp"
#include "Acts/Material/GridSurfaceMaterialFactory.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/AxisSpec.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/MultiAxisSpec.hpp"

#include <algorithm>
#include <array>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>

#include <TFile.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TKey.h>
#include <TList.h>
#include <TObject.h>
#include <TTree.h>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/finder.hpp>
#include <boost/algorithm/string/iter_find.hpp>

using namespace Acts;

namespace ActsPlugins {

namespace {

/// @brief Legacy `BinningOption` (open = 0, closed = 1), kept only to
/// translate axis boundary types read back from files written before the
/// GridSurfaceMaterial migration - it has no `Bound` counterpart.
enum LegacyBinningOption { legacyOpen = 0, legacyClosed = 1 };

/// @brief Translate a legacy `BinningOption` value into the corresponding
/// `AxisBoundaryType` - `Bound` never appears in legacy files.
AxisBoundaryType boundaryTypeFromLegacy(int legacyOption) {
  switch (legacyOption) {
    case LegacyBinningOption::legacyOpen:
      return AxisBoundaryType::Open;
    case LegacyBinningOption::legacyClosed:
      return AxisBoundaryType::Closed;
    default:
      throw std::invalid_argument(
          "RootMaterialMapIo: unknown legacy binning option " +
          std::to_string(legacyOption));
  }
}

using GridStorageKind = RootMaterialMapIo::Options::GridStorageKind;
using MaterialKind = RootMaterialMapIo::MaterialKind;

/// @brief The MaterialKind matching a GridSurfaceMaterial's storage backend
MaterialKind materialKindOf(const GridSurfaceMaterial::Storage& storage) {
  return std::visit(
      [](const auto& s) {
        using T = std::decay_t<decltype(s)>;
        if constexpr (std::is_same_v<T, GridSurfaceMaterial::Direct>) {
          return MaterialKind::Direct;
        } else if constexpr (std::is_same_v<T, GridSurfaceMaterial::Indexed>) {
          return MaterialKind::Indexed;
        } else {
          return MaterialKind::GloballyIndexed;
        }
      },
      storage);
}

/// @brief Build a GridSurfaceMaterial from a dense (regular-bin only)
/// material matrix - used for the legacy raw-histogram on-disk shape, which
/// carries no natural index/dedup structure.
std::shared_ptr<const GridSurfaceMaterial> buildFromMaterialPayload(
    const IAxis& axis0, const IAxis& axis1,
    const std::vector<std::vector<MaterialSlab>>& payload,
    GridStorageKind storageKind) {
  if (storageKind == GridStorageKind::Direct) {
    return GridSurfaceMaterialFactory::createDirect(axis0, axis1, payload);
  }
  // Indexed / GloballyIndexed: one unique slab per bin, in bin order
  std::vector<MaterialSlab> material;
  std::vector<std::vector<std::size_t>> indices(payload.size());
  for (std::size_t i0 = 0; i0 < payload.size(); ++i0) {
    indices[i0].resize(payload[i0].size());
    for (std::size_t i1 = 0; i1 < payload[i0].size(); ++i1) {
      indices[i0][i1] = material.size();
      material.push_back(payload[i0][i1]);
    }
  }
  if (storageKind == GridStorageKind::Indexed) {
    return GridSurfaceMaterialFactory::createIndexed(
        axis0, axis1, std::move(material), indices);
  }
  auto sharedMaterial =
      std::make_shared<std::vector<MaterialSlab>>(std::move(material));
  return GridSurfaceMaterialFactory::createGloballyIndexed(
      axis0, axis1, sharedMaterial, indices);
}

/// @brief Build a GridSurfaceMaterial from indices into an already read-in
/// shared material vector - used for both the new and the legacy indexed
/// on-disk shape, which are byte-for-byte identical.
std::shared_ptr<const GridSurfaceMaterial> buildFromIndexPayload(
    const IAxis& axis0, const IAxis& axis1,
    const std::vector<std::vector<std::size_t>>& indexPayload,
    const std::shared_ptr<std::vector<MaterialSlab>>& globalMaterial,
    GridStorageKind storageKind) {
  if (storageKind == GridStorageKind::GloballyIndexed) {
    return GridSurfaceMaterialFactory::createGloballyIndexed(
        axis0, axis1, globalMaterial, indexPayload);
  }
  if (storageKind == GridStorageKind::Direct) {
    std::vector<std::vector<MaterialSlab>> payload(indexPayload.size());
    for (std::size_t i0 = 0; i0 < indexPayload.size(); ++i0) {
      payload[i0].resize(indexPayload[i0].size());
      for (std::size_t i1 = 0; i1 < indexPayload[i0].size(); ++i1) {
        payload[i0][i1] = globalMaterial->at(indexPayload[i0][i1]);
      }
    }
    return GridSurfaceMaterialFactory::createDirect(axis0, axis1, payload);
  }
  // Indexed: compact the referenced rows into a locally owned vector,
  // remapping indices
  std::unordered_map<std::size_t, std::size_t> remap;
  std::vector<MaterialSlab> material;
  std::vector<std::vector<std::size_t>> localIndices(indexPayload.size());
  for (std::size_t i0 = 0; i0 < indexPayload.size(); ++i0) {
    localIndices[i0].resize(indexPayload[i0].size());
    for (std::size_t i1 = 0; i1 < indexPayload[i0].size(); ++i1) {
      std::size_t globalRow = indexPayload[i0][i1];
      auto [it, inserted] = remap.try_emplace(globalRow, material.size());
      if (inserted) {
        material.push_back(globalMaterial->at(globalRow));
      }
      localIndices[i0][i1] = it->second;
    }
  }
  return GridSurfaceMaterialFactory::createIndexed(
      axis0, axis1, std::move(material), localIndices);
}

}  // namespace

void RootMaterialMapIo::write(TFile& rFile, const GeometryIdentifier& geoID,
                              const ISurfaceMaterial& surfaceMaterial,
                              const Options& options) {
  /// Change to the file
  rFile.cd();

  // Homogeneous surface material writing - folded into the same
  // per-directory, shared-tree model as grid material: a trivial 1x1 Direct
  // grid whose MaterialKind marks it as Homogeneous, so the picture stays
  // consistent (GeometryIdentifier + index histogram + kind) regardless of
  // what kind of material a given surface carries or the order surfaces are
  // written in.
  auto homogeneousMaterial =
      dynamic_cast<const HomogeneousSurfaceMaterial*>(&surfaceMaterial);
  if (homogeneousMaterial != nullptr) {
    std::string tdName = surfaceDirectoryName(geoID, options);
    rFile.mkdir(tdName.c_str());
    rFile.cd(tdName.c_str());

    // The axes are never used to reconstruct a HomogeneousSurfaceMaterial on
    // read (only the MaterialKind is), so their range is an arbitrary but
    // valid placeholder. Only the one regular bin ever reaches disk (see
    // fillGridSurfaceMaterial), the under-/overflow entries just exist to
    // satisfy GridSurfaceMaterial's storage-size invariant.
    MultiAxisSpec2D binning(
        {AxisSpec::Equidistant(1, 0., 1., AxisBoundaryType::Bound),
         AxisSpec::Equidistant(1, 0., 1., AxisBoundaryType::Bound)});
    auto trivialMultiAxis = binning.buildMultiAxis();
    GridSurfaceMaterial::Direct storage(trivialMultiAxis->getNTotalBins(true),
                                        homogeneousMaterial->materialSlab());
    GridSurfaceMaterial temporaryGrid(std::move(binning), std::move(storage));

    writeGridAxisHistograms(temporaryGrid);
    writeMaterialKind(MaterialKind::Homogeneous);

    if (m_gTree == nullptr) {
      rFile.cd();
      m_gTree = new TTree(options.indexedMaterialTreeName.c_str(),
                          "Indexed Material Tree");
      connectForWrite(*m_gTree, m_indexedMaterialTreePayload);
      rFile.cd(tdName.c_str());
    }
    fillGridSurfaceMaterial(temporaryGrid);
    return;
  }

  // Binned surface material writing - legacy format, kept for callers still
  // producing BinnedSurfaceMaterial in memory
  auto bsMaterial =
      dynamic_cast<const BinnedSurfaceMaterial*>(&surfaceMaterial);
  if (bsMaterial != nullptr) {
    std::string tdName = surfaceDirectoryName(geoID, options);
    // create a new directory
    rFile.mkdir(tdName.c_str());
    rFile.cd(tdName.c_str());

    // Boundary condistions
    // Get the binning data
    auto& binningData = bsMaterial->binUtility().binningData();
    // 1-D or 2-D maps
    auto bins = static_cast<int>(binningData.size());
    auto fBins = static_cast<float>(bins);

    // The bin number information
    TH1F n(m_cfg.nBinsHistName.c_str(), "bins; bin", bins, -0.5, fBins - 0.5);

    // The binning value information
    TH1F v(m_cfg.axisDirHistName.c_str(), "binning values; bin", bins, -0.5,
           fBins - 0.5);

    // The binning option information
    TH1F o(m_cfg.axisBoundaryTypeHistName.c_str(), "binning options; bin", bins,
           -0.5, fBins - 0.5);

    // The binning option information - range min
    TH1F rmin(m_cfg.minRangeHistName.c_str(), "min; bin", bins, -0.5,
              fBins - 0.5);

    // The binning option information - range max
    TH1F rmax(m_cfg.maxRangeHistName.c_str(), "max; bin", bins, -0.5,
              fBins - 0.5);

    // Now fill the histogram content
    for (auto [b, bData] : enumerate(binningData)) {
      // Fill: nbins, value, option, min, max
      n.SetBinContent(static_cast<int>(b) + 1, static_cast<int>(bData.bins()));
      v.SetBinContent(static_cast<int>(b) + 1,
                      static_cast<int>(bData.binvalue));
      o.SetBinContent(static_cast<int>(b) + 1, static_cast<int>(bData.option));
      rmin.SetBinContent(static_cast<int>(b) + 1, bData.min);
      rmax.SetBinContent(static_cast<int>(b) + 1, bData.max);
    }
    n.Write();
    v.Write();
    o.Write();
    rmin.Write();
    rmax.Write();

    // If compressed writing is not enabled, write the binned surface material
    // as histograms
    if (!options.indexedMaterial) {
      fillBinnedSurfaceMaterial(*bsMaterial);
      return;
    }

    // Otherwise, write the binned surface material into the TTree
    if (m_gTree == nullptr) {
      // Back to file level
      rFile.cd();
      m_gTree = new TTree(options.indexedMaterialTreeName.c_str(),
                          "Indexed Material Tree");
      connectForWrite(*m_gTree, m_indexedMaterialTreePayload);
      // Back to the directory
      rFile.cd(tdName.c_str());
    }
    fillBinnedSurfaceMaterial(m_indexedMaterialTreePayload, *bsMaterial);
    return;
  }

  // Grid surface material writing - always through the shared indexed
  // material tree, regardless of the in-memory storage backend
  auto gridMaterial =
      dynamic_cast<const GridSurfaceMaterial*>(&surfaceMaterial);
  if (gridMaterial != nullptr) {
    std::string tdName = surfaceDirectoryName(geoID, options);
    rFile.mkdir(tdName.c_str());
    rFile.cd(tdName.c_str());

    writeGridAxisHistograms(*gridMaterial);
    writeMaterialKind(materialKindOf(gridMaterial->storage()));

    if (m_gTree == nullptr) {
      rFile.cd();
      m_gTree = new TTree(options.indexedMaterialTreeName.c_str(),
                          "Indexed Material Tree");
      connectForWrite(*m_gTree, m_indexedMaterialTreePayload);
      rFile.cd(tdName.c_str());
    }
    fillGridSurfaceMaterial(*gridMaterial);
    return;
  }
}

void RootMaterialMapIo::write(TFile& rFile,
                              const TrackingGeometryMaterial& detectorMaterial,
                              const Options& options) {
  const auto& [surfaceMaterials, volumeMaterials] = detectorMaterial;
  for (const auto& [geoID, sMaterial] : surfaceMaterials) {
    write(rFile, geoID, *sMaterial, options);
  }
  if (m_gTree != nullptr) {
    m_gTree->Write();
  }
}

std::string RootMaterialMapIo::surfaceDirectoryName(
    const GeometryIdentifier& geoID, const Options& options) const {
  std::string tdName = options.folderSurfaceNameBase;
  tdName += m_cfg.volumePrefix + std::to_string(geoID.volume());
  tdName += m_cfg.portalPrefix + std::to_string(geoID.boundary());
  tdName += m_cfg.layerPrefix + std::to_string(geoID.layer());
  tdName += m_cfg.passivePrefix + std::to_string(geoID.approach());
  tdName += m_cfg.sensitivePrefix + std::to_string(geoID.sensitive());
  return tdName;
}

void RootMaterialMapIo::connectForWrite(TTree& rTree,
                                        MaterialTreePayload& treePayload) {
  rTree.Branch(m_cfg.thicknessHistName.c_str(), &treePayload.ht);
  rTree.Branch(m_cfg.x0HistName.c_str(), &treePayload.hX0);
  rTree.Branch(m_cfg.l0HistName.c_str(), &treePayload.hL0);
  rTree.Branch(m_cfg.aHistName.c_str(), &treePayload.hA);
  rTree.Branch(m_cfg.zHistName.c_str(), &treePayload.hZ);
  rTree.Branch(m_cfg.rhoHistName.c_str(), &treePayload.hRho);
}

void RootMaterialMapIo::connectForRead(TTree& rTree,
                                       MaterialTreePayload& treePayload) {
  if (&treePayload == &m_homogenousMaterialTreePayload) {
    rTree.SetBranchAddress("hGeoId", &treePayload.hGeoId);
  }
  rTree.SetBranchAddress(m_cfg.thicknessHistName.c_str(), &treePayload.ht);
  rTree.SetBranchAddress(m_cfg.x0HistName.c_str(), &treePayload.hX0);
  rTree.SetBranchAddress(m_cfg.l0HistName.c_str(), &treePayload.hL0);
  rTree.SetBranchAddress(m_cfg.aHistName.c_str(), &treePayload.hA);
  rTree.SetBranchAddress(m_cfg.zHistName.c_str(), &treePayload.hZ);
  rTree.SetBranchAddress(m_cfg.rhoHistName.c_str(), &treePayload.hRho);
}

void RootMaterialMapIo::fillMaterialSlab(MaterialTreePayload& payload,
                                         const MaterialSlab& materialSlab) {
  payload.ht = materialSlab.thickness();
  payload.hX0 = materialSlab.material().X0();
  payload.hL0 = materialSlab.material().L0();
  payload.hA = materialSlab.material().Ar();
  payload.hZ = materialSlab.material().Z();
  payload.hRho = materialSlab.material().massDensity();
}

void RootMaterialMapIo::fillBinnedSurfaceMaterial(
    const BinnedSurfaceMaterial& bsMaterial) {
  auto bins0 = static_cast<int>(bsMaterial.binUtility().bins(0));
  auto bins1 = static_cast<int>(bsMaterial.binUtility().bins(1));
  auto fBins0 = static_cast<float>(bins0);
  auto fBins1 = static_cast<float>(bins1);

  TH2F t(m_cfg.thicknessHistName.c_str(), "thickness [mm] ;b0 ;b1", bins0, -0.5,
         fBins0 - 0.5, bins1, -0.5, fBins1 - 0.5);
  TH2F x0(m_cfg.x0HistName.c_str(), "X_{0} [mm] ;b0 ;b1", bins0, -0.5,
          fBins0 - 0.5, bins1, -0.5, fBins1 - 0.5);
  TH2F l0(m_cfg.l0HistName.c_str(), "#Lambda_{0} [mm] ;b0 ;b1", bins0, -0.5,
          fBins0 - 0.5, bins1, -0.5, fBins1 - 0.5);
  TH2F A(m_cfg.aHistName.c_str(), "X_{0} [mm] ;b0 ;b1", bins0, -0.5,
         fBins0 - 0.5, bins1, -0.5, fBins1 - 0.5);
  TH2F Z(m_cfg.zHistName.c_str(), "#Lambda_{0} [mm] ;b0 ;b1", bins0, -0.5,
         fBins0 - 0.5, bins1, -0.5, fBins1 - 0.5);
  TH2F rho(m_cfg.rhoHistName.c_str(), "#rho [g/mm^3] ;b0 ;b1", bins0, -0.5,
           fBins0 - 0.5, bins1, -0.5, fBins1 - 0.5);

  // Loop over the material matrix and fill the histograms
  const auto& materialMatrix = bsMaterial.fullMaterial();
  for (auto [b1, materialVector] : enumerate(materialMatrix)) {
    for (auto [b0, mat] : enumerate(materialVector)) {
      t.SetBinContent(static_cast<int>(b0) + 1, static_cast<int>(b1) + 1,
                      mat.thickness());
      x0.SetBinContent(static_cast<int>(b0) + 1, static_cast<int>(b1) + 1,
                       mat.material().X0());
      l0.SetBinContent(static_cast<int>(b0) + 1, static_cast<int>(b1) + 1,
                       mat.material().L0());
      A.SetBinContent(static_cast<int>(b0) + 1, static_cast<int>(b1) + 1,
                      mat.material().Ar());
      Z.SetBinContent(static_cast<int>(b0) + 1, static_cast<int>(b1) + 1,
                      mat.material().Z());
      rho.SetBinContent(static_cast<int>(b0) + 1, static_cast<int>(b1) + 1,
                        mat.material().massDensity());
    }
  }
  t.Write();
  x0.Write();
  l0.Write();
  A.Write();
  Z.Write();
  rho.Write();
}

void RootMaterialMapIo::fillBinnedSurfaceMaterial(
    MaterialTreePayload& payload, const BinnedSurfaceMaterial& bsMaterial) {
  std::size_t bins0 = bsMaterial.binUtility().bins(0);
  std::size_t bins1 = bsMaterial.binUtility().bins(1);

  TH2I idx(m_cfg.indexHistName.c_str(), "indices; bin0; bin1",
           static_cast<int>(bins0), -0.5, static_cast<float>(bins0) - 0.5,
           static_cast<int>(bins1), -0.5, static_cast<float>(bins1) - 0.5);
  // lLop over the material matrix, record the index and fill the indexed tree
  const auto& materialMatrix = bsMaterial.fullMaterial();
  for (auto [b1, materialVector] : enumerate(materialMatrix)) {
    for (auto [b0, mat] : enumerate(materialVector)) {
      idx.SetBinContent(static_cast<int>(b0) + 1, static_cast<int>(b1) + 1,
                        static_cast<float>(payload.index));
      payload.index++;
      fillMaterialSlab(payload, mat);
      m_gTree->Fill();
    }
  }
  idx.Write();
}

void RootMaterialMapIo::writeGridAxisHistograms(
    const GridSurfaceMaterial& gridMaterial) {
  const IMultiAxis2D& multiAxis = gridMaterial.multiAxis();

  // The bin number/direction/boundary-type/range information, one bin per
  // axis (always 2)
  TH1F n(m_cfg.nBinsHistName.c_str(), "bins; bin", 2, -0.5, 1.5);
  TH1F v(m_cfg.axisDirHistName.c_str(), "binning values; bin", 2, -0.5, 1.5);
  TH1F o(m_cfg.axisBoundaryTypeHistName.c_str(), "binning options; bin", 2,
         -0.5, 1.5);
  TH1F rmin(m_cfg.minRangeHistName.c_str(), "min; bin", 2, -0.5, 1.5);
  TH1F rmax(m_cfg.maxRangeHistName.c_str(), "max; bin", 2, -0.5, 1.5);

  for (std::size_t ia = 0; ia < 2; ++ia) {
    const IAxis& axis = multiAxis.getAxis(ia);
    auto ib = static_cast<int>(ia) + 1;
    n.SetBinContent(ib, static_cast<int>(axis.getNBins()));
    // -1 marks "no direction known", since AxisDirection has no such value
    v.SetBinContent(ib, axis.getDirection().has_value()
                            ? static_cast<int>(axis.getDirection().value())
                            : -1);
    o.SetBinContent(ib, static_cast<int>(axis.getBoundaryType()));
    rmin.SetBinContent(ib, axis.getMin());
    rmax.SetBinContent(ib, axis.getMax());
  }
  n.Write();
  v.Write();
  o.Write();
  rmin.Write();
  rmax.Write();
}

void RootMaterialMapIo::writeMaterialKind(MaterialKind kind) {
  TH1I k(m_cfg.materialKindHistName.c_str(), "material kind", 1, -0.5, 0.5);
  k.SetBinContent(1, static_cast<int>(kind));
  k.Write();
}

void RootMaterialMapIo::fillGridSurfaceMaterial(
    const GridSurfaceMaterial& gridMaterial) {
  const IMultiAxis2D& multiAxis = gridMaterial.multiAxis();
  IMultiAxis2D::LocalBins nBins = multiAxis.getNBins();
  auto n0 = static_cast<int>(nBins[0]);
  auto n1 = static_cast<int>(nBins[1]);

  TH2I idx(m_cfg.indexHistName.c_str(), "indices; bin0; bin1", n0, -0.5,
           static_cast<float>(n0) - 0.5, n1, -0.5,
           static_cast<float>(n1) - 0.5);

  // Fill the index histogram over the regular bins, translating a local bin
  // (i0, i1) to a global multi-axis bin and asking @p rowAt for the row in
  // the shared material tree that global bin should point to.
  auto fillIndexHist = [&](auto&& rowAt) {
    for (int i0 = 0; i0 < n0; ++i0) {
      for (int i1 = 0; i1 < n1; ++i1) {
        IMultiAxis2D::LocalBins lbin{static_cast<std::size_t>(i0) + 1,
                                     static_cast<std::size_t>(i1) + 1};
        std::size_t bin = multiAxis.getGlobalBinFromLocalBins(lbin);
        idx.SetBinContent(i0 + 1, i1 + 1, static_cast<float>(rowAt(bin)));
      }
    }
  };

  std::visit(
      [&]<typename T>(const T& storage) {
        if constexpr (std::is_same_v<T, GridSurfaceMaterial::Direct>) {
          // No shared identity to dedup against - one fresh row per bin
          fillIndexHist([&](std::size_t bin) {
            auto row = static_cast<std::size_t>(m_gTree->GetEntries());
            fillMaterialSlab(m_indexedMaterialTreePayload, storage.at(bin));
            m_gTree->Fill();
            return row;
          });
        } else if constexpr (std::is_same_v<T, GridSurfaceMaterial::Indexed>) {
          // Indexed::material is owned per-instance, never shared across
          // surfaces - always append fresh
          auto offset = static_cast<std::size_t>(m_gTree->GetEntries());
          for (const auto& slab : storage.material) {
            fillMaterialSlab(m_indexedMaterialTreePayload, slab);
            m_gTree->Fill();
          }
          fillIndexHist([&](std::size_t bin) {
            return offset + storage.indices.at(bin);
          });
        } else {
          // GloballyIndexed::material is the one canonical global vector
          // for this write session: written once, every surface referencing
          // it (by shared_ptr identity) reuses its row offset. A different
          // vector showing up as GloballyIndexed is a caller error - use
          // Indexed storage for a per-surface, non-shared vector instead.
          const void* identity = storage.material.get();
          if (m_globalMaterialIdentity == nullptr) {
            m_globalMaterialIdentity = identity;
            m_globalMaterialOffset =
                static_cast<std::size_t>(m_gTree->GetEntries());
            for (const auto& slab : *storage.material) {
              fillMaterialSlab(m_indexedMaterialTreePayload, slab);
              m_gTree->Fill();
            }
          } else if (identity != m_globalMaterialIdentity) {
            throw std::invalid_argument(
                "RootMaterialMapIo: a second, distinct GloballyIndexed "
                "material vector was encountered in this write session - "
                "only one canonical global material vector is supported per "
                "session; use Indexed storage instead for a per-surface, "
                "non-shared vector.");
          }
          fillIndexHist([&](std::size_t bin) {
            return m_globalMaterialOffset + storage.indices.at(bin);
          });
        }
      },
      gridMaterial.storage());

  idx.Write();
}

TrackingGeometryMaterial RootMaterialMapIo::read(TFile& rFile,
                                                 const Options& options) {
  TrackingGeometryMaterial detectorMaterial;

  auto& [surfaceMaterials, volumeMaterials] = detectorMaterial;

  auto homogeneousMaterialTree = dynamic_cast<TTree*>(
      rFile.Get(options.homogeneousMaterialTreeName.c_str()));

  // Read homogeneous material tree
  if (homogeneousMaterialTree != nullptr) {
    connectForRead(*homogeneousMaterialTree, m_homogenousMaterialTreePayload);
    for (int i = 0; i < homogeneousMaterialTree->GetEntries(); ++i) {
      homogeneousMaterialTree->GetEntry(i);
      GeometryIdentifier geoID(m_homogenousMaterialTreePayload.hGeoId);
      MaterialSlab materialSlab(
          Material::fromMassDensity(m_homogenousMaterialTreePayload.hX0,
                                    m_homogenousMaterialTreePayload.hL0,
                                    m_homogenousMaterialTreePayload.hA,
                                    m_homogenousMaterialTreePayload.hZ,
                                    m_homogenousMaterialTreePayload.hRho),
          m_homogenousMaterialTreePayload.ht);
      auto homogeneousMaterial =
          std::make_shared<HomogeneousSurfaceMaterial>(materialSlab);
      surfaceMaterials.try_emplace(geoID, homogeneousMaterial);
    }
  }

  // Read the shared indexed material tree, if there - once, in full, so it
  // can be handed out as a single shared vector to every GloballyIndexed
  // GridSurfaceMaterial reconstructed below
  auto indexedMaterialTree =
      dynamic_cast<TTree*>(rFile.Get(options.indexedMaterialTreeName.c_str()));
  std::shared_ptr<std::vector<MaterialSlab>> globalMaterial;
  if (indexedMaterialTree != nullptr) {
    connectForRead(*indexedMaterialTree, m_indexedMaterialTreePayload);
    globalMaterial = std::make_shared<std::vector<MaterialSlab>>();
    auto nEntries = indexedMaterialTree->GetEntries();
    globalMaterial->reserve(static_cast<std::size_t>(nEntries));
    for (Long64_t i = 0; i < nEntries; ++i) {
      indexedMaterialTree->GetEntry(i);
      const auto material = Material::fromMassDensity(
          m_indexedMaterialTreePayload.hX0, m_indexedMaterialTreePayload.hL0,
          m_indexedMaterialTreePayload.hA, m_indexedMaterialTreePayload.hZ,
          m_indexedMaterialTreePayload.hRho);
      globalMaterial->emplace_back(material, m_indexedMaterialTreePayload.ht);
    }
  }

  // Get the list of keys from the file
  TList* tlist = rFile.GetListOfKeys();
  auto tIter = tlist->MakeIterator();
  tIter->Reset();

  // Iterate over the keys in the file
  while (auto key = static_cast<TKey*>(tIter->Next())) {
    // Remember the directory
    std::string tdName(key->GetName());

    ACTS_VERBOSE("Processing directory: " << tdName);

    // volume
    std::vector<std::string> splitNames;
    iter_split(splitNames, tdName,
               boost::algorithm::first_finder(m_cfg.volumePrefix));
    // Surface Material
    if (splitNames[0] == options.folderSurfaceNameBase) {
      boost::split(splitNames, splitNames[1], boost::is_any_of("_"));
      GeometryIdentifier::Value volID = std::stoi(splitNames[0]);
      // boundary
      iter_split(splitNames, tdName,
                 boost::algorithm::first_finder(m_cfg.portalPrefix));
      boost::split(splitNames, splitNames[1], boost::is_any_of("_"));
      GeometryIdentifier::Value bouID = std::stoi(splitNames[0]);
      // layer
      iter_split(splitNames, tdName,
                 boost::algorithm::first_finder(m_cfg.layerPrefix));
      boost::split(splitNames, splitNames[1], boost::is_any_of("_"));
      GeometryIdentifier::Value layID = std::stoi(splitNames[0]);
      // approach
      iter_split(splitNames, tdName,
                 boost::algorithm::first_finder(m_cfg.passivePrefix));
      boost::split(splitNames, splitNames[1], boost::is_any_of("_"));
      GeometryIdentifier::Value appID = std::stoi(splitNames[0]);
      // sensitive
      iter_split(splitNames, tdName,
                 boost::algorithm::first_finder(m_cfg.sensitivePrefix));
      GeometryIdentifier::Value senID = std::stoi(splitNames[1]);

      // Reconstruct the geometry ID
      auto geoID = GeometryIdentifier()
                       .withVolume(volID)
                       .withBoundary(bouID)
                       .withLayer(layID)
                       .withApproach(appID)
                       .withSensitive(senID);

      ACTS_VERBOSE("GeometryIdentifier re-constructed as " << geoID);

      auto surfaceMaterial = readSurfaceMaterial(rFile, tdName, globalMaterial,
                                                 options.gridStorageKind);
      surfaceMaterials.try_emplace(geoID, surfaceMaterial);
    }
  }
  return detectorMaterial;
}

std::shared_ptr<const ISurfaceMaterial> RootMaterialMapIo::readSurfaceMaterial(
    TFile& rFile, const std::string& tdName,
    const std::shared_ptr<std::vector<MaterialSlab>>& globalMaterial,
    Options::GridStorageKind legacyStorageKind) {
  // The material kind marker: present only on directories written by the
  // current writer. Its presence also settles the "o" histogram encoding
  // (see boundaryTypeFromLegacy); its absence marks a legacy
  // BinnedSurfaceMaterial directory.
  std::string kindName = tdName + "/" + m_cfg.materialKindHistName;
  auto kindHist = dynamic_cast<TH1I*>(rFile.Get(kindName.c_str()));
  bool newFormat = kindHist != nullptr;

  // An index histogram means this surface's material lives in the shared
  // indexed material tree - true for every MaterialKind of the current
  // writer, and for legacy indexed BinnedSurfaceMaterial
  std::string indexName = tdName + "/" + m_cfg.indexHistName;
  auto ih = dynamic_cast<TH2I*>(rFile.Get(indexName.c_str()));

  if (newFormat &&
      static_cast<MaterialKind>(static_cast<int>(kindHist->GetBinContent(1))) ==
          MaterialKind::Homogeneous) {
    if (ih == nullptr || globalMaterial == nullptr) {
      ACTS_ERROR(
          "Failed to read homogeneous surface material - missing index "
          "histogram or shared indexed material tree for directory: "
          << tdName);
      return nullptr;
    }
    auto row = static_cast<std::size_t>(ih->GetBinContent(1, 1));
    return std::make_shared<const HomogeneousSurfaceMaterial>(
        globalMaterial->at(row));
  }

  // Construct the common names & get the common histograms
  std::string nName = tdName + "/" + m_cfg.nBinsHistName;
  std::string vName = tdName + "/" + m_cfg.axisDirHistName;
  std::string oName = tdName + "/" + m_cfg.axisBoundaryTypeHistName;
  std::string minName = tdName + "/" + m_cfg.minRangeHistName;
  std::string maxName = tdName + "/" + m_cfg.maxRangeHistName;
  // Get the histograms
  auto n = dynamic_cast<TH1F*>(rFile.Get(nName.c_str()));
  auto v = dynamic_cast<TH1F*>(rFile.Get(vName.c_str()));
  auto o = dynamic_cast<TH1F*>(rFile.Get(oName.c_str()));
  auto minh = dynamic_cast<TH1F*>(rFile.Get(minName.c_str()));
  auto maxh = dynamic_cast<TH1F*>(rFile.Get(maxName.c_str()));

  std::vector<const TH1*> hists{n, v, o, minh, maxh};
  if (std::ranges::any_of(hists,
                          [](const auto* hist) { return hist == nullptr; })) {
    ACTS_ERROR(
        "Failed to read all required histograms for grid surface "
        "material from file: "
        << rFile.GetName());
    return nullptr;
  }
  if (n->GetNbinsX() != 2) {
    ACTS_ERROR(
        "RootMaterialMapIo: only 2D grid surface material is "
        "supported, found "
        << n->GetNbinsX() << " axes in " << tdName);
    return nullptr;
  }

  std::vector<AxisSpec> specs;
  specs.reserve(2);
  for (int ib = 1; ib <= 2; ++ib) {
    auto nbins = static_cast<std::size_t>(n->GetBinContent(ib));
    auto rmin = static_cast<double>(minh->GetBinContent(ib));
    auto rmax = static_cast<double>(maxh->GetBinContent(ib));
    auto dirValue = static_cast<int>(v->GetBinContent(ib));
    std::optional<AxisDirection> direction;
    if (dirValue >= 0) {
      direction = static_cast<AxisDirection>(dirValue);
    }
    auto legacyOrNewOption = static_cast<int>(o->GetBinContent(ib));
    AxisBoundaryType boundaryType =
        newFormat ? static_cast<AxisBoundaryType>(legacyOrNewOption)
                  : boundaryTypeFromLegacy(legacyOrNewOption);
    specs.push_back(
        AxisSpec::Equidistant(nbins, rmin, rmax, boundaryType, direction));
  }
  std::array<AxisSpec, 2> axisSpecs{specs[0], specs[1]};
  MultiAxisSpec2D binning(axisSpecs);
  auto multiAxis = binning.buildMultiAxis();
  const IAxis& axis0 = multiAxis->getAxis(0);
  const IAxis& axis1 = multiAxis->getAxis(1);
  auto nbins0 = axis0.getNBins();
  auto nbins1 = axis1.getNBins();

  ACTS_VERBOSE("Reconstructed grid material binning for " << tdName);

  // A directory carrying a MaterialKind always reconstructs that kind,
  // ignoring legacyStorageKind - a legacy BinnedSurfaceMaterial directory
  // (no MaterialKind) falls back to it instead.
  GridStorageKind storageKind = legacyStorageKind;
  if (newFormat) {
    switch (static_cast<MaterialKind>(
        static_cast<int>(kindHist->GetBinContent(1)))) {
      case MaterialKind::Direct:
        storageKind = GridStorageKind::Direct;
        break;
      case MaterialKind::Indexed:
        storageKind = GridStorageKind::Indexed;
        break;
      case MaterialKind::GloballyIndexed:
        storageKind = GridStorageKind::GloballyIndexed;
        break;
      case MaterialKind::Homogeneous:
        // handled above, before the axis histograms were even read
        break;
    }
  }

  if (ih != nullptr) {
    if (globalMaterial == nullptr) {
      ACTS_ERROR(
          "Found an index histogram but no shared indexed material tree in "
          "file: "
          << rFile.GetName());
      return nullptr;
    }
    std::vector<std::vector<std::size_t>> indexPayload(
        nbins0, std::vector<std::size_t>(nbins1));
    for (std::size_t ib0 = 1; ib0 <= nbins0; ++ib0) {
      for (std::size_t ib1 = 1; ib1 <= nbins1; ++ib1) {
        indexPayload[ib0 - 1][ib1 - 1] = static_cast<std::size_t>(
            ih->GetBinContent(static_cast<int>(ib0), static_cast<int>(ib1)));
      }
    }
    return buildFromIndexPayload(axis0, axis1, indexPayload, globalMaterial,
                                 storageKind);
  }

  // Otherwise fall back to the legacy raw-histogram storage
  std::string tName = tdName + "/" + m_cfg.thicknessHistName;
  std::string x0Name = tdName + "/" + m_cfg.x0HistName;
  std::string l0Name = tdName + "/" + m_cfg.l0HistName;
  std::string aName = tdName + "/" + m_cfg.aHistName;
  std::string zName = tdName + "/" + m_cfg.zHistName;
  std::string rhoName = tdName + "/" + m_cfg.rhoHistName;

  auto t = dynamic_cast<TH2F*>(rFile.Get(tName.c_str()));
  auto x0 = dynamic_cast<TH2F*>(rFile.Get(x0Name.c_str()));
  auto l0 = dynamic_cast<TH2F*>(rFile.Get(l0Name.c_str()));
  auto A = dynamic_cast<TH2F*>(rFile.Get(aName.c_str()));
  auto Z = dynamic_cast<TH2F*>(rFile.Get(zName.c_str()));
  auto rho = dynamic_cast<TH2F*>(rFile.Get(rhoName.c_str()));

  std::vector<const TH1*> materialHists{t, x0, l0, A, Z, rho};
  if (std::ranges::any_of(materialHists,
                          [](const auto* hist) { return hist == nullptr; })) {
    ACTS_ERROR(
        "Failed to read grid surface material - neither an index histogram "
        "nor raw material histograms found for directory: "
        << tdName);
    return nullptr;
  }

  std::vector<std::vector<MaterialSlab>> payload(
      nbins0, std::vector<MaterialSlab>(nbins1, MaterialSlab::Nothing()));
  for (std::size_t ib0 = 1; ib0 <= nbins0; ++ib0) {
    for (std::size_t ib1 = 1; ib1 <= nbins1; ++ib1) {
      auto i0 = static_cast<int>(ib0);
      auto i1 = static_cast<int>(ib1);
      auto dt = static_cast<float>(t->GetBinContent(i0, i1));
      if (dt <= 0.) {
        continue;
      }
      auto dx0 = static_cast<float>(x0->GetBinContent(i0, i1));
      auto dl0 = static_cast<float>(l0->GetBinContent(i0, i1));
      auto da = static_cast<float>(A->GetBinContent(i0, i1));
      auto dz = static_cast<float>(Z->GetBinContent(i0, i1));
      auto drho = static_cast<float>(rho->GetBinContent(i0, i1));
      const auto material = Material::fromMassDensity(dx0, dl0, da, dz, drho);
      payload[ib0 - 1][ib1 - 1] = MaterialSlab(material, dt);
    }
  }
  return buildFromMaterialPayload(axis0, axis1, payload, storageKind);
}

}  // namespace ActsPlugins
