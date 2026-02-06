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
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Enumerate.hpp"

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

void ActsPlugins::RootMaterialMapIo::write(
    TFile& rFile, const GeometryIdentifier& geoID,
    const ISurfaceMaterial& surfaceMaterial, const Options& options) {
  /// Change to the file
  rFile.cd();

  // Homogeneous surface material writing into one tree
  auto homogeneousMaterial =
      dynamic_cast<const HomogeneousSurfaceMaterial*>(&surfaceMaterial);
  if (homogeneousMaterial != nullptr) {
    if (m_hTree == nullptr) {
      m_hTree = new TTree(options.homogeneousMaterialTreeName.c_str(),
                          "Homogeneous Material Tree");
      connectForWrite(*m_hTree, m_homogenousMaterialTreePayload);
    }
    fillMaterialSlab(m_homogenousMaterialTreePayload,
                     homogeneousMaterial->materialSlab());
    m_homogenousMaterialTreePayload.hGeoId = geoID.value();
    m_hTree->Fill();
    return;
  }

  // Binned surface material writing
  auto bsMaterial =
      dynamic_cast<const BinnedSurfaceMaterial*>(&surfaceMaterial);
  if (bsMaterial != nullptr) {
    // decode the geometryID
    const auto gvolID = geoID.volume();
    const auto gbouID = geoID.boundary();
    const auto glayID = geoID.layer();
    const auto gappID = geoID.approach();
    const auto gsenID = geoID.sensitive();
    // create the directory
    std::string tdName = options.folderSurfaceNameBase.c_str();
    tdName += m_cfg.volumePrefix + std::to_string(gvolID);
    tdName += m_cfg.portalPrefix + std::to_string(gbouID);
    tdName += m_cfg.layerPrefix + std::to_string(glayID);
    tdName += m_cfg.passivePrefix + std::to_string(gappID);
    tdName += m_cfg.sensitivePrefix + std::to_string(gsenID);
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
}

void ActsPlugins::RootMaterialMapIo::write(
    TFile& rFile, const TrackingGeometryMaterial& detectorMaterial,
    const Options& options) {
  const auto& [surfaceMaterials, volumeMaterials] = detectorMaterial;
  for (const auto& [geoID, sMaterial] : surfaceMaterials) {
    write(rFile, geoID, *sMaterial, options);
  }
  if (m_hTree != nullptr) {
    m_hTree->Write();
  }
  if (m_gTree != nullptr) {
    m_gTree->Write();
  }
}

void ActsPlugins::RootMaterialMapIo::connectForWrite(
    TTree& rTree, MaterialTreePayload& treePayload) {
  if (&treePayload == &m_homogenousMaterialTreePayload) {
    rTree.Branch("hGeoId", &treePayload.hGeoId);
  }
  rTree.Branch(m_cfg.thicknessHistName.c_str(), &treePayload.ht);
  rTree.Branch(m_cfg.x0HistName.c_str(), &treePayload.hX0);
  rTree.Branch(m_cfg.l0HistName.c_str(), &treePayload.hL0);
  rTree.Branch(m_cfg.aHistName.c_str(), &treePayload.hA);
  rTree.Branch(m_cfg.zHistName.c_str(), &treePayload.hZ);
  rTree.Branch(m_cfg.rhoHistName.c_str(), &treePayload.hRho);
}

void ActsPlugins::RootMaterialMapIo::connectForRead(
    TTree& rTree, MaterialTreePayload& treePayload) {
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

void ActsPlugins::RootMaterialMapIo::fillMaterialSlab(
    MaterialTreePayload& payload, const MaterialSlab& materialSlab) {
  payload.ht = materialSlab.thickness();
  payload.hX0 = materialSlab.material().X0();
  payload.hL0 = materialSlab.material().L0();
  payload.hA = materialSlab.material().Ar();
  payload.hZ = materialSlab.material().Z();
  payload.hRho = materialSlab.material().massDensity();
}

void ActsPlugins::RootMaterialMapIo::fillBinnedSurfaceMaterial(
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

void ActsPlugins::RootMaterialMapIo::fillBinnedSurfaceMaterial(
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

TrackingGeometryMaterial ActsPlugins::RootMaterialMapIo::read(
    TFile& rFile, const Options& options) {
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

  // Read the binned surface material, if there - connect it to the payload
  auto indexedMaterialTree =
      dynamic_cast<TTree*>(rFile.Get(options.indexedMaterialTreeName.c_str()));
  if (indexedMaterialTree != nullptr) {
    connectForRead(*indexedMaterialTree, m_indexedMaterialTreePayload);
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
      // The surface material to be read in for this
      std::shared_ptr<const ISurfaceMaterial> sMaterial = nullptr;

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

      auto texturedSurfaceMaterial =
          readTextureSurfaceMaterial(rFile, tdName, indexedMaterialTree);
      surfaceMaterials.try_emplace(geoID, texturedSurfaceMaterial);
    }
  }
  return detectorMaterial;
}

std::shared_ptr<const ISurfaceMaterial>
ActsPlugins::RootMaterialMapIo::readTextureSurfaceMaterial(
    TFile& rFile, const std::string& tdName, TTree* indexedMaterialTree) {
  std::shared_ptr<const ISurfaceMaterial> texturedSurfaceMaterial = nullptr;

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
        "Failed to read all required histograms for textured surface "
        "material from file: "
        << rFile.GetName());
    return nullptr;
  }

  // Now reconstruct the bin utilities
  BinUtility bUtility;
  for (int ib = 1; ib < n->GetNbinsX() + 1; ++ib) {
    auto nbins = static_cast<std::size_t>(n->GetBinContent(ib));
    auto val = static_cast<AxisDirection>(v->GetBinContent(ib));
    auto opt = static_cast<BinningOption>(o->GetBinContent(ib));
    auto rmin = static_cast<float>(minh->GetBinContent(ib));
    auto rmax = static_cast<float>(maxh->GetBinContent(ib));
    bUtility += BinUtility(nbins, rmin, rmax, opt, val);
  }
  ACTS_VERBOSE("Created " << bUtility);

  /// Draw from histogram only source
  if (indexedMaterialTree == nullptr) {
    // Construct the names for histogram type storage
    std::string tName = tdName + "/" + m_cfg.thicknessHistName;
    std::string x0Name = tdName + "/" + m_cfg.x0HistName;
    std::string l0Name = tdName + "/" + m_cfg.l0HistName;
    std::string aName = tdName + "/" + m_cfg.aHistName;
    std::string zName = tdName + "/" + m_cfg.zHistName;
    std::string rhoName = tdName + "/" + m_cfg.rhoHistName;

    // Get the histograms
    auto t = dynamic_cast<TH2F*>(rFile.Get(tName.c_str()));
    auto x0 = dynamic_cast<TH2F*>(rFile.Get(x0Name.c_str()));
    auto l0 = dynamic_cast<TH2F*>(rFile.Get(l0Name.c_str()));
    auto A = dynamic_cast<TH2F*>(rFile.Get(aName.c_str()));
    auto Z = dynamic_cast<TH2F*>(rFile.Get(zName.c_str()));
    auto rho = dynamic_cast<TH2F*>(rFile.Get(rhoName.c_str()));

    hists = {t, x0, l0, A, Z, rho};

    // Only go on when you have all histograms
    if (std::ranges::all_of(hists,
                            [](const auto* hist) { return hist != nullptr; })) {
      // Get the number of bins
      int nbins0 = t->GetNbinsX();
      int nbins1 = t->GetNbinsY();
      // The material matrix
      MaterialSlabMatrix materialMatrix(
          nbins1, MaterialSlabVector(nbins0, MaterialSlab::Nothing()));
      // Fill the matrix from the histogram content
      for (int ib0 = 1; ib0 <= nbins0; ++ib0) {
        for (int ib1 = 1; ib1 <= nbins1; ++ib1) {
          auto dt = static_cast<float>(t->GetBinContent(ib0, ib1));
          if (dt > 0.) {
            auto dx0 = static_cast<float>(x0->GetBinContent(ib0, ib1));
            auto dl0 = static_cast<float>(l0->GetBinContent(ib0, ib1));
            auto da = static_cast<float>(A->GetBinContent(ib0, ib1));
            auto dz = static_cast<float>(Z->GetBinContent(ib0, ib1));
            auto drho = static_cast<float>(rho->GetBinContent(ib0, ib1));
            // Create material properties
            const auto material =
                Material::fromMassDensity(dx0, dl0, da, dz, drho);
            materialMatrix[ib1 - 1][ib0 - 1] = MaterialSlab(material, dt);
          }
        }
      }  // Construct the binned material with the right bin utility
      texturedSurfaceMaterial = std::make_shared<const BinnedSurfaceMaterial>(
          bUtility, std::move(materialMatrix));
    }
  } else {
    // Construct the names for histogram type storage
    std::string indexName = tdName + "/" + m_cfg.indexHistName;
    // Get the histograms
    auto ih = dynamic_cast<TH2I*>(rFile.Get(indexName.c_str()));
    if (ih != nullptr) {
      // Get the number of bins
      int nbins0 = ih->GetNbinsX();
      int nbins1 = ih->GetNbinsY();
      // The material matrix
      MaterialSlabMatrix materialMatrix(
          nbins1, MaterialSlabVector(nbins0, MaterialSlab::Nothing()));
      // Fill the matrix from the tree entries
      for (int ib0 = 1; ib0 <= nbins0; ++ib0) {
        for (int ib1 = 1; ib1 <= nbins1; ++ib1) {
          auto idx = static_cast<int>(ih->GetBinContent(ib0, ib1));
          indexedMaterialTree->GetEntry(idx);
          const auto material = Material::fromMassDensity(
              m_indexedMaterialTreePayload.hX0,
              m_indexedMaterialTreePayload.hL0, m_indexedMaterialTreePayload.hA,
              m_indexedMaterialTreePayload.hZ,
              m_indexedMaterialTreePayload.hRho);
          materialMatrix[ib1 - 1][ib0 - 1] =
              MaterialSlab(material, m_indexedMaterialTreePayload.ht);
        }
      }  // Construct the binned material with the right bin utility
      texturedSurfaceMaterial = std::make_shared<const BinnedSurfaceMaterial>(
          bUtility, std::move(materialMatrix));
    }
  }

  return texturedSurfaceMaterial;
}
