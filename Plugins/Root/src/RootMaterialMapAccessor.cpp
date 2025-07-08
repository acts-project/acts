// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Root/RootMaterialMapAccessor.hpp"

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

void Acts::RootMaterialMapAccessor::write(
    TFile& rFile, const GeometryIdentifier& geoID,
    const ISurfaceMaterial& surfaceMaterial) {
  /// Change to the file
  rFile.cd();

  // Homogeneous surface material writing into one tree
  auto homogeneousMaterial =
      dynamic_cast<const HomogeneousSurfaceMaterial*>(&surfaceMaterial);
  if (homogeneousMaterial != nullptr) {
    if (m_hTree == nullptr) {
      m_hTree = new TTree(m_cfg.homogeneousMaterialTree.c_str(),
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
    std::string tdName = m_cfg.folderSurfaceNameBase.c_str();
    tdName += m_cfg.voltag + std::to_string(gvolID);
    tdName += m_cfg.boutag + std::to_string(gbouID);
    tdName += m_cfg.laytag + std::to_string(glayID);
    tdName += m_cfg.apptag + std::to_string(gappID);
    tdName += m_cfg.sentag + std::to_string(gsenID);
    // create a new directory
    rFile.mkdir(tdName.c_str());
    rFile.cd(tdName.c_str());

    // Boundary condistions
    // Get the binning data
    auto& binningData = bsMaterial->binUtility().binningData();
    // 1-D or 2-D maps
    std::size_t binningBins = binningData.size();

    // The bin number information
    auto n = new TH1F(m_cfg.ntag.c_str(), "bins; bin", binningBins, -0.5,
                      binningBins - 0.5);

    // The binning value information
    auto v = new TH1F(m_cfg.vtag.c_str(), "binning values; bin", binningBins,
                      -0.5, binningBins - 0.5);

    // The binning option information
    auto o = new TH1F(m_cfg.otag.c_str(), "binning options; bin", binningBins,
                      -0.5, binningBins - 0.5);

    // The binning option information - range min
    auto rmin = new TH1F(m_cfg.mintag.c_str(), "min; bin", binningBins, -0.5,
                         binningBins - 0.5);

    // The binning option information - range max
    auto rmax = new TH1F(m_cfg.maxtag.c_str(), "max; bin", binningBins, -0.5,
                         binningBins - 0.5);

    // Now fill the histogram content
    std::size_t b = 1;
    for (auto bData : binningData) {
      // Fill: nbins, value, option, min, max
      n->SetBinContent(b, static_cast<int>(binningData[b - 1].bins()));
      v->SetBinContent(b, static_cast<int>(binningData[b - 1].binvalue));
      o->SetBinContent(b, static_cast<int>(binningData[b - 1].option));
      rmin->SetBinContent(b, binningData[b - 1].min);
      rmax->SetBinContent(b, binningData[b - 1].max);
      ++b;
    }
    n->Write();
    v->Write();
    o->Write();
    rmin->Write();
    rmax->Write();

    // If compressed writing is not enabled, write the binned surface material
    // as histograms
    if (!m_cfg.indexedMaterial) {
      fillBinnedSurfaceMaterial(*bsMaterial);
      return;
    }

    // Otherwise, write the binned surface material into the TTree
    if (m_gTree == nullptr) {
      // Back to file level
      rFile.cd();
      m_gTree = new TTree(m_cfg.indexMaterialTreeName.c_str(),
                          "Indexed Material Tree");
      connectForWrite(*m_gTree, m_indexedMaterialTreePayload);
      // Back to the directory
      rFile.cd(tdName.c_str());
    }
    fillBinnedSurfaceMaterial(m_indexedMaterialTreePayload, *bsMaterial);
    return;
  }

  // auto gridMaterial = dynamic_cast<const
  // IGridSurfaceMaterial*>(surfaceMaterial); if (gridMaterial != nullptr) {
  //  writeGridMaterial(rFile, gctx, geoID, *gridMaterial);
  // return;
  //}

  return;
}

void Acts::RootMaterialMapAccessor::write(
    TFile& rFile, const DetectorMaterialMaps& DetectorMaterialMaps) {
  auto [surfaceMaterials, volumeMaterials] = DetectorMaterialMaps;
  for (const auto& [geoID, sMaterial] : surfaceMaterials) {
    write(rFile, geoID, *sMaterial);
  }
}

void Acts::RootMaterialMapAccessor::connectForWrite(
    TTree& rTree, MaterialTreePayload& treePayload) {
  if (&treePayload == &m_homogenousMaterialTreePayload) {
    rTree.Branch("hGeoId", &treePayload.hGeoId);
  }
  rTree.Branch(m_cfg.ttag.c_str(), &treePayload.ht);
  rTree.Branch(m_cfg.x0tag.c_str(), &treePayload.hX0);
  rTree.Branch(m_cfg.l0tag.c_str(), &treePayload.hL0);
  rTree.Branch(m_cfg.atag.c_str(), &treePayload.hA);
  rTree.Branch(m_cfg.ztag.c_str(), &treePayload.hZ);
  rTree.Branch(m_cfg.rhotag.c_str(), &treePayload.hRho);
}

void Acts::RootMaterialMapAccessor::connectForRead(
    TTree& rTree, MaterialTreePayload& treePayload) {
  if (&treePayload == &m_homogenousMaterialTreePayload) {
    rTree.SetBranchAddress("hGeoId", &treePayload.hGeoId);
  }
  rTree.SetBranchAddress(m_cfg.ttag.c_str(), &treePayload.ht);
  rTree.SetBranchAddress(m_cfg.x0tag.c_str(), &treePayload.hX0);
  rTree.SetBranchAddress(m_cfg.l0tag.c_str(), &treePayload.hL0);
  rTree.SetBranchAddress(m_cfg.atag.c_str(), &treePayload.hA);
  rTree.SetBranchAddress(m_cfg.ztag.c_str(), &treePayload.hZ);
  rTree.SetBranchAddress(m_cfg.rhotag.c_str(), &treePayload.hRho);
}

void Acts::RootMaterialMapAccessor::fillMaterialSlab(
    MaterialTreePayload& payload, const MaterialSlab& materialSlab) {
  payload.ht = materialSlab.thickness();
  payload.hX0 = materialSlab.material().X0();
  payload.hL0 = materialSlab.material().L0();
  payload.hA = materialSlab.material().Ar();
  payload.hZ = materialSlab.material().Z();
  payload.hRho = materialSlab.material().massDensity();
}

void Acts::RootMaterialMapAccessor::fillBinnedSurfaceMaterial(
    const BinnedSurfaceMaterial& bsMaterial) {
  std::size_t bins0 = bsMaterial.binUtility().bins(0);
  std::size_t bins1 = bsMaterial.binUtility().bins(1);

  TH2F t(m_cfg.ttag.c_str(), "thickness [mm] ;b0 ;b1", bins0, -0.5, bins0 - 0.5,
         bins1, -0.5, bins1 - 0.5);
  TH2F x0(m_cfg.x0tag.c_str(), "X_{0} [mm] ;b0 ;b1", bins0, -0.5, bins0 - 0.5,
          bins1, -0.5, bins1 - 0.5);
  TH2F l0(m_cfg.l0tag.c_str(), "#Lambda_{0} [mm] ;b0 ;b1", bins0, -0.5,
          bins0 - 0.5, bins1, -0.5, bins1 - 0.5);
  TH2F A(m_cfg.atag.c_str(), "X_{0} [mm] ;b0 ;b1", bins0, -0.5, bins0 - 0.5,
         bins1, -0.5, bins1 - 0.5);
  TH2F Z(m_cfg.ztag.c_str(), "#Lambda_{0} [mm] ;b0 ;b1", bins0, -0.5,
         bins0 - 0.5, bins1, -0.5, bins1 - 0.5);
  TH2F rho(m_cfg.rhotag.c_str(), "#rho [g/mm^3] ;b0 ;b1", bins0, -0.5,
           bins0 - 0.5, bins1, -0.5, bins1 - 0.5);

  // loop over the material and fill
  const auto& materialMatrix = bsMaterial.fullMaterial();
  for (auto [b1, materialVector] : enumerate(materialMatrix)) {
    for (auto [b0, mat] : enumerate(materialVector)) {
      t.SetBinContent(b0 + 1, b1 + 1, mat.thickness());
      x0.SetBinContent(b0 + 1, b1 + 1, mat.material().X0());
      l0.SetBinContent(b0 + 1, b1 + 1, mat.material().L0());
      A.SetBinContent(b0 + 1, b1 + 1, mat.material().Ar());
      Z.SetBinContent(b0 + 1, b1 + 1, mat.material().Z());
      rho.SetBinContent(b0 + 1, b1 + 1, mat.material().massDensity());
    }
  }
  t.Write();
  x0.Write();
  l0.Write();
  A.Write();
  Z.Write();
  rho.Write();
}

void Acts::RootMaterialMapAccessor::fillBinnedSurfaceMaterial(
    MaterialTreePayload& payload, const BinnedSurfaceMaterial& bsMaterial) {
  std::size_t bins0 = bsMaterial.binUtility().bins(0);
  std::size_t bins1 = bsMaterial.binUtility().bins(1);

  TH2I idx(m_cfg.itag.c_str(), "indices; bin0; bin1", bins0, -0.5, bins0 - 0.5,
           bins1, -0.5, bins1 - 0.5);
  // loop over the material and fill
  const auto& materialMatrix = bsMaterial.fullMaterial();
  for (auto [b1, materialVector] : enumerate(materialMatrix)) {
    for (auto [b0, mat] : enumerate(materialVector)) {
      idx.SetBinContent(b0 + 1, b1 + 1, payload.index++);
      fillMaterialSlab(payload, mat);
      m_gTree->Fill();
    }
  }
  idx.Write();
}

Acts::DetectorMaterialMaps Acts::RootMaterialMapAccessor::read(TFile& rFile) {
  Acts::DetectorMaterialMaps DetectorMaterialMaps;

  auto& [surfaceMaterials, volumeMaterials] = DetectorMaterialMaps;

  TTree* homogeneousMaterialTree =
      dynamic_cast<TTree*>(rFile.Get(m_cfg.homogeneousMaterialTree.c_str()));

  // Read homogeneous material tree
  if (homogeneousMaterialTree != nullptr) {
    connectForRead(*homogeneousMaterialTree, m_homogenousMaterialTreePayload);
    for (int i = 0; i < homogeneousMaterialTree->GetEntries(); ++i) {
      homogeneousMaterialTree->GetEntry(i);
      Acts::GeometryIdentifier geoID(m_homogenousMaterialTreePayload.hGeoId);
      Acts::MaterialSlab materialSlab(
          Acts::Material::fromMassDensity(m_homogenousMaterialTreePayload.hX0,
                                          m_homogenousMaterialTreePayload.hL0,
                                          m_homogenousMaterialTreePayload.hA,
                                          m_homogenousMaterialTreePayload.hZ,
                                          m_homogenousMaterialTreePayload.hRho),
          m_homogenousMaterialTreePayload.ht);
      auto homogeneousMaterial =
          std::make_shared<Acts::HomogeneousSurfaceMaterial>(materialSlab);
      surfaceMaterials.emplace(geoID, homogeneousMaterial);
    }
  }

  // Read the binned surface material, if there - connect it to the payload
  TTree* indexedMaterialTree =
      dynamic_cast<TTree*>(rFile.Get(m_cfg.indexMaterialTreeName.c_str()));
  if (indexedMaterialTree != nullptr) {
    connectForRead(*indexedMaterialTree, m_indexedMaterialTreePayload);
  }

  // Get the list of keys from the file
  TList* tlist = rFile.GetListOfKeys();
  auto tIter = tlist->MakeIterator();
  tIter->Reset();

  // Iterate over the keys in the file
  while (TKey* key = static_cast<TKey*>(tIter->Next())) {
    // Remember the directory
    std::string tdName(key->GetName());

    ACTS_VERBOSE("Processing directory: " << tdName);

    // volume
    std::vector<std::string> splitNames;
    iter_split(splitNames, tdName,
               boost::algorithm::first_finder(m_cfg.voltag));
    // Surface Material
    if (splitNames[0] == m_cfg.folderSurfaceNameBase) {
      // The surface material to be read in for this
      std::shared_ptr<const Acts::ISurfaceMaterial> sMaterial = nullptr;

      boost::split(splitNames, splitNames[1], boost::is_any_of("_"));
      Acts::GeometryIdentifier::Value volID = std::stoi(splitNames[0]);
      // boundary
      iter_split(splitNames, tdName,
                 boost::algorithm::first_finder(m_cfg.boutag));
      boost::split(splitNames, splitNames[1], boost::is_any_of("_"));
      Acts::GeometryIdentifier::Value bouID = std::stoi(splitNames[0]);
      // layer
      iter_split(splitNames, tdName,
                 boost::algorithm::first_finder(m_cfg.laytag));
      boost::split(splitNames, splitNames[1], boost::is_any_of("_"));
      Acts::GeometryIdentifier::Value layID = std::stoi(splitNames[0]);
      // approach
      iter_split(splitNames, tdName,
                 boost::algorithm::first_finder(m_cfg.apptag));
      boost::split(splitNames, splitNames[1], boost::is_any_of("_"));
      Acts::GeometryIdentifier::Value appID = std::stoi(splitNames[0]);
      // sensitive
      iter_split(splitNames, tdName,
                 boost::algorithm::first_finder(m_cfg.sentag));
      Acts::GeometryIdentifier::Value senID = std::stoi(splitNames[1]);

      // Reconstruct the geometry ID
      auto geoID = Acts::GeometryIdentifier()
                       .withVolume(volID)
                       .withBoundary(bouID)
                       .withLayer(layID)
                       .withApproach(appID)
                       .withSensitive(senID);

      auto texturedSurfaceMaterial =
          readTextureSurfaceMaterial(rFile, tdName, indexedMaterialTree);
      surfaceMaterials.emplace(geoID, texturedSurfaceMaterial);
    }
  }
  return DetectorMaterialMaps;
}

std::shared_ptr<const Acts::ISurfaceMaterial>
Acts::RootMaterialMapAccessor::readTextureSurfaceMaterial(
    TFile& rFile, const std::string& tdName, TTree* indexedMaterialTree) {
  std::shared_ptr<const Acts::ISurfaceMaterial> texturedSurfaceMaterial =
      nullptr;

  // Construct the common names & get the common histograms
  std::string nName = tdName + "/" + m_cfg.ntag;
  std::string vName = tdName + "/" + m_cfg.vtag;
  std::string oName = tdName + "/" + m_cfg.otag;
  std::string minName = tdName + "/" + m_cfg.mintag;
  std::string maxName = tdName + "/" + m_cfg.maxtag;
  // Get the histograms
  TH1F* n = dynamic_cast<TH1F*>(rFile.Get(nName.c_str()));
  TH1F* v = dynamic_cast<TH1F*>(rFile.Get(vName.c_str()));
  TH1F* o = dynamic_cast<TH1F*>(rFile.Get(oName.c_str()));
  TH1F* min = dynamic_cast<TH1F*>(rFile.Get(minName.c_str()));
  TH1F* max = dynamic_cast<TH1F*>(rFile.Get(maxName.c_str()));

  // Now reconstruct the bin untilities
  BinUtility bUtility;
  for (int ib = 1; ib < n->GetNbinsX() + 1; ++ib) {
    std::size_t nbins = static_cast<std::size_t>(n->GetBinContent(ib));
    auto val = static_cast<AxisDirection>(v->GetBinContent(ib));
    auto opt = static_cast<BinningOption>(o->GetBinContent(ib));
    float rmin = min->GetBinContent(ib);
    float rmax = max->GetBinContent(ib);
    bUtility += Acts::BinUtility(nbins, rmin, rmax, opt, val);
  }
  ACTS_VERBOSE("Created " << bUtility);

  /// Draw from histogram only source
  if (indexedMaterialTree == nullptr) {
    // Construct the names for histogram type storage
    std::string tName = tdName + "/" + m_cfg.ttag;
    std::string x0Name = tdName + "/" + m_cfg.x0tag;
    std::string l0Name = tdName + "/" + m_cfg.l0tag;
    std::string aName = tdName + "/" + m_cfg.atag;
    std::string zName = tdName + "/" + m_cfg.ztag;
    std::string rhoName = tdName + "/" + m_cfg.rhotag;

    // Get the histograms
    TH2F* t = dynamic_cast<TH2F*>(rFile.Get(tName.c_str()));
    TH2F* x0 = dynamic_cast<TH2F*>(rFile.Get(x0Name.c_str()));
    TH2F* l0 = dynamic_cast<TH2F*>(rFile.Get(l0Name.c_str()));
    TH2F* A = dynamic_cast<TH2F*>(rFile.Get(aName.c_str()));
    TH2F* Z = dynamic_cast<TH2F*>(rFile.Get(zName.c_str()));
    TH2F* rho = dynamic_cast<TH2F*>(rFile.Get(rhoName.c_str()));

    std::vector<const TH1*> hists{n, v, o, min, max, t, x0, l0, A, Z, rho};

    // Only go on when you have all histograms
    if (std::ranges::all_of(hists,
                            [](const auto* hist) { return hist != nullptr; })) {
      // Get the number of bins
      int nbins0 = t->GetNbinsX();
      int nbins1 = t->GetNbinsY();

      // The material matrix
      MaterialSlabMatrix materialMatrix(
          nbins1, MaterialSlabVector(nbins0, MaterialSlab::Nothing()));

      // Fill the matrix first
      for (int ib0 = 1; ib0 <= nbins0; ++ib0) {
        for (int ib1 = 1; ib1 <= nbins1; ++ib1) {
          double dt = t->GetBinContent(ib0, ib1);
          if (dt > 0.) {
            double dx0 = x0->GetBinContent(ib0, ib1);
            double dl0 = l0->GetBinContent(ib0, ib1);
            double da = A->GetBinContent(ib0, ib1);
            double dz = Z->GetBinContent(ib0, ib1);
            double drho = rho->GetBinContent(ib0, ib1);
            // Create material properties
            const auto material =
                Acts::Material::fromMassDensity(dx0, dl0, da, dz, drho);
            materialMatrix[ib1 - 1][ib0 - 1] = Acts::MaterialSlab(material, dt);
          }
        }
      }  // Construct the binned material with the right bin utility
      texturedSurfaceMaterial = std::make_shared<const BinnedSurfaceMaterial>(
          bUtility, std::move(materialMatrix));
    }
  } else {
    // Construct the names for histogram type storage
    std::string indexName = tdName + "/" + m_cfg.itag;
    // Get the histograms
    TH2I* ih = dynamic_cast<TH2I*>(rFile.Get(indexName.c_str()));

    if (ih != nullptr) {
      // Get the number of bins
      int nbins0 = ih->GetNbinsX();
      int nbins1 = ih->GetNbinsY();

      // The material matrix
      MaterialSlabMatrix materialMatrix(
          nbins1, MaterialSlabVector(nbins0, MaterialSlab::Nothing()));

      // Fill the matrix first
      for (int ib0 = 1; ib0 <= nbins0; ++ib0) {
        for (int ib1 = 1; ib1 <= nbins1; ++ib1) {
          int idx = static_cast<int>(ih->GetBinContent(ib0, ib1));
          indexedMaterialTree->GetEntry(idx);
          double dt = m_indexedMaterialTreePayload.ht;
          double dx0 = m_indexedMaterialTreePayload.hX0;
          double dl0 = m_indexedMaterialTreePayload.hL0;
          double da = m_indexedMaterialTreePayload.hA;
          double dz = m_indexedMaterialTreePayload.hZ;
          double drho = m_indexedMaterialTreePayload.hRho;
          // Create material properties
          const auto material =
              Acts::Material::fromMassDensity(dx0, dl0, da, dz, drho);
          materialMatrix[ib1 - 1][ib0 - 1] = Acts::MaterialSlab(material, dt);
        }
      }  // Construct the binned material with the right bin utility
      texturedSurfaceMaterial = std::make_shared<const BinnedSurfaceMaterial>(
          bUtility, std::move(materialMatrix));
    }
  }

  return texturedSurfaceMaterial;
}
