// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Detector/detail/KdtSurfacesProvider.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/ActSVG/IndexedSurfacesSvgConverter.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <array>
#include <memory>
#include <optional>
#include <vector>

namespace ActsExamples {

/// @brief  Helper method to read the surfaces from file
///
/// @param fname The file name
///
/// @return a vector of created surfaces for the KDT tree creation
std::vector<std::shared_ptr<Acts::Surface>> readSurfacesFromJson(
    const std::string& fname);

/// The surface indexing struct for KDT driven layer structure creation
template <size_t kDIM>
struct SurfaceIndexing {
  /// The KDT structure definition
  using KDT = Acts::Experimental::detail::KdtSurfaces<
      kDIM, 100u, Acts::Experimental::PolyhedronReferenceGenerator>;

 public:
  /// Constructor from arguments
  ///
  /// @param fname the file name of the json file to be read in
  /// @param bValues the vector of binning values for the KDT construction
  ///
  SurfaceIndexing(const std::string& fname,
                  const std::vector<Acts::BinningValue>& bValues = {
                      Acts::binZ, Acts::binR}) {
    if (bValues.size() != kDIM) {
      throw std::invalid_argument(
          "SurfaceIndexing: wrong number of binning values");
    }
    // Fill the binning Values
    for (auto [ibv, bv] : Acts::enumerate(bValues)) {
      m_binningValues[ibv] = bv;
    }
    // Get the surfaces
    auto surfaces = readSurfacesFromJson(fname);
    m_surfacesKDT = std::make_shared<KDT>(
        KDT(Acts::GeometryContext(), surfaces, m_binningValues));
  }

  /// pyton: axis binning value, axis option, axis type, nbins, boundaries,
  /// expansion
  using Binning = std::tuple<std::string, std::string, std::string,
                             unsigned int, std::vector<float>, unsigned int>;

  /// python: support description, values, type, split
  using Support =
      std::tuple<std::array<Acts::ActsScalar, 5u>, std::string, unsigned int>;

  // Translate the binning - ACTS version
  using LayerBinning =
      typename Acts::Experimental::LayerStructureBuilder::Binning;

  // Translate the support - ACTS class version
  using LayerSupport =
      typename Acts::Experimental::LayerStructureBuilder::Support;

  /// Inspect the layer - written for an easy python binidng
  ///
  /// @param name of the layer structure
  /// @param qRange is the query range for the Surfaces
  /// @param representation of the layer structure
  /// @param binnings the surface binning description
  /// @param supports the supports to be added
  /// @param drawOptions the display/draw options
  ///
  void createIndexing(
      const std::string& name,
      const std::array<std::array<Acts::ActsScalar, 2u>, kDIM> qRange,
      const std::vector<Binning>& binnings,
      const std::vector<Support>& supports,
      Acts::Svg::IndexedSurfacesConverter::Options drawOptions) {
    // Translate into an Extent
    Acts::Extent lExtent;
    for (const auto& [iq, qr] : Acts::enumerate(qRange)) {
      lExtent.set(m_binningValues[iq], qr[0u], qr[1u]);
    }

    // Translate the binnings
    std::vector<LayerBinning> lBinnings;
    for (const auto& bn : binnings) {
      // Translate the number of bins
      unsigned int nBins = std::get<3u>(bn);
      // Translate the boundaries
      const std::vector<float>& edges = std::get<4u>(bn);
      float eMin = edges.front();
      float eMax = edges.back();
      // Translate the value
      Acts::BinningValue bValue = Acts::binPhi;
      std::string bstring = std::get<0u>(bn);
      if (bstring == "binX" or bstring == "x") {
        bValue = Acts::binX;
      } else if (bstring == "binY" or bstring == "y") {
        bValue = Acts::binY;
      } else if (bstring == "binZ" or bstring == "z") {
        bValue = Acts::binZ;
      } else if (bstring == "binR" or bstring == "r") {
        bValue = Acts::binR;
      } else if (bstring == "binPhi" or bstring == "phi") {
        bValue = Acts::binPhi;
      } else {
        throw std::invalid_argument(
            "SurfaceBinning: binning type not supported.");
      }
      // Translate the option
      Acts::BinningOption bOption = Acts::open;
      std::string ostring = std::get<1u>(bn);
      if (ostring == "closed") {
        bOption = Acts::closed;
      }
      // Translate the type
      std::string tstring = std::get<2u>(bn);
      // Expansion
      unsigned int expansion = std::get<5u>(bn);
      if (tstring == "arbitrary" or tstring == "variable") {
        // Set the different type
        lBinnings.push_back(
            {Acts::BinningData(bOption, bValue, edges), expansion});
      } else {
        lBinnings.push_back(
            {Acts::BinningData(bOption, bValue, nBins, eMin, eMax), expansion});
      }
    }

    // Translate the layer supports
    std::vector<LayerSupport> lSupports;
    if (not supports.empty()) {
      std::vector<Acts::BinningValue> sConstraints;
      for (const auto& bv : m_binningValues) {
        sConstraints.push_back(bv);
      }
      // Run over the supports
      for (const auto& ls : supports) {
        std::string type = std::get<1u>(ls);
        // Translate the representation
        Acts::Surface::SurfaceType sType = Acts::Surface::SurfaceType::Other;
        if (type == "cylinder") {
          sType = Acts::Surface::SurfaceType::Cylinder;
        } else if (type == "disc" or type == "disk") {
          sType = Acts::Surface::SurfaceType::Disc;
        } else if (type == "plane") {
          sType = Acts::Surface::SurfaceType::Plane;
        } else {
          throw std::invalid_argument(
              "Inspectors::SurfaceIndexing: representation type not "
              "recognized, "
              "please use 'cylinder', 'disc', or 'plane'.");
        }
        std::array<Acts::ActsScalar, 5u> sValues = std::get<0u>(ls);
        unsigned int sSplits = std::get<2u>(ls);
        lSupports.push_back(
            LayerSupport{sValues, sType, sConstraints, sSplits});
      }
    }
    // Call the ACTS based implementations
    createIndexingImpl(name, lExtent, lBinnings, lSupports, drawOptions);
  }

  /// Inspect the layer - view implementation using ACTS classes
  ///
  /// @param name of the layer structure
  /// @param lExtent the query range for the surfaces as an extent
  /// @param lBinnings the surface binning description
  /// @param lSupports the supports to be added
  /// @param drawOptions the display/draw options
  ///
  void createIndexingImpl(
      const std::string& name, const Acts::Extent& lExtent,
      const std::vector<LayerBinning>& lBinnings,
      const std::vector<LayerSupport>& lSupports,
      Acts::Svg::IndexedSurfacesConverter::Options drawOptions) {
    // A test context for this
    Acts::GeometryContext tContext;

    Acts::Experimental::detail::KdtSurfacesProvider<kDIM> selectedSurfaces;
    selectedSurfaces.kdt = m_surfacesKDT;
    selectedSurfaces.region = lExtent;

    // Configure the layer structure builder
    Acts::Experimental::LayerStructureBuilder::Config lsConfig;
    lsConfig.auxilliary =
        std::string("*** Building ") + name + std::string(" ***");
    lsConfig.surfaces = selectedSurfaces;
    lsConfig.binnings = lBinnings;
    lsConfig.supports = lSupports;

    auto builder = Acts::Experimental::LayerStructureBuilder(
        lsConfig, Acts::getDefaultLogger(name, Acts::Logging::VERBOSE));

    auto [surfaces, volumes, surfacesUpdator, volumeUpdator] =
        builder.create(tContext);

    // The displaying
    auto pIndexedStructure = Acts::Svg::IndexedSurfacesConverter::convert(
        tContext, surfaces, surfacesUpdator, drawOptions);
    auto pIndexedView = Acts::Svg::View::xy(pIndexedStructure, name);
    Acts::Svg::toFile({pIndexedView}, pIndexedView._id + ".svg");
  }

 private:
  /// The KDT surface structure to be created
  std::shared_ptr<KDT> m_surfacesKDT = nullptr;
  /// Query binning of the KDT surface structure
  std::array<Acts::BinningValue, kDIM> m_binningValues;
};

/// Typify for cylinders
class CylindricalDetectorIndexing : public SurfaceIndexing<2u> {
 public:
  /// Constructor with @param fname filename
  CylindricalDetectorIndexing(const std::string& fname);

  /// Inspect the layer - written for an easy python binidng
  ///
  /// @param name of the layer structure
  /// @param qRange is the query range for the Surfaces
  /// @param representation of the layer structure
  /// @param binnings the surface binning description
  /// @param supports the supports to be added
  /// @param drawOptions the display/draw options
  ///
  void inspect(const std::string& name,
               const std::array<std::array<Acts::ActsScalar, 2u>, 2u> qRange,
               const std::vector<SurfaceIndexing<2u>::Binning>& binnings,
               const std::vector<SurfaceIndexing<2u>::Support>& supports);
};

}  // namespace ActsExamples
