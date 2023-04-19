// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/LayerStructureKDT.hpp"
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
  using KDT = Acts::Experimental::SurfacesKDT<
      kDIM, 100u, Acts::Experimental::PolyhedronReferenceGenerator>;

 public:
  /// The KDT surface structure to be created
  std::shared_ptr<KDT> surfacesKDT = nullptr;
  /// Query binning of the KDT surface structure
  std::array<Acts::BinningValue, kDIM> binningValues;

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
      binningValues[ibv] = bv;
    }
    // Get the surfaces
    auto surfaces = readSurfacesFromJson(fname);
    surfacesKDT = std::make_shared<KDT>(
        KDT(Acts::GeometryContext(), surfaces, binningValues));
  }

  /// pyton: axis binning value, axis option, axis type, nbins, boundaries,
  /// expansion
  using Binning = std::tuple<std::string, std::string, std::string,
                             unsigned int, std::vector<float>, unsigned int>;

  /// python: axis support
  using Support = std::tuple<std::array<Acts::ActsScalar, 5u>, unsigned int>;

  // Translate the binning - ACTS version
  using LayerBinning =
      typename Acts::Experimental::LayerStructureKDT<kDIM>::Binning;

  // Translate the support - ACTS class version
  using LayerSupport =
      typename Acts::Experimental::LayerStructureKDT<kDIM>::Support;

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
      const std::string& representation, const std::vector<Binning>& binnings,
      const std::vector<Support>& supports,
      Acts::Svg::IndexedSurfacesConverter::Options drawOptions) {
    // Translate into an Extent
    Acts::Extent lExtent;
    for (const auto& [iq, qr] : Acts::enumerate(qRange)) {
      lExtent.set(binningValues[iq], qr[0u], qr[1u]);
    }
    // Translate the representation
    Acts::Surface::SurfaceType lRep = Acts::Surface::SurfaceType::Other;
    if (representation == "cylinder") {
      lRep = Acts::Surface::SurfaceType::Cylinder;
    } else if (representation == "disc" or representation == "disk") {
      lRep = Acts::Surface::SurfaceType::Disc;
    } else if (representation == "plane") {
      lRep = Acts::Surface::SurfaceType::Plane;
    } else {
      throw std::invalid_argument(
          "Inspectors::SurfaceIndexing: representation type not recognized, "
          "please use 'cylinder', 'disc', or 'plane'.");
    }

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

    // The layer supports
    std::vector<LayerSupport> lSupports;
    for (const auto& ls : supports) {
      lSupports.push_back({std::get<0u>(ls), std::nullopt, std::get<1u>(ls)});
    }

    // Call the ACTS based implementations
    createIndexingImpl(name, lExtent, lRep, lBinnings, lSupports, drawOptions);
  }

  /// Inspect the layer - view implementation using ACTS classes
  ///
  /// @param name of the layer structure
  /// @param lExtent the query range for the surfaces as an extent
  /// @param lRepresentation the representation type of this
  /// @param lBinnings the surface binning description
  /// @param lSupports the supports to be added
  /// @param drawOptions the display/draw options
  ///
  void createIndexingImpl(
      const std::string& name, const Acts::Extent& lExtent,
      Acts::Surface::SurfaceType lRepresentation,
      const std::vector<LayerBinning>& lBinnings,
      const std::vector<LayerSupport>& lSupports,
      Acts::Svg::IndexedSurfacesConverter::Options drawOptions) {
    // Create the layer structure
    Acts::Experimental::LayerStructureKDT<kDIM> lStructure;
    lStructure.surfacesKDT = surfacesKDT;
    lStructure.layerExtent = lExtent;
    lStructure.representation = lRepresentation;
    lStructure.surfaceBinning = lBinnings;
    lStructure.layerSupports = lSupports;

    // The surfaces should be filled and the updator ready
    auto drawContext = Acts::GeometryContext();
    auto [surfaces, updator] = lStructure.create(drawContext);

    // The displaying
    auto pIndexedStructure = Acts::Svg::IndexedSurfacesConverter::convert(
        drawContext, surfaces, updator, drawOptions);
    auto pIndexedView = Acts::Svg::View::xy(pIndexedStructure, name);
    Acts::Svg::toFile({pIndexedView}, pIndexedView._id + ".svg");
  }
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
               const std::string& representation,
               const std::vector<SurfaceIndexing<2u>::Binning>& binnings,
               const std::vector<SurfaceIndexing<2u>::Support>& supports);
};

}  // namespace ActsExamples