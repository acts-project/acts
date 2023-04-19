// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Detector/LayerStructureKDT.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/CylindricalTrackingGeometry.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <memory>
#include <vector>

using namespace Acts;
using namespace Acts::Test;
using namespace Acts::Experimental;

GeometryContext tContext;
CylindricalTrackingGeometry cGeometry = CylindricalTrackingGeometry(tContext);

namespace {
/// Helper mehtod that allows to use the already existing testing
/// infrastructure with the new const-correct detector design
///
std::vector<std::shared_ptr<Acts::Surface>> unpackSurfaces(
    const std::vector<const Acts::Surface*>& surfaces) {
  std::vector<std::shared_ptr<Acts::Surface>> uSurfaces;
  uSurfaces.reserve(surfaces.size());
  for (const auto& s : surfaces) {
    Surface* ncs = const_cast<Surface*>(s);
    uSurfaces.push_back(ncs->getSharedPtr());
  }
  return uSurfaces;
}

std::vector<std::shared_ptr<Acts::Surface>> pixelSurfaces(
    CylindricalTrackingGeometry::DetectorStore& dStore) {
  // The surfaces for the KDTree structure
  std::vector<std::shared_ptr<Acts::Surface>> pixelSurfaces;
  // Fill Discs
  std::vector<Acts::ActsScalar> pixelDiscs = {-800., -700., -600.,
                                              600.,  700.,  800.};
  for (const auto& z : pixelDiscs) {
    auto rSurfaces = cGeometry.surfacesRing(dStore, 6.4, 12.4, 36., 0.125, 0.,
                                            55., z, 2., 22u);
    auto urSurfaces = unpackSurfaces(rSurfaces);
    pixelSurfaces.insert(pixelSurfaces.end(), urSurfaces.begin(),
                         urSurfaces.end());
  }
  // Fill Barrels
  std::vector<Acts::ActsScalar> pixelBarrels = {32., 72., 116., 172.};
  std::vector<std::pair<int, int>> pixelBinning = {
      {16, 14}, {32, 14}, {52, 14}, {78, 14}};
  for (auto [ir, r] : enumerate(pixelBarrels)) {
    auto cSurfaces = cGeometry.surfacesCylinder(dStore, 8.4, 36., 0.15, 0.145,
                                                r, 3., 2., pixelBinning[ir]);

    auto ucSurfaces = unpackSurfaces(cSurfaces);
    pixelSurfaces.insert(pixelSurfaces.end(), ucSurfaces.begin(),
                         ucSurfaces.end());
  }
  return pixelSurfaces;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(Detector)

// Test only the KDT infrastructure
BOOST_AUTO_TEST_CASE(SurfacesKDT_RZ) {
  // Detector store
  CylindricalTrackingGeometry::DetectorStore dStore;
  auto pSurfaces = pixelSurfaces(dStore);
  // Count the number of surfacees
  size_t refNumber = 6u * 22u + 14u * (16u + 32u + 52u + 78u);
  BOOST_CHECK(pSurfaces.size() == refNumber);

  SurfacesKDT<> skdt(tContext, pSurfaces, {binZ, binR});

  // query: Negative disc 3, it should yield 22 surfaces
  Extent end3;
  end3.set(binZ, -820, -780);
  end3.set(binR, 0., 200.);
  auto nd3 = skdt.surfaces(end3);
  BOOST_CHECK(nd3.size() == 22u);

  // query: 2nd Pixel barrel
  Extent ba1;
  ba1.set(binZ, -580, 580);
  ba1.set(binR, 60., 80.);
  auto b1 = skdt.surfaces(ba1);
  refNumber = 32u * 14u;
  BOOST_CHECK(b1.size() == refNumber);
}

// Test the creation of the Ring
BOOST_AUTO_TEST_CASE(LayerStructureKDT_creationRing) {
  // Detector store
  CylindricalTrackingGeometry::DetectorStore dStore;
  auto pSurfaces = pixelSurfaces(dStore);

  using KDT = SurfacesKDT<2u, 100u, PolyhedronReferenceGenerator>;

  // Count the number of surfacees
  PolyhedronReferenceGenerator pGenerator;
  auto sKdt =
      std::make_shared<KDT>(KDT(tContext, pSurfaces, {binZ, binR}, pGenerator));

  // query: Negative disc 3, it should yield 22 surfaces - No support
  LayerStructureKDT<> lEnd3;
  lEnd3.surfacesKDT = sKdt;
  lEnd3.layerExtent.set(binZ, -820, -780);
  lEnd3.layerExtent.set(binR, 0., 200.);
  lEnd3.surfaceBinning = {LayerStructureKDT<>::Binning{
      Acts::BinningData(Acts::closed, Acts::binPhi, 22u, -M_PI, M_PI), 1u}};

  // The surfaces should be filled and the updator ready
  auto [surfaces, updator] = lEnd3.create(tContext);
  BOOST_TEST(surfaces.size() == 22u);
  BOOST_TEST(updator.connected());
}

// Test the creation of the Cylinder
BOOST_AUTO_TEST_CASE(LayerStructureKDT_creationCylinder) {
  // Detector store
  CylindricalTrackingGeometry::DetectorStore dStore;
  auto pSurfaces = pixelSurfaces(dStore);

  using KDT = SurfacesKDT<2u, 100u, PolyhedronReferenceGenerator>;

  // Count the number of surfacees
  PolyhedronReferenceGenerator pGenerator;
  auto sKdt =
      std::make_shared<KDT>(KDT(tContext, pSurfaces, {binZ, binR}, pGenerator));

  // query: Negative disc 3, it should yield 22 surfaces - No support
  LayerStructureKDT<> lCyl1;
  lCyl1.surfacesKDT = sKdt;

  using LayerBinning = LayerStructureKDT<>::Binning;

  lCyl1.layerExtent.set(binZ, -580, 580);
  lCyl1.layerExtent.set(binR, 60., 80.);
  lCyl1.surfaceBinning = {
      LayerBinning{Acts::BinningData(Acts::open, Acts::binZ, 14u, -480, 480),
                   1u},
      LayerBinning{
          Acts::BinningData(Acts::closed, Acts::binPhi, 32u, -M_PI, M_PI), 1u}};
  // We will add an outside support cylinder
  using LayerSupport = LayerStructureKDT<>::Support;
  lCyl1.layerSupports = {LayerSupport{{15., 10., 10., 0., 0.}}};

  // The surfaces should be filled and the updator ready
  auto [surfaces, updator] = lCyl1.create(tContext);
  BOOST_TEST(surfaces.size() == 32u * 14u + 1u);
  BOOST_TEST(updator.connected());
}

BOOST_AUTO_TEST_SUITE_END()
