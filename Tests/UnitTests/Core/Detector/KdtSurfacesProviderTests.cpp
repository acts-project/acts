// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Detector/detail/KdtSurfacesProvider.hpp"
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
/// Helper method that allows to use the already existing testing
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
BOOST_AUTO_TEST_CASE(KdtSurfacesProvider) {
  // Detector store
  CylindricalTrackingGeometry::DetectorStore dStore;
  auto pSurfaces = pixelSurfaces(dStore);
  // Count the number of surfacees
  size_t refNumber = 6u * 22u + 14u * (16u + 32u + 52u + 78u);
  BOOST_CHECK(pSurfaces.size() == refNumber);

  using KDTS = Acts::Experimental::detail::KdtSurfaces<>;
  auto skdt = std::make_shared<KDTS>(KDTS(tContext, pSurfaces, {binZ, binR}));

  // query: Negative disc 3, it should yield 22 surfaces
  Acts::Experimental::detail::KdtSurfacesProvider<> end3;
  end3.kdt = skdt;
  end3.region.set(binZ, -820, -780);
  end3.region.set(binR, 0., 200.);
  auto nd3 = end3();
  BOOST_CHECK(nd3.size() == 22u);

  // query: 2nd Pixel barrel
  Acts::Experimental::detail::KdtSurfacesProvider<> ba1;
  ba1.kdt = skdt;
  ba1.region.set(binZ, -580, 580);
  ba1.region.set(binR, 60., 80.);
  auto b1 = ba1();
  refNumber = 32u * 14u;
  BOOST_CHECK(b1.size() == refNumber);
}

BOOST_AUTO_TEST_SUITE_END()
