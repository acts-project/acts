// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Detector/PortalGenerators.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/NavigationDelegates.hpp"
#include "Acts/Navigation/NavigationState.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <memory>

namespace Acts {
namespace Experimental {

/// a simple link to volume struct
class LinkToVolumeImpl : public INavigationDelegate {
 public:
  std::shared_ptr<DetectorVolume> dVolume = nullptr;

  /// Constructor from volume
  LinkToVolumeImpl(std::shared_ptr<DetectorVolume> dv)
      : dVolume(std::move(dv)) {}

  /// @return the link to the contained volume
  /// @note the parameters are ignored
  void link(const GeometryContext& /*gctx*/, NavigationState& nState) const {
    nState.currentVolume = dVolume.get();
  }
};

}  // namespace Experimental
}  // namespace Acts

/// Unpack to shared - simply to test the getSharedPtr mechanism
///
/// @tparam referenced_type is the type of the referenced object
///
/// @param rt is the referenced object
///
/// @returns a shared pointer
template <typename referenced_type>
std::shared_ptr<referenced_type> unpackToShared(referenced_type& rt) {
  return rt.getSharedPtr();
}

using namespace Acts::Experimental;

// A test context
Acts::GeometryContext tContext;

BOOST_AUTO_TEST_SUITE(Detector)

BOOST_AUTO_TEST_CASE(PortalTest) {
  auto dTransform = Acts::Transform3::Identity();
  auto pGenerator = defaultPortalGenerator();
  auto volumeA = DetectorVolumeFactory::construct(
      pGenerator, tContext, "dummyA", dTransform,
      std::make_unique<Acts::CuboidVolumeBounds>(1, 1, 1),
      tryAllPortalsAndSurfaces());
  auto volumeB = DetectorVolumeFactory::construct(
      pGenerator, tContext, "dummyB", dTransform,
      std::make_unique<Acts::CuboidVolumeBounds>(1, 1, 1),
      tryAllPortalsAndSurfaces());

  // A rectangle bound surface
  auto rectangle = std::make_shared<Acts::RectangleBounds>(10., 100.);
  auto surface =
      Acts::Surface::makeShared<Acts::PlaneSurface>(dTransform, rectangle);

  // Create a portal out of it
  auto portalA = Portal::makeShared(surface);

  BOOST_CHECK(&(portalA->surface()) == surface.get());

  portalA->assignGeometryId(Acts::GeometryIdentifier{5});
  BOOST_CHECK(portalA->surface().geometryId() == Acts::GeometryIdentifier{5});

  BOOST_CHECK(portalA == unpackToShared<Portal>(*portalA));
  BOOST_CHECK(portalA == unpackToShared<const Portal>(*portalA));

  // Create a links to volumes
  auto linkToAImpl = std::make_unique<const LinkToVolumeImpl>(volumeA);
  DetectorVolumeUpdator linkToA;
  linkToA.connect<&LinkToVolumeImpl::link>(std::move(linkToAImpl));
  portalA->assignDetectorVolumeUpdator(Acts::Direction::Positive,
                                       std::move(linkToA), {volumeA});

  auto attachedDetectorVolumes = portalA->attachedDetectorVolumes();
  BOOST_CHECK(attachedDetectorVolumes[0u].empty());
  BOOST_CHECK(attachedDetectorVolumes[1u].size() == 1u);
  BOOST_CHECK(attachedDetectorVolumes[1u][0u] == volumeA);

  NavigationState nState;
  nState.position = Acts::Vector3(0., 0., 0.);
  nState.direction = Acts::Vector3(0., 0., 1.);
  // The next volume in positive should be volume A
  portalA->updateDetectorVolume(tContext, nState);
  BOOST_CHECK(nState.currentVolume == volumeA.get());
  // negative should yield nullptr
  nState.direction = Acts::Vector3(0., 0., -1.);
  portalA->updateDetectorVolume(tContext, nState);
  BOOST_CHECK(nState.currentVolume == nullptr);

  auto portalB = Portal::makeShared(surface);
  DetectorVolumeUpdator linkToB;
  auto linkToBImpl = std::make_unique<const LinkToVolumeImpl>(volumeB);
  linkToB.connect<&LinkToVolumeImpl::link>(std::move(linkToBImpl));
  portalB->assignDetectorVolumeUpdator(Acts::Direction::Negative,
                                       std::move(linkToB), {volumeB});

  // Reverse: positive volume nullptr, negative volume volumeB
  nState.direction = Acts::Vector3(0., 0., 1.);
  portalB->updateDetectorVolume(tContext, nState);
  BOOST_CHECK(nState.currentVolume == nullptr);
  nState.direction = Acts::Vector3(0., 0., -1.);
  portalB->updateDetectorVolume(tContext, nState);
  BOOST_CHECK(nState.currentVolume == volumeB.get());

  // Now fuse the portals together, both links valid
  portalA->fuse(portalB);
  nState.direction = Acts::Vector3(0., 0., 1.);
  portalA->updateDetectorVolume(tContext, nState);
  BOOST_CHECK(nState.currentVolume == volumeA.get());
  nState.direction = Acts::Vector3(0., 0., -1.);
  portalA->updateDetectorVolume(tContext, nState);
  BOOST_CHECK(nState.currentVolume == volumeB.get());

  // Portal A is now identical to portal B
  BOOST_CHECK(portalA == portalB);

  // An invalid fusing setup
  auto linkToAIImpl = std::make_unique<const LinkToVolumeImpl>(volumeA);
  auto linkToBIImpl = std::make_unique<const LinkToVolumeImpl>(volumeB);

  auto portalAI = Portal::makeShared(surface);
  DetectorVolumeUpdator linkToAI;
  linkToAI.connect<&LinkToVolumeImpl::link>(std::move(linkToAIImpl));
  portalAI->assignDetectorVolumeUpdator(Acts::Direction::Positive,
                                        std::move(linkToAI), {volumeA});

  auto portalBI = Portal::makeShared(surface);
  DetectorVolumeUpdator linkToBI;
  linkToBI.connect<&LinkToVolumeImpl::link>(std::move(linkToBIImpl));
  portalBI->assignDetectorVolumeUpdator(Acts::Direction::Positive,
                                        std::move(linkToBI), {volumeB});

  BOOST_CHECK_THROW(portalAI->fuse(portalBI), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
