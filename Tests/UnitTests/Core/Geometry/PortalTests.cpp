// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/NavigationDelegates.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

#include <memory>

namespace Acts {
namespace Experimental {
/// Define a dummy detector volume
class DetectorVolume {};

/// a simple link to volume struct
class LinkToVolumeImpl : public IDelegateImpl {
 public:
  std::shared_ptr<DetectorVolume> dVolume = nullptr;

  /// Constructor from volume
  LinkToVolumeImpl(std::shared_ptr<DetectorVolume> dv) : dVolume(dv) {}

  /// @return the link to the contained volume
  /// @note the parameters are ignored
  const DetectorVolume* link(const GeometryContext&, const Vector3&,
                             const Vector3&) const {
    return dVolume.get();
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

BOOST_AUTO_TEST_SUITE(Experimental)

auto volumeA = std::make_shared<DetectorVolume>();
auto volumeB = std::make_shared<DetectorVolume>();
auto dTransform = Acts::Transform3::Identity();

BOOST_AUTO_TEST_CASE(PortalTest) {
  // A rectangle bound surface
  auto rectangle = std::make_shared<Acts::RectangleBounds>(10., 100.);
  auto surface =
      Acts::Surface::makeShared<Acts::PlaneSurface>(dTransform, rectangle);

  // Create a portal out of it
  auto portalA = Portal::makeShared(surface);

  BOOST_TEST(&(portalA->surface()), surface.get());

  portalA->assignGeometryId(Acts::GeometryIdentifier{5});
  BOOST_CHECK(portalA->surface().geometryId() == Acts::GeometryIdentifier{5});

  BOOST_CHECK(portalA == unpackToShared<Portal>(*portalA));
  BOOST_CHECK(portalA == unpackToShared<const Portal>(*portalA));

  // Create a links to volumes
  auto linkToAImpl = std::make_shared<LinkToVolumeImpl>(volumeA);
  auto linkToBImpl = std::make_shared<LinkToVolumeImpl>(volumeB);

  DetectorVolumeLink linkToA;
  linkToA.connect<&LinkToVolumeImpl::link>(linkToAImpl.get());
  ManagedDetectorVolumeLink mLinkToA{std::move(linkToA), linkToAImpl};
  portalA->updateVolumeLink(Acts::NavigationDirection::Forward,
                            std::move(mLinkToA), {volumeA});

  auto attachedVolumes = portalA->attachedVolumes();
  BOOST_CHECK(attachedVolumes[0u].size() == 0u);
  BOOST_CHECK(attachedVolumes[1u].size() == 1u);
  BOOST_CHECK(attachedVolumes[1u][0u] == volumeA);

  // The next volume in forward should be volume A
  auto forwardVolume = portalA->nextVolume(tContext, Acts::Vector3(0., 0., 0.),
                                           Acts::Vector3(0., 0., 1.));
  BOOST_CHECK(forwardVolume == volumeA.get());
  auto backwardVolume = portalA->nextVolume(tContext, Acts::Vector3(0., 0., 0.),
                                            Acts::Vector3(0., 0., -1.));
  // The backwar volume should be nullptr
  BOOST_CHECK(backwardVolume == nullptr);

  auto portalB = Portal::makeShared(surface);
  DetectorVolumeLink linkToB;
  linkToB.connect<&LinkToVolumeImpl::link>(linkToBImpl.get());
  ManagedDetectorVolumeLink mLinkToB{std::move(linkToB), linkToBImpl};
  portalB->updateVolumeLink(Acts::NavigationDirection::Backward,
                            std::move(mLinkToB), {volumeB});

  // Reverse: forwad volume nullptr, backward volume volumeB
  forwardVolume = portalB->nextVolume(tContext, Acts::Vector3(0., 0., 0.),
                                      Acts::Vector3(0., 0., 1.));
  BOOST_CHECK(forwardVolume == nullptr);
  backwardVolume = portalB->nextVolume(tContext, Acts::Vector3(0., 0., 0.),
                                       Acts::Vector3(0., 0., -1.));
  BOOST_CHECK(backwardVolume == volumeB.get());

  // Now fuse the portals together
  portalA->fuse(portalB);
  forwardVolume = portalA->nextVolume(tContext, Acts::Vector3(0., 0., 0.),
                                      Acts::Vector3(0., 0., 1.));
  BOOST_CHECK(forwardVolume == volumeA.get());
  backwardVolume = portalA->nextVolume(tContext, Acts::Vector3(0., 0., 0.),
                                       Acts::Vector3(0., 0., -1.));
  BOOST_CHECK(backwardVolume == volumeB.get());

  // Portal A is now identical to portal B
  BOOST_CHECK(portalA == portalB);

  // An invalid fusing setup
  auto linkToAIImpl = std::make_shared<LinkToVolumeImpl>(volumeA);
  auto linkToBIImpl = std::make_shared<LinkToVolumeImpl>(volumeB);

  auto portalAI = Portal::makeShared(surface);
  DetectorVolumeLink linkToAI;
  linkToAI.connect<&LinkToVolumeImpl::link>(linkToAIImpl.get());
  ManagedDetectorVolumeLink mLinkToAI{std::move(linkToAI), linkToAIImpl};
  portalAI->updateVolumeLink(Acts::NavigationDirection::Forward,
                             std::move(mLinkToAI), {volumeA});

  auto portalBI = Portal::makeShared(surface);
  DetectorVolumeLink linkToBI;
  linkToBI.connect<&LinkToVolumeImpl::link>(linkToBIImpl.get());
  ManagedDetectorVolumeLink mLinkToBI{std::move(linkToBI), linkToBIImpl};
  portalBI->updateVolumeLink(Acts::NavigationDirection::Forward,
                             std::move(mLinkToBI), {volumeB});

  BOOST_CHECK_THROW(portalAI->fuse(portalBI), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
