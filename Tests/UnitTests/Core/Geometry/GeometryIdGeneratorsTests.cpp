// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/DetectorVolume.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/NavigationState.hpp"
#include "Acts/Geometry/detail/DetectorVolumeUpdators.hpp"
#include "Acts/Geometry/detail/GeometryIdGenerators.hpp"
#include "Acts/Geometry/detail/PortalGenerators.hpp"
#include "Acts/Geometry/detail/PortalHelper.hpp"
#include "Acts/Geometry/detail/SurfaceCandidatesUpdators.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <array>
#include <memory>
#include <vector>

namespace Acts {

/// A mockup detector element
class TestDetectorElement : public DetectorElementBase {
 public:
  Transform3 m_transform = Transform3::Identity();
  std::shared_ptr<PlaneSurface> m_surface = nullptr;

  /// @brief  Constructor of a dummy test element
  TestDetectorElement(const Transform3& trf = Transform3::Identity())
      : m_transform(trf) {
    m_surface = Surface::makeShared<PlaneSurface>(
        std::make_shared<RectangleBounds>(10., 10.), *this);
  }

  /// Return the transform for the Element proxy mechanism
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  const Transform3& transform(const GeometryContext& /*gctx*/) const final {
    return m_transform;
  };

  /// Return surface representation const access
  const Surface& surface() const final { return *(m_surface.get()); }

  /// Return surface representation, non- const access
  Surface& surface() { return *(m_surface.get()); }

  /// Returns the thickness of the module
  /// @return double that indicates the thickness of the module
  double thickness() const final { return 1.; }
};

}  // namespace Acts

Acts::GeometryContext tContext;

BOOST_AUTO_TEST_SUITE(Experimental)

namespace {

std::vector<Acts::TestDetectorElement> centralStore = {};

/// @brief  This generator a vector of mockup volumes
///
///  | A | B | C | D  | E | F |
///  |   | s |   | sp | p | v |
///
/// A -> empty volume
/// B -> volume with sensitive surface
/// C -> empty volume
/// D -> volume with sensitive and passive surface
/// E -> volume with passive surface
/// F -> volume with sub volume
///
/// @return
std::array<std::shared_ptr<Acts::Experimental::DetectorVolume>, 6u>
mockupVolumes() {
  auto portals = Acts::Experimental::detail::defaultPortalGenerator();

  Acts::ActsScalar hLength = 100.;

  std::array<std::shared_ptr<Acts::Experimental::DetectorVolume>, 6u> mVolumes =
      {nullptr, nullptr, nullptr, nullptr, nullptr, nullptr};
  std::shared_ptr<Acts::Experimental::DetectorVolume> lastVolume = nullptr;

  std::vector<std::string> volumeNames = {"A", "Bs", "C", "Dsp", "Ep", "Fv"};

  for (auto [in, n] : Acts::enumerate(volumeNames)) {
    // Position the volume
    Acts::Transform3 vTransform = Acts::Transform3::Identity();
    vTransform.pretranslate(2 * in * Acts::Vector3(0., 0., hLength));
    // Volume bounds
    auto vBounds =
        std::make_unique<Acts::CuboidVolumeBounds>(hLength, hLength, hLength);
    std::vector<std::shared_ptr<Acts::Surface>> surfaces = {};
    std::vector<std::shared_ptr<Acts::Experimental::DetectorVolume>> volumes =
        {};
    // Create constituents
    // - sensitive
    if (n.find("s") != std::string::npos) {
      centralStore.push_back(Acts::TestDetectorElement(vTransform));
      surfaces.push_back(centralStore.back().surface().getSharedPtr());
    }
    // - passive
    if (n.find("p") != std::string::npos) {
      surfaces.push_back(Acts::Surface::makeShared<Acts::PlaneSurface>(
          vTransform.pretranslate(Acts::Vector3(0., 0., 0.5 * hLength)),
          std::make_shared<Acts::RectangleBounds>(10., 10.)));
    }
    // - sub volumes
    if (n.find("v") != std::string::npos) {
      auto innerVolumeBounds = std::make_unique<Acts::CuboidVolumeBounds>(
          0.5 * hLength, 0.5 * hLength, 0.5 * hLength);

      // Create the inner box
      auto innerVolume = Acts::Experimental::DetectorVolumeFactory::construct(
          portals, tContext, "InnerVolume", vTransform,
          std::move(innerVolumeBounds),
          Acts::Experimental::detail::allPortals());

      volumes.push_back(innerVolume);
    }
    // Create the volume
    auto volumeBounds =
        std::make_unique<Acts::CuboidVolumeBounds>(hLength, hLength, hLength);

    auto volume = Acts::Experimental::DetectorVolumeFactory::construct(
        portals, tContext, n, vTransform, std::move(volumeBounds), surfaces,
        volumes, Acts::Experimental::detail::allPortals());
    // Fuse if necessary
    if (lastVolume != nullptr) {
      lastVolume->portalPtrs()[1u]->fuse(volume->portalPtrs()[0u]);
      volume->updatePortal(volume->portalPtrs()[0u], 0u);
    }

    // To make it realistic, fuse the volumes
    mVolumes[in] = volume;
    lastVolume = volume;
  }
  return mVolumes;
}

template <typename restricter>
void runRestrictedTest(
    restricter rGen,
    const std::array<std::shared_ptr<Acts::Experimental::DetectorVolume>, 6u>&
        vols) {
  auto [A, B, C, D, E, F] = vols;

  // For the rest of the volume, we do a standard recursive generation
  Acts::Experimental::detail::VolumeCounter vCounterAll{true, 1u};
  // Unset checker
  Acts::Experimental::detail::UnsetIdChecker uic;
  // Duplicate checker
  Acts::Experimental::detail::DuplicateIdChecker dic;

  // Final generator
  // - first the sensitive one
  // - then the rest
  // - then uset checker
  // - then duplicate checker
  Acts::Experimental::detail::ChainedGeometryIdGenerator<
      decltype(rGen), decltype(vCounterAll), decltype(uic), decltype(dic)>
      idGenerator(std::tie(rGen, vCounterAll, uic, dic));

  // That how they would be called in sequence
  idGenerator.generateIds(*A);
  idGenerator.generateIds(*B);
  idGenerator.generateIds(*C);
  idGenerator.generateIds(*D);
  idGenerator.generateIds(*E);
  idGenerator.generateIds(*F);

  // A should volumeID = 2
  BOOST_CHECK(A->geometryId().volume() == 2u);
  // B and D shhould have volumeID = 1, but increasing layer IDs
  BOOST_CHECK(B->geometryId().volume() == 1u);
  BOOST_CHECK(B->geometryId().layer() == 1u);
  BOOST_CHECK(D->geometryId().volume() == 1u);
  BOOST_CHECK(D->geometryId().layer() == 2u);
  // C should volumeID = 3
  BOOST_CHECK(C->geometryId().volume() == 3u);
  // E should volumeID = 4
  BOOST_CHECK(E->geometryId().volume() == 4u);
  // F should volumeID = 5
  BOOST_CHECK(F->geometryId().volume() == 5u);

  // The sensitive surface of B - sensitive one
  BOOST_CHECK(B->surfaces()[0u]->geometryId().sensitive() == 1u);
  // The passive surface of D - passive one
  BOOST_CHECK(D->surfaces()[0u]->geometryId().sensitive() == 1u);
  BOOST_CHECK(D->surfaces()[1u]->geometryId().passive() == 1u);
}

}  // namespace

BOOST_AUTO_TEST_CASE(LayeredDetectorTests_VolumeIdentified) {
  auto [A, B, C, D, E, F] = mockupVolumes();
  // Volume B and D have sensitive surfaces and will be counted
  // as layers of volume 1, their portals and their sensitive and passive will
  // just be counted per volume
  Acts::Experimental::detail::LayerCounter lCounter1{
      Acts::GeometryIdentifier().setVolume(1)};
  Acts::Experimental::detail::PortalCounter poCounter1;
  Acts::Experimental::detail::SensitiveCounter sCounter1;
  Acts::Experimental::detail::PassiveCounter pCounter1;
  // Chain these to gether
  Acts::Experimental::detail::ChainedGeometryIdGenerator<
      decltype(lCounter1), decltype(poCounter1), decltype(sCounter1),
      decltype(pCounter1)>
      chGenerator1(std::tie(lCounter1, poCounter1, sCounter1, pCounter1));
  // These are restricted to volumes B, D
  Acts::Experimental::detail::VolumeRestrictedIdGenerator<decltype(
      chGenerator1)>
      vrIdGenerator(chGenerator1, {B.get(), D.get()});
  runRestrictedTest<decltype(vrIdGenerator)>(vrIdGenerator, {A, B, C, D, E, F});
}

BOOST_AUTO_TEST_CASE(LayeredDetectorTests_NameIdentified) {
  auto [A, B, C, D, E, F] = mockupVolumes();
  // Volume B and D have sensitive surfaces and will be counted
  // as layers of volume 1, their portals and their sensitive and passive will
  // just be counted per volume
  Acts::Experimental::detail::LayerCounter lCounter1{
      Acts::GeometryIdentifier().setVolume(1)};
  Acts::Experimental::detail::PortalCounter poCounter1;
  Acts::Experimental::detail::SensitiveCounter sCounter1;
  Acts::Experimental::detail::PassiveCounter pCounter1;
  // Chain these to gether
  Acts::Experimental::detail::ChainedGeometryIdGenerator<
      decltype(lCounter1), decltype(poCounter1), decltype(sCounter1),
      decltype(pCounter1)>
      chGenerator1(std::tie(lCounter1, poCounter1, sCounter1, pCounter1));
  // These are restricted to volumes B, D
  Acts::Experimental::detail::NameRestrictedIdGenerator<decltype(chGenerator1)>
      vrIdGenerator(chGenerator1, "s");
  runRestrictedTest<decltype(vrIdGenerator)>(vrIdGenerator, {A, B, C, D, E, F});
}

BOOST_AUTO_TEST_CASE(DuplicateAndUnsetIdCheckerTest) {
  auto [A, B, C, D, E, F] = mockupVolumes();
  A->assignGeometryId(Acts::GeometryIdentifier().setVolume(3u));
  C->assignGeometryId(Acts::GeometryIdentifier().setVolume(3u));

  // Duplicate checker
  Acts::Experimental::detail::DuplicateIdChecker dic;
  dic.generateIds(*A);
  BOOST_CHECK_THROW(dic.generateIds(*C), std::runtime_error);

  Acts::Experimental::detail::UnsetIdChecker uic;
  BOOST_CHECK_THROW(uic.generateIds(*B), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
