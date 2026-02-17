// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/DiamondVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Result.hpp"

#include <array>
#include <memory>
#include <utility>

namespace {
/// Checks whether the two transforms are the same
/// @param a: First transform to check
/// @param b: Second transform to check
inline bool isSame(const Acts::Transform3& a, const Acts::Transform3& b) {
  const Acts::Transform3 c = a * b.inverse();
  if (c.translation().norm() > Acts::s_onSurfaceTolerance) {
    return false;
  }
  for (std::size_t d = 0; d < 3; ++d) {
    const Acts::Vector3 e = Acts::Vector3::Unit(d);
    if (std::abs(e.dot(c * e) - 1.) > Acts::s_onSurfaceTolerance) {
      return false;
    }
  }
  return true;
}
}  // namespace

namespace Acts {
class PlanarBounds;
}  // namespace Acts

using namespace Acts;
using namespace UnitLiterals;

namespace ActsTests {

/// @class AlignmentContext
struct AlignmentContext {
  /// Return the geometry context
  Acts::GeometryContext getContext() const {
    return Acts::GeometryContext{this};
  }
  /// We have 2 different transforms
  using StoreType_t = std::array<Transform3, 2>;
  std::shared_ptr<const StoreType_t> detElementAlignment = nullptr;
  // Two transforms for the local -> global trf of the volume
  std::shared_ptr<StoreType_t> volumeLocToGlobAlign = nullptr;
  // Two transforms of the global -> local trf of the volume
  std::shared_ptr<StoreType_t> volumeGlobToLocAlign = nullptr;
  // List of portal transforms
  std::vector<StoreType_t> portalAlignments{};

  /// Context index
  unsigned int alignmentIndex{0};

  /// Default constructor
  AlignmentContext() = default;

  /// Constructor with Store and context index
  explicit AlignmentContext(
      std::shared_ptr<const std::array<Transform3, 2>> aStore,
      unsigned int aIndex = 0)
      : detElementAlignment(std::move(aStore)), alignmentIndex(aIndex) {}
};

/// @class AlignableDetectorElement
///
/// This is a lightweight type of detector element,
/// it simply implements the base class.
class AlignableDetectorElement : public SurfacePlacementBase {
 public:
  // Deleted default constructor
  AlignableDetectorElement() = delete;

  /// Constructor for single sided detector element
  /// - bound to a Plane Surface
  ///
  /// @param transform is the transform that element the layer in 3D frame
  /// @param pBounds is the planar bounds for the planar detector element
  /// @param thickness is the module thickness
  AlignableDetectorElement(std::shared_ptr<const Transform3> transform,
                           const std::shared_ptr<const PlanarBounds>& pBounds,
                           double thickness)
      : m_elementTransform(std::move(transform)) {
    m_elementSurface = Surface::makeShared<PlaneSurface>(pBounds, *this);
    m_elementSurface->assignThickness(thickness);
  }

  ///  Destructor
  ~AlignableDetectorElement() override = default;

  /// Return local to global transform associated with this identifier
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  ///
  /// @note this is called from the surface().transform() in the PROXY mode
  const Transform3& localToGlobalTransform(
      const GeometryContext& gctx) const override;

  /// Return surface associated with this detector element
  const Surface& surface() const override;

  /// Non-const access to the surface associated with this detector element
  Surface& surface() override;

  /// Is the detector element a sensitive element
  bool isSensitive() const override { return true; }

 private:
  /// the transform for positioning in 3D space
  std::shared_ptr<const Transform3> m_elementTransform;
  /// the surface represented by it
  std::shared_ptr<Surface> m_elementSurface{nullptr};
};

inline const Transform3& AlignableDetectorElement::localToGlobalTransform(
    const GeometryContext& gctx) const {
  const auto* alignContext = gctx.get<const AlignmentContext*>();
  if (alignContext != nullptr && alignContext->detElementAlignment != nullptr &&
      alignContext->alignmentIndex < 2) {
    return alignContext->detElementAlignment->at(alignContext->alignmentIndex);
  }
  return (*m_elementTransform);
}

inline const Surface& AlignableDetectorElement::surface() const {
  return *m_elementSurface;
}

inline Surface& AlignableDetectorElement::surface() {
  return *m_elementSurface;
}

class AlignableVolumePlacement : public VolumePlacementBase {
 public:
  ///
  explicit AlignableVolumePlacement(const Transform3& volTrf)
      : m_locToGlob{volTrf} {}

  const Transform3& localToGlobalTransform(
      const GeometryContext& gctx) const override {
    const auto* alignContext = gctx.get<const AlignmentContext*>();
    if (alignContext != nullptr &&
        alignContext->volumeLocToGlobAlign != nullptr &&
        alignContext->alignmentIndex < 2) {
      return alignContext->volumeLocToGlobAlign->at(
          alignContext->alignmentIndex);
    }
    return m_locToGlob;
  }

  const Transform3& globalToLocalTransform(
      const GeometryContext& gctx) const override {
    const auto* alignContext = gctx.get<const AlignmentContext*>();
    if (alignContext != nullptr &&
        alignContext->volumeGlobToLocAlign != nullptr &&
        alignContext->alignmentIndex < 2) {
      return alignContext->volumeGlobToLocAlign->at(
          alignContext->alignmentIndex);
    }
    return m_globToLoc;
  }

  const Transform3& portalLocalToGlobal(
      const GeometryContext& gctx, const std::size_t portalIdx) const override {
    const auto* alignContext = gctx.get<const AlignmentContext*>();
    if (alignContext != nullptr && alignContext->alignmentIndex < 2 &&
        alignContext->portalAlignments.size() > portalIdx) {
      return alignContext
          ->portalAlignments[portalIdx][alignContext->alignmentIndex];
    }
    assert(portalIdx < m_portalTrfs.size());
    return m_portalTrfs[portalIdx];
  }

  void setAlignmentDelta(AlignmentContext& context, const Transform3& delta,
                         const std::size_t idx) {
    using StoreType_t = AlignmentContext::StoreType_t;
    assert(idx < 2);

    if (context.volumeLocToGlobAlign == nullptr) {
      context.volumeLocToGlobAlign = std::make_shared<StoreType_t>();
      context.volumeGlobToLocAlign = std::make_shared<StoreType_t>();
    }
    context.volumeLocToGlobAlign->at(idx) = m_locToGlob * delta;
    context.volumeGlobToLocAlign->at(idx) =
        context.volumeLocToGlobAlign->at(idx).inverse();

    context.portalAlignments.resize(nPortalPlacements());
    for (std::size_t portal = 0ul; portal < nPortalPlacements(); ++portal) {
      context.portalAlignments[portal][idx] =
          alignPortal(context.getContext(), portal);
    }
  }

  void makePortalsAlignable(const GeometryContext& gctx,
                            const std::vector<std::shared_ptr<RegularSurface>>&
                                portalsToAlign) override {
    VolumePlacementBase::makePortalsAlignable(gctx, portalsToAlign);
    for (std::size_t portal = 0ul; portal < portalsToAlign.size(); ++portal) {
      m_portalTrfs.push_back(alignPortal(gctx, portal));
    }
  }

 private:
  Transform3 m_locToGlob{Transform3::Identity()};
  Transform3 m_globToLoc{m_locToGlob.inverse()};
  std::vector<Transform3> m_portalTrfs{};
};

}  // namespace ActsTests

using namespace ActsTests;

BOOST_AUTO_TEST_SUITE(GeometrySuite);

/// Unit test for creating compliant/non-compliant Surface object
BOOST_AUTO_TEST_CASE(AlignmentContextTests) {
  // The nominal and alignments
  Vector3 nominalCenter(0., 0., 0.);
  Vector3 negativeCenter(0., 0., -1.);
  Vector3 positiveCenter(0., 0., 1.);

  // Checkpoints
  Vector3 onNominal(3., 3., 0.);
  Vector3 onNegative(3., 3., -1.);
  Vector3 onPositive(3., 3., 1.);

  // Local position
  Vector2 localPosition(3., 3.);

  // A position placeholder and dummy momentum
  Vector3 globalPosition(100_cm, 100_cm, 100_cm);
  Vector3 dummyMomentum(4., 4., 4.);

  Transform3 negativeTransform = Transform3::Identity();
  negativeTransform.translation() = negativeCenter;

  Transform3 positiveTransform = Transform3::Identity();
  positiveTransform.translation() = positiveCenter;

  std::array<Transform3, 2> alignmentArray = {negativeTransform,
                                              positiveTransform};

  auto detElementAlignment =
      std::make_shared<const std::array<Transform3, 2>>(alignmentArray);

  // The detector element at nominal position
  AlignableDetectorElement alignedElement(
      std::make_shared<const Transform3>(Transform3::Identity()),
      std::make_shared<const RectangleBounds>(100_cm, 100_cm), 1_mm);

  const auto& alignedSurface = alignedElement.surface();

  // The alignment contexts
  AlignmentContext alignDefault{};
  AlignmentContext alignNegative{detElementAlignment, 0};
  AlignmentContext alignPositive{detElementAlignment, 1};

  GeometryContext defaultContext{alignDefault.getContext()};
  GeometryContext negativeContext{alignNegative.getContext()};
  GeometryContext positiveContext{alignPositive.getContext()};

  // Test the transforms
  BOOST_CHECK(isSame(alignedSurface.localToGlobalTransform(defaultContext),
                     Transform3::Identity()));
  BOOST_CHECK(isSame(alignedSurface.localToGlobalTransform(negativeContext),
                     negativeTransform));
  BOOST_CHECK(isSame(alignedSurface.localToGlobalTransform(positiveContext),
                     positiveTransform));

  // Test the centers
  BOOST_CHECK_EQUAL(alignedSurface.center(defaultContext), nominalCenter);
  BOOST_CHECK_EQUAL(alignedSurface.center(negativeContext), negativeCenter);
  BOOST_CHECK_EQUAL(alignedSurface.center(positiveContext), positiveCenter);

  // Test OnSurface
  BOOST_CHECK(
      alignedSurface.isOnSurface(defaultContext, onNominal, dummyMomentum));
  BOOST_CHECK(
      alignedSurface.isOnSurface(negativeContext, onNegative, dummyMomentum));
  BOOST_CHECK(
      alignedSurface.isOnSurface(positiveContext, onPositive, dummyMomentum));

  // Test local to Global and vice versa
  globalPosition = alignedSurface.localToGlobal(defaultContext, localPosition,
                                                dummyMomentum);
  BOOST_CHECK_EQUAL(globalPosition, onNominal);
  localPosition =
      alignedSurface.globalToLocal(defaultContext, onNominal, dummyMomentum)
          .value();
  BOOST_CHECK_EQUAL(localPosition, Vector2(3., 3.));

  globalPosition = alignedSurface.localToGlobal(negativeContext, localPosition,
                                                dummyMomentum);
  BOOST_CHECK_EQUAL(globalPosition, onNegative);
  localPosition =
      alignedSurface.globalToLocal(negativeContext, onNegative, dummyMomentum)
          .value();
  BOOST_CHECK_EQUAL(localPosition, Vector2(3., 3.));

  globalPosition = alignedSurface.localToGlobal(positiveContext, localPosition,
                                                dummyMomentum);
  BOOST_CHECK_EQUAL(globalPosition, onPositive);
  localPosition =
      alignedSurface.globalToLocal(positiveContext, onPositive, dummyMomentum)
          .value();
  BOOST_CHECK_EQUAL(localPosition, Vector2(3., 3.));
}

BOOST_AUTO_TEST_CASE(AlignVolumeTests) {
  Transform3 volTrf1{Transform3::Identity()};
  volTrf1.translation() = Vector3{100._mm, 200._mm, 300._mm};
  AlignableVolumePlacement volumePlacement{volTrf1};

  /// define the diamond volume bounds
  constexpr double halfX1 = 10.0_cm;
  constexpr double halfX2 = 12.0_cm;
  constexpr double halfX3 = 12.0_cm;
  constexpr double halfY1 = 5.0_cm;
  constexpr double halfY2 = 10.0_cm;
  constexpr double halfZ = 2.0_cm;

  /// Create an alignable volume
  Volume alignedVol1{volumePlacement,
                     std::make_shared<DiamondVolumeBounds>(
                         halfX1, halfX2, halfX3, halfY1, halfY2, halfZ)};

  BOOST_CHECK_EQUAL(alignedVol1.isAlignable(), true);
  BOOST_CHECK_THROW(alignedVol1.setTransform(volTrf1), std::runtime_error);
  const AlignmentContext defaultContext{};
  // Check that the transform of the alignable volume is defined by the
  // placement

  BOOST_CHECK_EQUAL(
      &alignedVol1.localToGlobalTransform(defaultContext.getContext()),
      &volumePlacement.localToGlobalTransform(defaultContext.getContext()));

  BOOST_CHECK_EQUAL(
      &alignedVol1.globalToLocalTransform(defaultContext.getContext()),
      &volumePlacement.globalToLocalTransform(defaultContext.getContext()));

  BOOST_CHECK(isSame(
      volumePlacement.localToGlobalTransform(defaultContext.getContext()),
      volTrf1));
  // Then fetch the oriented surfaces
  std::vector<OrientedSurface> orientedSurfaces =
      alignedVol1.volumeBounds().orientedSurfaces(
          alignedVol1.localToGlobalTransform(defaultContext.getContext()));

  for (const OrientedSurface& surface : orientedSurfaces) {
    BOOST_CHECK_EQUAL(surface.surface->isAlignable(), false);
    BOOST_CHECK_EQUAL(surface.surface->isSensitive(), false);
  }

  std::vector<std::shared_ptr<RegularSurface>> portalSurfaces{};
  std::ranges::transform(
      orientedSurfaces, std::back_inserter(portalSurfaces),
      [](const OrientedSurface& surface) { return surface.surface; });

  volumePlacement.makePortalsAlignable(defaultContext.getContext(),
                                       portalSurfaces);

  // Reset the surfaces
  orientedSurfaces = alignedVol1.volumeBounds().orientedSurfaces(
      alignedVol1.localToGlobalTransform(defaultContext.getContext()));

  for (std::size_t portal = 0ul; portal < portalSurfaces.size(); ++portal) {
    // The portal mut be alignable
    BOOST_CHECK_EQUAL(portalSurfaces[portal]->isAlignable(), true);
    // But not sensitive
    BOOST_CHECK_EQUAL(portalSurfaces[portal]->isSensitive(), false);
    // Then check that the surface are placed at the same spot
    BOOST_CHECK(isSame(orientedSurfaces[portal].surface->localToGlobalTransform(
                           defaultContext.getContext()),
                       portalSurfaces[portal]->localToGlobalTransform(
                           defaultContext.getContext())));
    // Also check that their bounds are the same
    BOOST_CHECK(orientedSurfaces[portal].surface->bounds() ==
                portalSurfaces[portal]->bounds());
  }

  AlignmentContext alignedContext{};
  Transform3 rotationDelta{};
  rotationDelta.linear() = AngleAxis3{25._degree, Vector3{0., 1., 0.}}.matrix();

  volumePlacement.setAlignmentDelta(alignedContext, rotationDelta, 0);

  BOOST_CHECK(isSame(
      volumePlacement.localToGlobalTransform(alignedContext.getContext()),
      volTrf1 * rotationDelta));

  orientedSurfaces = alignedVol1.volumeBounds().orientedSurfaces(
      alignedVol1.localToGlobalTransform(alignedContext.getContext()));

  for (std::size_t portal = 0ul; portal < portalSurfaces.size(); ++portal) {
    BOOST_CHECK(isSame(orientedSurfaces[portal].surface->localToGlobalTransform(
                           alignedContext.getContext()),
                       portalSurfaces[portal]->localToGlobalTransform(
                           alignedContext.getContext())));
  }

  // Ensure that the bound values can no longer be updated
  BOOST_CHECK_THROW(
      alignedVol1.assignVolumeBounds(std::make_shared<DiamondVolumeBounds>(
          halfX1, 2. * halfX2, halfX3, halfY1, halfY2, halfZ)),
      std::runtime_error);
  // But equivalent bounds can be pushed
  BOOST_CHECK_NO_THROW(
      alignedVol1.assignVolumeBounds(std::make_shared<DiamondVolumeBounds>(
          halfX1, halfX2, halfX3, halfY1, halfY2, halfZ)));
}

BOOST_AUTO_TEST_SUITE_END();
