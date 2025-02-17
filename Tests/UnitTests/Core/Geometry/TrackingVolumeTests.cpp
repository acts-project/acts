// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/PortalLinkBase.hpp"
#include "Acts/Geometry/TrackingGeometryVisitor.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include "LayerStub.hpp"

using namespace Acts::UnitLiterals;

namespace Acts::Test {

BOOST_AUTO_TEST_SUITE(Geometry)
BOOST_AUTO_TEST_SUITE(TrackingVolumeTests)

std::size_t countVolumes(const TrackingVolume& tv) {
  std::size_t count = 0;
  tv.visitVolumes([&count](const auto&) { ++count; });
  return count;
}

BOOST_AUTO_TEST_CASE(TrackigVolumeChildren) {
  auto cylBounds =
      std::make_shared<Acts::CylinderVolumeBounds>(10_mm, 20_mm, 100_mm);

  TrackingVolume tv{Transform3::Identity(), cylBounds};

  BOOST_CHECK(tv.volumes().empty());
  BOOST_CHECK_EQUAL(countVolumes(tv), 1);

  auto& child1 = tv.addVolume(
      std::make_unique<TrackingVolume>(Transform3::Identity(), cylBounds));

  BOOST_CHECK_EQUAL(tv.volumes().size(), 1);

  auto it = tv.volumes().begin();
  static_assert(std::is_same_v<decltype(*it), TrackingVolume&>);

  const auto& tvConst = tv;
  auto cit = tvConst.volumes().begin();
  static_assert(std::is_same_v<decltype(*cit), const TrackingVolume&>);

  BOOST_CHECK_EQUAL(&*it, &child1);

  BOOST_CHECK_EQUAL(countVolumes(tv), 2);

  tv.addVolume(
      std::make_unique<TrackingVolume>(Transform3::Identity(), cylBounds));

  BOOST_CHECK_EQUAL(countVolumes(tv), 3);
}

namespace {
void testVisitor(auto visitor) {
  BOOST_CHECK(!visitor.m_surfaceCalled);
  BOOST_CHECK(!visitor.m_volumeCalled);
  BOOST_CHECK(!visitor.m_boundaryCalled);
  BOOST_CHECK(!visitor.m_layerCalled);
  BOOST_CHECK(!visitor.m_portalCalled);

  auto surface = Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<RectangleBounds>(10_mm, 10_mm));

  visitor.visitSurface(*surface);
  BOOST_CHECK(visitor.m_surfaceCalled);

  auto cylBounds =
      std::make_shared<Acts::CylinderVolumeBounds>(10_mm, 20_mm, 100_mm);

  TrackingVolume tv{Transform3::Identity(), cylBounds};
  visitor.visitVolume(tv);

  BOOST_CHECK(visitor.m_volumeCalled);

  BoundarySurfaceT<TrackingVolume> boundary{surface, &tv, nullptr};
  visitor.visitBoundarySurface(boundary);
  BOOST_CHECK(visitor.m_boundaryCalled);

  LayerStub layer{nullptr};
  visitor.visitLayer(layer);
  BOOST_CHECK(visitor.m_layerCalled);

  Portal portal{Direction::AlongNormal(), surface, tv};
  visitor.visitPortal(portal);
  BOOST_CHECK(visitor.m_portalCalled);

  visitor.m_volumeCalled = false;
  tv.apply(visitor);
  BOOST_CHECK(visitor.m_volumeCalled);
}
}  // namespace

BOOST_AUTO_TEST_CASE(Visit) {
  struct Visitor : public TrackingGeometryVisitor {
    void visitSurface(const Surface& /*surface*/) override {
      m_surfaceCalled = true;
    }

    void visitVolume(const TrackingVolume& /*volume*/) override {
      m_volumeCalled = true;
    }

    void visitBoundarySurface(
        const BoundarySurfaceT<TrackingVolume>& /*boundary*/) override {
      m_boundaryCalled = true;
    }

    void visitLayer(const Layer& /*layer*/) override { m_layerCalled = true; }

    void visitPortal(const Portal& /*portal*/) override {
      m_portalCalled = true;
    }

    bool m_surfaceCalled = false;
    bool m_volumeCalled = false;
    bool m_boundaryCalled = false;
    bool m_layerCalled = false;
    bool m_portalCalled = false;
  };

  testVisitor(Visitor{});
}

BOOST_AUTO_TEST_CASE(VisitMutable) {
  struct MutableVisitor : public TrackingGeometryMutableVisitor {
    void visitSurface(Surface& /*surface*/) override { m_surfaceCalled = true; }

    void visitVolume(TrackingVolume& /*volume*/) override {
      m_volumeCalled = true;
    }

    void visitBoundarySurface(
        BoundarySurfaceT<TrackingVolume>& /*boundary*/) override {
      m_boundaryCalled = true;
    }

    void visitLayer(Layer& /*layer*/) override { m_layerCalled = true; }

    void visitPortal(Portal& /*portal*/) override { m_portalCalled = true; }

    bool m_surfaceCalled = false;
    bool m_volumeCalled = false;
    bool m_boundaryCalled = false;
    bool m_layerCalled = false;
    bool m_portalCalled = false;
  };

  testVisitor(MutableVisitor{});
}

template <typename Callable>
class LambdaVisitor : public TrackingGeometryVisitor {
 public:
  explicit LambdaVisitor(Callable callable) : m_callable(std::move(callable)) {}

  void visitSurface(const Surface& surface) override {
    if constexpr (std::is_invocable_v<Callable, const Surface&>) {
      m_callable(surface);
    }
  }

  void visitVolume(const TrackingVolume& volume) override {
    if constexpr (std::is_invocable_v<Callable, const TrackingVolume&>) {
      m_callable(volume);
    }
  }

  void visitBoundarySurface(
      const BoundarySurfaceT<TrackingVolume>& boundary) override {
    if constexpr (std::is_invocable_v<
                      Callable, const BoundarySurfaceT<TrackingVolume>&>) {
      m_callable(boundary);
    }
  }

  void visitLayer(const Layer& layer) override {
    if constexpr (std::is_invocable_v<Callable, const Layer&>) {
      m_callable(layer);
    }
  }

  void visitPortal(const Portal& portal) override {
    if constexpr (std::is_invocable_v<Callable, const Portal&>) {
      m_callable(portal);
    }
  }

 private:
  Callable m_callable;
};

template <typename Callable>
class LambdaMutableVisitor : public TrackingGeometryMutableVisitor {
 public:
  explicit LambdaMutableVisitor(Callable callable)
      : m_callable(std::move(callable)) {}

  void visitSurface(Surface& surface) override {
    if constexpr (std::is_invocable_v<Callable, Surface&>) {
      m_callable(surface);
    }
  }

  void visitVolume(TrackingVolume& volume) override {
    if constexpr (std::is_invocable_v<Callable, TrackingVolume&>) {
      m_callable(volume);
    }
  }

  void visitBoundarySurface(
      BoundarySurfaceT<TrackingVolume>& boundary) override {
    if constexpr (std::is_invocable_v<Callable,
                                      BoundarySurfaceT<TrackingVolume>&>) {
      m_callable(boundary);
    }
  }

  void visitLayer(Layer& layer) override {
    if constexpr (std::is_invocable_v<Callable, Layer&>) {
      m_callable(layer);
    }
  }

  void visitPortal(Portal& portal) override {
    if constexpr (std::is_invocable_v<Callable, Portal&>) {
      m_callable(portal);
    }
  }

 private:
  Callable m_callable;
};

BOOST_AUTO_TEST_CASE(VisitLambda) {
  bool surfaceCalled = false;
  bool volumeCalled = false;
  bool boundaryCalled = false;
  bool layerCalled = false;
  bool portalCalled = false;

  auto visitor = LambdaVisitor{overloaded{
      [&](const Surface& /*surface*/) { surfaceCalled = true; },
      [&](const TrackingVolume& /*volume*/) { volumeCalled = true; },
      [&](const BoundarySurfaceT<TrackingVolume>& /*boundary*/) {
        boundaryCalled = true;
      },
      [&](const Layer& /*layer*/) { layerCalled = true; },
      [&](const Portal& /*portal*/) { portalCalled = true; },
  }};

  auto surface = Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<RectangleBounds>(10_mm, 10_mm));

  visitor.visitSurface(*surface);
  BOOST_CHECK(surfaceCalled);

  auto cylBounds =
      std::make_shared<Acts::CylinderVolumeBounds>(10_mm, 20_mm, 100_mm);

  TrackingVolume tv{Transform3::Identity(), cylBounds};
  visitor.visitVolume(tv);

  BOOST_CHECK(volumeCalled);

  BoundarySurfaceT<TrackingVolume> boundary{surface, &tv, nullptr};
  visitor.visitBoundarySurface(boundary);
  BOOST_CHECK(boundaryCalled);

  LayerStub layer{nullptr};
  visitor.visitLayer(layer);
  BOOST_CHECK(layerCalled);

  Portal portal{Direction::AlongNormal(), surface, tv};
  visitor.visitPortal(portal);
  BOOST_CHECK(portalCalled);

  volumeCalled = false;
  tv.apply(visitor);
  BOOST_CHECK(volumeCalled);
}

BOOST_AUTO_TEST_CASE(VisitLambdaMutable) {
  bool surfaceCalled = false;
  bool volumeCalled = false;
  bool boundaryCalled = false;
  bool layerCalled = false;
  bool portalCalled = false;

  auto visitor = LambdaMutableVisitor{overloaded{
      [&](Surface& /*surface*/) { surfaceCalled = true; },
      [&](TrackingVolume& /*volume*/) { volumeCalled = true; },
      [&](BoundarySurfaceT<TrackingVolume>& /*boundary*/) {
        boundaryCalled = true;
      },
      [&](Layer& /*layer*/) { layerCalled = true; },
      [&](Portal& /*portal*/) { portalCalled = true; },
  }};

  auto surface = Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<RectangleBounds>(10_mm, 10_mm));

  visitor.visitSurface(*surface);
  BOOST_CHECK(surfaceCalled);

  auto cylBounds =
      std::make_shared<Acts::CylinderVolumeBounds>(10_mm, 20_mm, 100_mm);

  TrackingVolume tv{Transform3::Identity(), cylBounds};
  visitor.visitVolume(tv);

  BOOST_CHECK(volumeCalled);

  BoundarySurfaceT<TrackingVolume> boundary{surface, &tv, nullptr};
  visitor.visitBoundarySurface(boundary);
  BOOST_CHECK(boundaryCalled);

  LayerStub layer{nullptr};
  visitor.visitLayer(layer);
  BOOST_CHECK(layerCalled);

  Portal portal{Direction::AlongNormal(), surface, tv};
  visitor.visitPortal(portal);
  BOOST_CHECK(portalCalled);

  volumeCalled = false;
  tv.apply(visitor);
  BOOST_CHECK(volumeCalled);
}

template <typename Callable>
consteval bool callableWithMutable() {
  if constexpr (std::is_invocable_v<Callable, Surface&>) {
    return true;
  }
  if constexpr (std::is_invocable_v<Callable, TrackingVolume&>) {
    return true;
  }
  if constexpr (std::is_invocable_v<Callable,
                                    BoundarySurfaceT<TrackingVolume>&>) {
    return true;
  }
  if constexpr (std::is_invocable_v<Callable, Layer&>) {
    return true;
  }
  if constexpr (std::is_invocable_v<Callable, Portal&>) {
    return true;
  }
  return false;
}

template <typename Callable>
consteval bool callableWithConst() {
  if constexpr (std::is_invocable_v<Callable, const Surface&>) {
    return true;
  }
  if constexpr (std::is_invocable_v<Callable, const TrackingVolume&>) {
    return true;
  }
  if constexpr (std::is_invocable_v<Callable,
                                    const BoundarySurfaceT<TrackingVolume>&>) {
    return true;
  }
  if constexpr (std::is_invocable_v<Callable, const Layer&>) {
    return true;
  }
  if constexpr (std::is_invocable_v<Callable, const Portal&>) {
    return true;
  }
  return false;
}

template <typename Callable>
consteval bool callableWithAny() {
  return std::is_invocable_v<Callable, Surface&> ||
         std::is_invocable_v<Callable, const Surface&> ||
         std::is_invocable_v<Callable, TrackingVolume&> ||
         std::is_invocable_v<Callable, const TrackingVolume&> ||
         std::is_invocable_v<Callable, BoundarySurfaceT<TrackingVolume>&> ||
         std::is_invocable_v<Callable,
                             const BoundarySurfaceT<TrackingVolume>&> ||
         std::is_invocable_v<Callable, Layer&> ||
         std::is_invocable_v<Callable, const Layer&> ||
         std::is_invocable_v<Callable, Portal&> ||
         std::is_invocable_v<Callable, const Portal&>;
}

template <typename Callable>
  requires(callableWithAny<Callable>())
class AnyLambda : public std::conditional_t<callableWithConst<Callable>(),
                                            LambdaVisitor<Callable>,
                                            LambdaMutableVisitor<Callable>> {
 public:
  using Base =
      std::conditional_t<callableWithConst<Callable>(), LambdaVisitor<Callable>,
                         LambdaMutableVisitor<Callable>>;

  explicit AnyLambda(Callable callable) : Base(std::move(callable)) {}
};

BOOST_AUTO_TEST_CASE(DeductBaseClass) {
  static_assert(std::is_base_of_v<TrackingGeometryMutableVisitor,
                                  decltype(AnyLambda{overloaded{
                                      [&](Surface& /*surface*/) {},
                                  }})>);

  static_assert(std::is_base_of_v<TrackingGeometryVisitor,
                                  decltype(AnyLambda{overloaded{
                                      [&](const Surface& /*surface*/) {},
                                  }})>);

  // Test with TrackingVolume
  static_assert(std::is_base_of_v<TrackingGeometryMutableVisitor,
                                  decltype(AnyLambda{overloaded{
                                      [&](TrackingVolume& /*volume*/) {},
                                  }})>);

  static_assert(std::is_base_of_v<TrackingGeometryVisitor,
                                  decltype(AnyLambda{overloaded{
                                      [&](const TrackingVolume& /*volume*/) {},
                                  }})>);

  // Test with BoundarySurface
  static_assert(std::is_base_of_v<
                TrackingGeometryMutableVisitor,
                decltype(AnyLambda{overloaded{
                    [&](BoundarySurfaceT<TrackingVolume>& /*boundary*/) {},
                }})>);

  static_assert(
      std::is_base_of_v<
          TrackingGeometryVisitor,
          decltype(AnyLambda{overloaded{
              [&](const BoundarySurfaceT<TrackingVolume>& /*boundary*/) {},
          }})>);

  // Test with Layer
  static_assert(std::is_base_of_v<TrackingGeometryMutableVisitor,
                                  decltype(AnyLambda{overloaded{
                                      [&](Layer& /*layer*/) {},
                                  }})>);

  static_assert(std::is_base_of_v<TrackingGeometryVisitor,
                                  decltype(AnyLambda{overloaded{
                                      [&](const Layer& /*layer*/) {},
                                  }})>);

  // Test with Portal
  static_assert(std::is_base_of_v<TrackingGeometryMutableVisitor,
                                  decltype(AnyLambda{overloaded{
                                      [&](Portal& /*portal*/) {},
                                  }})>);

  static_assert(std::is_base_of_v<TrackingGeometryVisitor,
                                  decltype(AnyLambda{overloaded{
                                      [&](const Portal& /*portal*/) {},
                                  }})>);

  // Test with multiple types
  static_assert(std::is_base_of_v<TrackingGeometryMutableVisitor,
                                  decltype(AnyLambda{overloaded{
                                      [&](Surface& /*surface*/) {},
                                      [&](TrackingVolume& /*volume*/) {},
                                      [&](Layer& /*layer*/) {},
                                  }})>);

  static_assert(std::is_base_of_v<TrackingGeometryVisitor,
                                  decltype(AnyLambda{overloaded{
                                      [&](const Surface& /*surface*/) {},
                                      [&](const TrackingVolume& /*volume*/) {},
                                      [&](const Layer& /*layer*/) {},
                                  }})>);

  // Test callableWithAny
  static_assert(callableWithAny<decltype([](Surface& /*surface*/) {})>());
  static_assert(callableWithAny<decltype([](const Surface& /*surface*/) {})>());
  static_assert(callableWithAny<decltype([](Portal& /*portal*/) {})>());
  static_assert(callableWithAny<decltype([](const Portal& /*portal*/) {})>());
  static_assert(callableWithAny<decltype([](TrackingVolume& /*volume*/) {})>());
  static_assert(!callableWithAny<decltype([](int /*x*/) {})>());
  static_assert(!callableWithAny<decltype([](void* /*ptr*/) {})>());

  // Test mixed const/non-const overloads
  static_assert(callableWithAny<decltype(overloaded{
                    [](Surface& /*surface*/) {},
                    [](const TrackingVolume& /*volume*/) {},
                })>());

  // Test with unrelated overloads
  static_assert(!callableWithAny<decltype(overloaded{
                    [](int /*x*/) {},
                    [](double /*y*/) {},
                    [](const std::string& /*s*/) {},
                })>());

  // Test with mix of geometry and unrelated overloads
  static_assert(callableWithAny<decltype(overloaded{
                    [](int /*x*/) {},
                    [](Surface& /*surface*/) {},
                    [](const std::string& /*s*/) {},
                })>());

  static_assert(callableWithAny<decltype(overloaded{
                    [](void* /*ptr*/) {},
                    [](const Portal& /*portal*/) {},
                    [](double /*d*/) {},
                })>());
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
