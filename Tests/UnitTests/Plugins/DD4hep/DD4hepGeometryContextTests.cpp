// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/TransformStore.hpp"
#include "Acts/Plugins/DD4hep/DD4hepGeometryContext.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {
/// Mockup class for the DD4hepDetectorElement
class DD4hepDetectorElement {
 public:
  /// Placeholder for a surface
  std::shared_ptr<Acts::Surface> m_surface = nullptr;

  const Acts::Surface& surface() const {
    if (m_surface) {
      return *m_surface;
    }
    throw std::runtime_error("Surface not set");
  }
};

// Mockup delegate
struct DD4hepAlignmentStore {
  explicit DD4hepAlignmentStore(Acts::TransformStoreGeometryId transformStore)
      : m_transformStore(std::move(transformStore)) {}

  Acts::TransformStoreGeometryId m_transformStore;
  /// Return the contextual transform for a given surface (from detector
  /// element)
  /// @param detElem the dd4hep detector element
  /// @return a Transform3 pointer if found, otherwise nullptr
  const Acts::Transform3* call(const DD4hepDetectorElement& detElem) const {
    // Mockup implementation
    return m_transformStore.contextualTransform(detElem.surface());
  }
};

}  // namespace Acts

BOOST_AUTO_TEST_SUITE(DD4hepPlugin)

BOOST_AUTO_TEST_CASE(DD4hepGeometryContext) {
  auto cylinder0 = Acts::Surface::makeShared<Acts::CylinderSurface>(
      Acts::Transform3::Identity(), 1.0, 2.0);
  cylinder0->assignGeometryId(Acts::GeometryIdentifier().withVolume(1));
  auto contextualTransform0 = Acts::Transform3::Identity();
  contextualTransform0.translation() = Acts::Vector3(1.0, 2.0, 3.0);

  auto cylinder1 = Acts::Surface::makeShared<Acts::CylinderSurface>(
      Acts::Transform3::Identity(), 2.0, 3.0);
  cylinder1->assignGeometryId(Acts::GeometryIdentifier().withVolume(2));

  Acts::TransformStoreGeometryId transformStore(
      {{cylinder0->geometryId(), contextualTransform0}});

  Acts::DD4hepDetectorElement detElem0;
  detElem0.m_surface = cylinder0;

  Acts::DD4hepDetectorElement detElem1;
  detElem1.m_surface = cylinder1;

  // These delegates are created by the framework eventually that reads
  // alignment information from some sort of database and releates it to
  // intervals of validity.
  Acts::DD4hepAlignmentStore alignmentStore(transformStore);

  Acts::DD4hepGeometryContext::Alignment alignmentDelegate;
  alignmentDelegate.connect<&Acts::DD4hepAlignmentStore::call>(&alignmentStore);

  // Create the geometry context
  Acts::DD4hepGeometryContext geometryContext(alignmentDelegate);

  // Check if the contextual transform is available
  const Acts::Transform3* transform0 =
      geometryContext.contextualTransform(detElem0);
  BOOST_CHECK(transform0 != nullptr);
  BOOST_CHECK(transform0->isApprox(contextualTransform0));

  // Check with a different surface
  const Acts::Transform3* transform1 =
      geometryContext.contextualTransform(detElem1);
  BOOST_CHECK(transform1 == nullptr);

  // Test an u-connected delegate
  Acts::DD4hepGeometryContext::Alignment unconnectedDelegate;
  Acts::DD4hepGeometryContext unconnectedContext(unconnectedDelegate);

  BOOST_CHECK(unconnectedContext.contextualTransform(detElem0) == nullptr);
  BOOST_CHECK(unconnectedContext.contextualTransform(detElem1) == nullptr);
}

BOOST_AUTO_TEST_SUITE_END()
