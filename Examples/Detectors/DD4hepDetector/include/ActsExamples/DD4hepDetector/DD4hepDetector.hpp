// This file is part of the Acts project.
//
// Copyright (C) 2018-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorStructure.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepGeometryService.hpp"

#include <memory>
#include <tuple>
#include <utility>
#include <vector>

namespace dd4hep {
class Detector;
}  // namespace dd4hep

namespace Acts {
class TrackingGeometry;
class IMaterialDecorator;
class DD4hepFieldAdapter;
}  // namespace Acts

namespace Acts::Examples {
class Detector;
}  // namespace Acts::Examples

namespace ActsExamples {
class IContextDecorator;
}  // namespace ActsExamples

namespace ActsExamples::DD4hep {

class DD4hepDetector {
 public:
  /// @brief The context decorators
  using ContextDecorators =
      std::vector<std::shared_ptr<ActsExamples::IContextDecorator>>;

  /// @brief The tracking geometry
  using TrackingGeometryPtr = std::shared_ptr<const Acts::TrackingGeometry>;

  /// @brief The detector geometry
  using DetectorPtr = std::shared_ptr<const Acts::Experimental::Detector>;

  /// @brief Default constructor
  DD4hepDetector() = default;

  /// @brief Constructor from geometry service
  /// @param geometryService the geometry service
  explicit DD4hepDetector(
      std::unique_ptr<DD4hepGeometryService> geometryService);

  /// @brief Build the tracking geometry from the DD4hep geometry
  ///
  /// @param config is the configuration of the geometry service
  /// @param mdecorator is the material decorator provided
  ///
  /// @return a pair of tracking geometry and context decorators
  std::pair<TrackingGeometryPtr, ContextDecorators> finalize(
      DD4hepGeometryService::Config config,
      std::shared_ptr<const Acts::IMaterialDecorator> mdecorator);

  /// @brief Build the detector from the DD4hep geometry
  ///
  /// @param gctx is the geometry context
  /// @param options is the options struct for the building process
  ///
  /// @note the lifetime of the detector store has to exceed that of the
  ///      detector object as the converted surfaces point back to the
  ///      detector elements
  ///
  /// @return a tuple of detector, context decorators, and the element store
  std::tuple<DetectorPtr, ContextDecorators, Acts::DD4hepDetectorElement::Store>
  finalize(
      const Acts::GeometryContext& gctx,
      const Acts::Experimental::DD4hepDetectorStructure::Options& options = {});

  /// @brief Access to the DD4hep field
  /// @return a shared pointer to the DD4hep field
  std::shared_ptr<Acts::DD4hepFieldAdapter> field() const;

  dd4hep::Detector& dd4hepDetector() const;

  void free();

 private:
  std::unique_ptr<DD4hepGeometryService> m_geometryService = nullptr;
};

}  // namespace ActsExamples::DD4hep
