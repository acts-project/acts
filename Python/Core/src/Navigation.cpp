// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>
#include <stdexcept>
#include <utility>

#include <boost/core/demangle.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;

namespace ActsPython {

namespace Test {
class DetectorElementStub : public SurfacePlacementBase {
 public:
  DetectorElementStub() = default;
  const Transform3& localToGlobalTransform(
      const GeometryContext& /*gctx*/) const override {
    return m_transform;
  }

  /// Is the detector element a sensitive element
  bool isSensitive() const override { return true; }
  /// Return surface representation - const return pattern
  const Surface& surface() const override {
    throw std::runtime_error("Not implemented");
  }

  /// Non-const return pattern
  Surface& surface() override { throw std::runtime_error("Not implemented"); }

 private:
  Transform3 m_transform{Tranform3::Identity()};
};

}  // namespace Test

/// @brief Add the navigation bindings to a module.
/// @param m the module to add the bindings to
void addNavigation(py::module_& m) {
  {
    auto tryAll =
        py::class_<TryAllNavigationPolicy>(m, "TryAllNavigationPolicy");
    using Config = TryAllNavigationPolicy::Config;
    auto c = py::class_<Config>(tryAll, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, portals, sensitives);
  }

  py::class_<NavigationPolicyFactory, std::shared_ptr<NavigationPolicyFactory>>(
      m, "NavigationPolicyFactory")
      // only to mirror the C++ API
      .def_static("make", []() { return NavigationPolicyFactory{}; })
      .def("add",
           [](NavigationPolicyFactory* self, const py::object& cls) {
             auto mod = py::module_::import("acts");
             if (py::object o = mod.attr("TryAllNavigationPolicy"); cls.is(o)) {
               return std::move(*self).template add<TryAllNavigationPolicy>();
             } else {
               throw std::invalid_argument(
                   "Unknown navigation policy class: " +
                   cls.attr("__name__").cast<std::string>());
             }
           })

      .def("add",
           [](NavigationPolicyFactory* self, const py::object& cls,
              const SurfaceArrayNavigationPolicy::Config& config) {
             auto mod = py::module_::import("acts");
             if (py::object o = mod.attr("SurfaceArrayNavigationPolicy");
                 !cls.is(o)) {
               throw std::invalid_argument(
                   "Unknown navigation policy class: " +
                   cls.attr("__name__").cast<std::string>());
             }

             return std::move(*self).template add<SurfaceArrayNavigationPolicy>(
                 config);
           })

      .def("add",
           [](NavigationPolicyFactory* self, const py::object& cls,
              const TryAllNavigationPolicy::Config& config) {
             auto mod = py::module_::import("acts");
             if (py::object o = mod.attr("TryAllNavigationPolicy");
                 !cls.is(o)) {
               throw std::invalid_argument(
                   "Unknown navigation policy class: " +
                   cls.attr("__name__").cast<std::string>());
             }

             return std::move(*self).template add<TryAllNavigationPolicy>(
                 config);
           })

      .def("_buildTest", [](NavigationPolicyFactory* self) {
        auto vol1 = std::make_shared<TrackingVolume>(
            Transform3::Identity(),
            std::make_shared<CylinderVolumeBounds>(30, 40, 100));
        vol1->setVolumeName("TestVolume");

        auto detElem = std::make_unique<Test::DetectorElementStub>();

        auto surface = Surface::makeShared<CylinderSurface>(
            Transform3::Identity(), std::make_shared<CylinderBounds>(30, 40));
        surface->assignSurfacePlacement(*detElem);

        vol1->addSurface(std::move(surface));

        std::unique_ptr<INavigationPolicy> result =
            self->build(GeometryContext::dangerouslyDefaultConstruct(), *vol1,
                        *getDefaultLogger("Test", Logging::VERBOSE));
      });

  {
    auto saPolicy = py::class_<SurfaceArrayNavigationPolicy>(
        m, "SurfaceArrayNavigationPolicy");

    using LayerType = SurfaceArrayNavigationPolicy::LayerType;
    py::enum_<LayerType>(saPolicy, "LayerType")
        .value("Cylinder", LayerType::Cylinder)
        .value("Disc", LayerType::Disc)
        .value("Plane", LayerType::Plane);

    using Config = SurfaceArrayNavigationPolicy::Config;
    auto c = py::class_<Config>(saPolicy, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, layerType, bins);
  }
}

}  // namespace ActsPython
