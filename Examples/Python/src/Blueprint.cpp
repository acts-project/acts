// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/BlueprintNode.hpp"
#include "Acts/Geometry/CylinderContainerBlueprintNode.hpp"
#include "Acts/Geometry/CylinderVolumeStack.hpp"
#include "Acts/Geometry/MaterialDesignatorBlueprintNode.hpp"
#include "Acts/Geometry/StaticBlueprintNode.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace Acts::Python {
void addBlueprint(Context& ctx) {
  auto m = ctx.get("main");

  struct AddCylinderContainerHelper {
    BlueprintNode& parent;
    BinningValue direction;
    std::optional<std::string> name;

    std::shared_ptr<CylinderContainerBlueprintNode> operator()(
        py::function callback) {
      auto cylinder = std::make_shared<CylinderContainerBlueprintNode>(
          name.value_or(callback.attr("__name__").cast<std::string>()),
          direction);
      parent.addChild(cylinder);
      callback(cylinder);
      return cylinder;
    }

    std::shared_ptr<CylinderContainerBlueprintNode> enter() {
      if (!name.has_value()) {
        throw std::invalid_argument("Name is required in context manager");
      }
      auto cylinder = std::make_shared<CylinderContainerBlueprintNode>(
          name.value(), direction);
      parent.addChild(cylinder);
      return cylinder;
    }

    void exit(py::object, py::object, py::object) {}
  };

  py::class_<AddCylinderContainerHelper>(m, "_AddCylinderContainerHelper")
      .def("__call__", &AddCylinderContainerHelper::operator())
      .def("__enter__", &AddCylinderContainerHelper::enter)
      .def("__exit__", &AddCylinderContainerHelper::exit);

  struct AddMaterialDesignatorHelper {
    BlueprintNode& parent;

    std::shared_ptr<BlueprintNode> operator()(py::function callback) {
      auto material = std::make_shared<MaterialDesignatorBlueprintNode>();
      parent.addChild(material);
      callback(material);
      return material;
    }

    std::shared_ptr<BlueprintNode> enter() {
      auto material = std::make_shared<MaterialDesignatorBlueprintNode>();
      parent.addChild(material);
      return material;
    }

    void exit(py::object, py::object, py::object) {}
  };

  py::class_<AddMaterialDesignatorHelper>(m, "_AddMaterialDesignatorHelper")
      .def("__call__", &AddMaterialDesignatorHelper::operator())
      .def("__enter__", &AddMaterialDesignatorHelper::enter)
      .def("__exit__", &AddMaterialDesignatorHelper::exit);

  struct AddStaticVolumeHelper {
    BlueprintNode& parent;
    Transform3 transform;
    std::shared_ptr<const VolumeBounds> bounds;
    std::optional<std::string> name;

    std::shared_ptr<StaticBlueprintNode> operator()(py::function callback) {
      std::string callbackName = callback.attr("__name__").cast<std::string>();
      auto node = std::make_shared<Acts::StaticBlueprintNode>(
          std::make_unique<Acts::TrackingVolume>(transform, bounds,
                                                 name.value_or(callbackName)));
      parent.addChild(node);
      callback(node);
      return node;
    }

    std::shared_ptr<StaticBlueprintNode> enter() {
      if (!name.has_value()) {
        throw std::invalid_argument("Name is required in context manager");
      }
      auto node = std::make_shared<Acts::StaticBlueprintNode>(
          std::make_unique<Acts::TrackingVolume>(transform, bounds,
                                                 name.value()));
      parent.addChild(node);
      return node;
    }

    void exit(py::object, py::object, py::object) {}
  };

  py::class_<AddStaticVolumeHelper>(m, "_AddStaticVolumeHelper")
      .def("__call__", &AddStaticVolumeHelper::operator())
      .def("__enter__", &AddStaticVolumeHelper::enter)
      .def("__exit__", &AddStaticVolumeHelper::exit);

  auto blueprintNode =
      py::class_<BlueprintNode, std::shared_ptr<BlueprintNode>>(m,
                                                                "BlueprintNode")
          .def("__str__",
               [](const BlueprintNode& self) {
                 std::stringstream ss;
                 self.toStream(ss);
                 return ss.str();
               })
          .def(
              "addStaticVolume",
              [](BlueprintNode& self, const Transform3& transform,
                 std::shared_ptr<const VolumeBounds> bounds,
                 const std::string& name) {
                auto node = std::make_shared<Acts::StaticBlueprintNode>(
                    std::make_unique<Acts::TrackingVolume>(transform, bounds,
                                                           name));
                self.addChild(node);
                return node;
              },
              py::arg("transform"), py::arg("bounds"),
              py::arg("name") = "undefined")

          .def(
              "StaticVolume",
              [](BlueprintNode& self, const Transform3& transform,
                 std::shared_ptr<const VolumeBounds> bounds,
                 std::optional<std::string> name = std::nullopt) {
                return AddStaticVolumeHelper{self, transform, bounds,
                                             name.value_or("undefined")};
              },
              py::arg("transform"), py::arg("bounds"),
              py::arg("name") = std::nullopt)

          .def(
              "addCylinderContainer",
              [](BlueprintNode& self, const std::string& name,
                 BinningValue direction) {
                auto cylinder =
                    std::make_shared<CylinderContainerBlueprintNode>(name,
                                                                     direction);
                self.addChild(cylinder);
                return cylinder;
              },
              py::arg("name"), py::arg("direction"))

          .def(
              "CylinderContainer",
              [](BlueprintNode& self, BinningValue direction,
                 std::optional<std::string> name = std::nullopt) {
                return AddCylinderContainerHelper{self, direction, name};
              },
              py::arg("direction"), py::arg("name") = std::nullopt)

          .def("Material",
               [](BlueprintNode& self) {
                 return AddMaterialDesignatorHelper{self};
               })

          .def("addChild", &BlueprintNode::addChild)
          .def_property_readonly("children",
                                 py::overload_cast<>(&BlueprintNode::children))
          .def_property_readonly("name", &BlueprintNode::name)
          .def("graphviz",
               [](BlueprintNode& self, py::object fh) {
                 std::stringstream ss;
                 self.graphviz(ss);
                 fh.attr("write")(ss.str());
               })
          .def("visualize", &BlueprintNode::visualize)
          .def(
              "build",
              [](BlueprintNode& self, Logging::Level level) {
                return self.build(*getDefaultLogger("Blueprint", level));
              },
              py::arg("level") = Logging::INFO);

  py::class_<BlueprintNode::MutableChildRange>(blueprintNode,
                                               "MutableChildRange")
      .def(
          "__iter__",
          [](BlueprintNode::MutableChildRange& self) {
            return py::make_iterator(self.begin(), self.end());
          },
          py::keep_alive<0, 1>())
      .def(
          "__getitem__",
          [](BlueprintNode::MutableChildRange& self,
             int i) -> Acts::BlueprintNode& {
            if (i < 0) {
              i += self.size();
            }
            return self.at(i);
          },
          py::return_value_policy::reference_internal)
      .def("__len__", [](const BlueprintNode::MutableChildRange& self) {
        return self.size();
      });

  py::class_<Acts::StaticBlueprintNode, Acts::BlueprintNode,
             std::shared_ptr<Acts::StaticBlueprintNode>>(m,
                                                         "StaticBlueprintNode")
      .def(py::init([](const Transform3& transform,
                       std::shared_ptr<const VolumeBounds> bounds,
                       const std::string& name) {
             return std::make_shared<Acts::StaticBlueprintNode>(
                 std::make_unique<Acts::TrackingVolume>(transform, bounds,
                                                        name));
           }),
           py::arg("transform"), py::arg("bounds"),
           py::arg("name") = "undefined");

  auto cylNode =
      py::class_<Acts::CylinderContainerBlueprintNode, Acts::BlueprintNode,
                 std::shared_ptr<Acts::CylinderContainerBlueprintNode>>(
          m, "CylinderContainerBlueprintNode")
          .def(py::init<const std::string&, BinningValue,
                        CylinderVolumeStack::AttachmentStrategy,
                        CylinderVolumeStack::ResizeStrategy>(),
               py::arg("name"), py::arg("direction"),
               py::arg("attachmentStrategy") =
                   CylinderVolumeStack::AttachmentStrategy::Gap,
               py::arg("resizeStrategy") =
                   CylinderVolumeStack::ResizeStrategy::Gap)
          .def_property(
              "attachmentStrategy",
              [](Acts::CylinderContainerBlueprintNode& self) {
                return self.attachmentStrategy();
              },
              [](Acts::CylinderContainerBlueprintNode& self,
                 CylinderVolumeStack::AttachmentStrategy strategy) {
                self.setAttachmentStrategy(strategy);
              })
          .def_property(
              "resizeStrategy",
              [](Acts::CylinderContainerBlueprintNode& self) {
                return self.resizeStrategy();
              },
              [](Acts::CylinderContainerBlueprintNode& self,
                 CylinderVolumeStack::ResizeStrategy strategy) {
                self.setResizeStrategy(strategy);
              })
          .def_property(
              "direction",
              [](Acts::CylinderContainerBlueprintNode& self) {
                return self.direction();
              },
              [](Acts::CylinderContainerBlueprintNode& self,
                 BinningValue value) { self.setDirection(value); });
}

}  // namespace Acts::Python
