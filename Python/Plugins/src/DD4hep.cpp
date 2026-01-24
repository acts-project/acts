// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsPlugins/DD4hep/DD4hepDetectorElement.hpp"
#include "ActsPlugins/DD4hep/DD4hepFieldAdapter.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>
#include <string>

#include <DD4hep/DetElement.h>
#include <DD4hep/Fields.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(ActsPluginsPythonBindingsDD4hep, dd4hep) {
  using namespace Acts;
  using namespace ActsPlugins;
  using namespace ActsPython;

  // Basic bindings
  {
    // The DD4hep Detector Element
    py::class_<dd4hep::DetElement, std::shared_ptr<dd4hep::DetElement>>(
        dd4hep, "DD4hepDetElement");

    // The Acts:: glue detector element
    py::class_<DD4hepDetectorElement, DetectorElementBase,
               std::shared_ptr<DD4hepDetectorElement>>(dd4hep,
                                                       "DD4hepDetectorElement");

    py::class_<DD4hepFieldAdapter, MagneticFieldProvider,
               std::shared_ptr<DD4hepFieldAdapter>>(dd4hep,
                                                    "DD4hepFieldAdapter");
  }

  // Helper method
  {
    dd4hep.def(
        "createDD4hepIdGeoIdMap",
        [](const TrackingGeometry& tGeometry)
            -> std::map<DD4hepDetectorElement::DD4hepVolumeID,
                        GeometryIdentifier> {
          // The surface visitor
          struct DD4hepIdGrabber {
            std::map<DD4hepDetectorElement::DD4hepVolumeID, GeometryIdentifier>
                dd4hepIdGeoIdMap;

            void operator()(const Surface* surface) {
              const auto* dde = surface->associatedDetectorElement();
              const auto* dd4hepDetElement =
                  dynamic_cast<const DD4hepDetectorElement*>(dde);
              // Check if it is valid
              if (dd4hepDetElement != nullptr) {
                dd4hep::DDSegmentation::VolumeID dd4hepID =
                    dd4hepDetElement->sourceElement().volumeID();
                auto geoID = surface->geometryId();
                dd4hepIdGeoIdMap[dd4hepID] = geoID;
              }
            }
          };

          // Create an instance
          DD4hepIdGrabber dd4hepIdGrabber;
          // Visit the surfaces & return what you have
          tGeometry.visitSurfaces(dd4hepIdGrabber);
          return dd4hepIdGrabber.dd4hepIdGeoIdMap;
        });
  }
}
