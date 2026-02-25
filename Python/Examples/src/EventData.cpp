// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/Track.hpp"
#include "ActsPython/Utilities/WhiteBoardTypeRegistry.hpp"

#include <pybind11/pybind11.h>
namespace py = pybind11;

using namespace ActsExamples;

namespace ActsPython {

void addEventData(py::module& mex) {
  auto constTrackProxy =
      py::class_<ConstTrackProxy>(mex, "ConstTrackProxy")
          .def_property_readonly("index", &ConstTrackProxy::index)
          .def_property_readonly("tipIndex", &ConstTrackProxy::tipIndex)
          .def_property_readonly("stemIndex", &ConstTrackProxy::stemIndex)
          .def_property_readonly("referenceSurface",
                                 &ConstTrackProxy::referenceSurface)
          .def_property_readonly("hasReferenceSurface",
                                 &ConstTrackProxy::hasReferenceSurface)
          // Convert the parameters to a BoundVector so they're accessible by
          // value
          .def_property_readonly("parameters",
                                 [](const ConstTrackProxy& self) {
                                   return Acts::BoundVector{self.parameters()};
                                 })
          // Convert the covariance to a BoundMatrix so they're accessible by
          // value
          .def_property_readonly("covariance",
                                 [](const ConstTrackProxy& self) {
                                   return Acts::BoundMatrix{self.covariance()};
                                 })
          .def_property_readonly("particleHypothesis",
                                 &ConstTrackProxy::particleHypothesis);

  auto constTrackContainer =
      py::class_<ConstTrackContainer>(mex, "ConstTrackContainer")
          .def("__len__", &ConstTrackContainer::size)
          .def("__iter__",
               [](const ConstTrackContainer& self) {
                 return py::make_iterator(self.begin(), self.end());
               })
      // .def("__getitem__", [](const TrackContainer& self, size_t index) {
      //   return self.getTrack(index);
      // })
      ;

  WhiteBoardRegistry::registerClass(constTrackContainer);
}

}  // namespace ActsPython
