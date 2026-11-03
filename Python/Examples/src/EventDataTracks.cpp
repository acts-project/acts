// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/AnyTrackProxy.hpp"
#include "Acts/EventData/AnyTrackStateProxy.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "ActsExamples/EventData/Track.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "EventDataHelpers.hpp"

namespace py = pybind11;

using namespace ActsExamples;
using namespace ActsPython::EventDataHelpers;

namespace {

/// Bind the surface shared by the concrete and Any track proxies: all four
/// (`ConstTrackProxy`/`TrackProxy`/`AnyConstTrackProxy`/`AnyMutableTrackProxy`)
/// publicly inherit the same `TrackProxyCommon` base.
template <typename Cls>
void bindTrackProxyCommon(Cls& cls) {
  using Proxy = typename Cls::type;
  cls.def_property_readonly("index", &Proxy::index)
      .def_property_readonly("referenceSurface", &Proxy::referenceSurface)
      .def_property_readonly("hasReferenceSurface", &Proxy::hasReferenceSurface)
      .def_property_readonly("particleHypothesis", &Proxy::particleHypothesis)
      .def_property_readonly("nTrackStates", &Proxy::nTrackStates)
      .def_property_readonly("isForwardLinked", &Proxy::isForwardLinked)
      .def_property_readonly("loc0", &Proxy::loc0)
      .def_property_readonly("loc1", &Proxy::loc1)
      .def_property_readonly("phi", &Proxy::phi)
      .def_property_readonly("theta", &Proxy::theta)
      .def_property_readonly("time", &Proxy::time)
      .def_property_readonly("qOverP", &Proxy::qOverP)
      .def_property_readonly("charge", &Proxy::charge)
      .def_property_readonly("absoluteMomentum", &Proxy::absoluteMomentum)
      .def_property_readonly("transverseMomentum", &Proxy::transverseMomentum)
      .def_property_readonly("direction", &Proxy::direction)
      .def_property_readonly("momentum", &Proxy::momentum)
      .def_property_readonly("fourMomentum", &Proxy::fourMomentum)
      .def(
          "hasColumn",
          [](const Proxy& self, const std::string& key) {
            return self.hasColumn(Acts::hashStringDynamic(key));
          },
          py::arg("key"))
      .def_property_readonly("trackStates",
                             py::cpp_function(
                                 [](Proxy& self) {
                                   auto range = self.trackStates();
                                   return py::make_iterator(range.begin(),
                                                            range.end());
                                 },
                                 py::keep_alive<0, 1>()))
      .def_property_readonly("trackStatesReversed",
                             py::cpp_function(
                                 [](Proxy& self) {
                                   auto range = self.trackStatesReversed();
                                   return py::make_iterator(range.begin(),
                                                            range.end());
                                 },
                                 py::keep_alive<0, 1>()))
      .def_property_readonly(
          "innermostTrackState",
          [](Proxy& self) { return self.innermostTrackState(); })
      .def_property_readonly("outermostTrackState", [](Proxy& self) {
        return self.outermostTrackState();
      });

  // tipIndex/stemIndex/nMeasurements/nHoles/nOutliers/nSharedHits/nSplitHits/
  // chi2/nDoF are overloaded on const (by value) vs. mutable (by reference);
  // the accessor lambdas below select the overload via self's constness.
  bindScalar(cls, "tipIndex",
             [](auto& self) -> decltype(auto) { return self.tipIndex(); });
  bindScalar(cls, "stemIndex",
             [](auto& self) -> decltype(auto) { return self.stemIndex(); });
  bindScalar(cls, "nMeasurements",
             [](auto& self) -> decltype(auto) { return self.nMeasurements(); });
  bindScalar(cls, "nHoles",
             [](auto& self) -> decltype(auto) { return self.nHoles(); });
  bindScalar(cls, "nOutliers",
             [](auto& self) -> decltype(auto) { return self.nOutliers(); });
  bindScalar(cls, "nSharedHits",
             [](auto& self) -> decltype(auto) { return self.nSharedHits(); });
  bindScalar(cls, "nSplitHits",
             [](auto& self) -> decltype(auto) { return self.nSplitHits(); });
  bindScalar(cls, "chi2",
             [](auto& self) -> decltype(auto) { return self.chi2(); });
  bindScalar(cls, "nDoF",
             [](auto& self) -> decltype(auto) { return self.nDoF(); });
  bindEigen(cls, "parameters",
            [](auto& self) -> decltype(auto) { return self.parameters(); });
  bindEigen(cls, "covariance",
            [](auto& self) -> decltype(auto) { return self.covariance(); });
}

/// Bind the extra surface only the concrete `ConstTrackProxy`/`TrackProxy`
/// have: `AnyTrackProxy` is a value-mutate handle, not a track-building API,
/// so it lacks even the `referenceSurface`/`particleHypothesis` setters.
template <typename Cls>
void bindTrackProxyVector(Cls& cls) {
  using Proxy = typename Cls::type;
  cls.def_property_readonly("createParametersAtReference",
                            &Proxy::createParametersAtReference);

  if constexpr (!Proxy::ReadOnly) {
    cls.def_property("referenceSurface", &Proxy::referenceSurface,
                     &Proxy::setReferenceSurface)
        .def_property("particleHypothesis", &Proxy::particleHypothesis,
                      &Proxy::setParticleHypothesis)
        .def("appendTrackState", &Proxy::appendTrackState,
             py::arg("mask") = Acts::TrackStatePropMask::All,
             py::keep_alive<0, 1>())
        .def("linkForward", &Proxy::linkForward)
        .def("reverseTrackStates", &Proxy::reverseTrackStates,
             py::arg("invertJacobians") = false)
        // copyFrom*() are member function templates, so explicit
        // instantiation per source type replaces what would otherwise need a
        // lambda (a template has no single address to bind).
        .def("copyFrom", &Proxy::template copyFrom<ConstTrackProxy>)
        .def("copyFrom", &Proxy::template copyFrom<TrackProxy>)
        .def("copyFromWithoutStates",
             &Proxy::template copyFromWithoutStates<ConstTrackProxy>)
        .def("copyFromWithoutStates",
             &Proxy::template copyFromWithoutStates<TrackProxy>)
        .def("copyFromShallow",
             &Proxy::template copyFromShallow<ConstTrackProxy>)
        .def("copyFromShallow", &Proxy::template copyFromShallow<TrackProxy>);
  }
}

}  // namespace

namespace ActsPython {

void addEventDataTracks(py::module& mex) {
  py::class_<ConstTrackProxy> constTrackProxy(mex, "ConstTrackProxy");
  bindTrackProxyCommon(constTrackProxy);
  bindTrackProxyVector(constTrackProxy);

  py::class_<TrackProxy> trackProxy(mex, "TrackProxy");
  bindTrackProxyCommon(trackProxy);
  bindTrackProxyVector(trackProxy);

  py::classh<Acts::AnyConstTrackProxy> anyConstTrack(mex, "AnyConstTrackProxy");
  bindTrackProxyCommon(anyConstTrack);
  anyConstTrack.def(py::init<ConstTrackProxy>());
  py::implicitly_convertible<ConstTrackProxy, Acts::AnyConstTrackProxy>();

  py::classh<Acts::AnyMutableTrackProxy> anyMutableTrack(
      mex, "AnyMutableTrackProxy");
  bindTrackProxyCommon(anyMutableTrack);
  anyMutableTrack.def(py::init<TrackProxy>());
  py::implicitly_convertible<TrackProxy, Acts::AnyMutableTrackProxy>();
}

}  // namespace ActsPython
