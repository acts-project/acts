// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/Types.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SeedSpacePointSelection.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsPython/Utilities/ProxyTether.hpp"
#include "ActsPython/Utilities/WhiteBoardRegistry.hpp"

#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

// Prevent stl.h's type_caster<std::vector<T>> from matching, which would
// break py::cast<std::unique_ptr<T>> needed by WhiteBoardRegistry (full
// specializations take priority over stl.h's partial one).
PYBIND11_MAKE_OPAQUE(ActsExamples::ClusterContainer)
PYBIND11_MAKE_OPAQUE(ActsExamples::ProtoTrackContainer)
// MeasurementSimHitsMap == SimHitMeasurementsMap at the C++ level (both are
// flat_multimap<std::uint32_t, std::uint32_t>), so only one MAKE_OPAQUE is
// needed.
PYBIND11_MAKE_OPAQUE(ActsExamples::MeasurementSimHitsMap)
PYBIND11_MAKE_OPAQUE(ActsExamples::MeasurementParticlesMap)
PYBIND11_MAKE_OPAQUE(ActsExamples::ParticleMeasurementsMap)

namespace py = pybind11;

using namespace ActsExamples;

namespace {

template <typename Map>
auto bindFlatMultimap(py::module& m, const char* name) {
  // Register the iterator type before the __iter__ binding that returns it.
  const std::string iterName = std::string(name) + "Iterator";
  ActsPython::bindIteratorTether<Map>(m, iterName.c_str());

  auto cls = py::classh<Map>(m, name)
                 .def(py::init<>())
                 .def("__len__", &Map::size)
                 .def("__iter__",
                      [](py::object self) {
                        const auto& map = self.cast<const Map&>();
                        return ActsPython::IteratorTether<Map>{
                            self, py::make_iterator(map.begin(), map.end())};
                      })
                 .def("__contains__",
                      [](const Map& self, const typename Map::key_type& key) {
                        return self.find(key) != self.end();
                      })
                 .def(
                     "valuesFor",
                     [](const Map& self, const typename Map::key_type& key) {
                       auto [first, last] = self.equal_range(key);
                       std::vector<typename Map::mapped_type> result;
                       for (auto it = first; it != last; ++it) {
                         result.push_back(it->second);
                       }
                       return result;
                     },
                     py::arg("key"))
                 .def(
                     "insert",
                     [](Map& self, const typename Map::key_type& key,
                        const typename Map::mapped_type& value) {
                       self.emplace(key, value);
                     },
                     py::arg("key"), py::arg("value"));
  ActsPython::WhiteBoardRegistry::registerClass(cls);
  return cls;
}

template <typename T>
void bindIndexMultimapPair(py::module& m, const char* forwardName,
                           const char* inverseName) {
  // Bind inverse first so its Python type is registered before being used as
  // the return type of .invert() on the forward map.
  bindFlatMultimap<ActsExamples::InverseMultimap<T>>(m, inverseName);
  auto fwd = bindFlatMultimap<ActsExamples::IndexMultimap<T>>(m, forwardName);
  fwd.def("invert", [](const ActsExamples::IndexMultimap<T>& self) {
    return ActsExamples::invertIndexMultimap(self);
  });
}

}  // namespace

namespace ActsPython {

// The event data bindings are split across several translation units to
// bound the peak compiler memory of each; they are registered here in
// dependency order.
void addEventDataTrackStates(py::module& mex);
void addEventDataTracks(py::module& mex);
void addEventDataTrackContainers(py::module& mex);
void addEventDataMeasurements(py::module& mex);

void addEventData(py::module& mex) {
  // TrackStatePropMask is bound first: several TrackProxy/TrackStateProxy
  // methods below use it as a default argument, which requires the type
  // caster to already be registered when those .def() calls execute.
  py::enum_<Acts::TrackStatePropMask>(mex, "TrackStatePropMask")
      .value("None_", Acts::TrackStatePropMask::None)
      .value("Predicted", Acts::TrackStatePropMask::Predicted)
      .value("Filtered", Acts::TrackStatePropMask::Filtered)
      .value("Smoothed", Acts::TrackStatePropMask::Smoothed)
      .value("Jacobian", Acts::TrackStatePropMask::Jacobian)
      .value("Calibrated", Acts::TrackStatePropMask::Calibrated)
      .value("All", Acts::TrackStatePropMask::All)
      .def("__or__", [](Acts::TrackStatePropMask a,
                        Acts::TrackStatePropMask b) { return a | b; })
      .def("__and__", [](Acts::TrackStatePropMask a,
                         Acts::TrackStatePropMask b) { return a & b; });

  addEventDataTrackStates(mex);
  addEventDataTracks(mex);
  addEventDataTrackContainers(mex);

  mex.attr("kTrackIndexInvalid") = Acts::kTrackIndexInvalid;

  py::bind_vector<ProtoTrack>(mex, "ProtoTrack");

  auto protoTrackContainer =
      py::bind_vector<ProtoTrackContainer, py::smart_holder>(
          mex, "ProtoTrackContainer");
  WhiteBoardRegistry::registerClass(protoTrackContainer);

  py::class_<Cluster>(mex, "Cluster")
      .def(py::init<>())
      .def_readwrite("sizeLoc0", &Cluster::sizeLoc0)
      .def_readwrite("sizeLoc1", &Cluster::sizeLoc1)
      .def_readwrite("globalPosition", &Cluster::globalPosition)
      .def_readwrite("localDirection", &Cluster::localDirection)
      .def_readwrite("lengthDirection", &Cluster::lengthDirection)
      .def_readwrite("localEta", &Cluster::localEta)
      .def_readwrite("localPhi", &Cluster::localPhi)
      .def_readwrite("globalEta", &Cluster::globalEta)
      .def_readwrite("globalPhi", &Cluster::globalPhi)
      .def_readwrite("etaAngle", &Cluster::etaAngle)
      .def_readwrite("phiAngle", &Cluster::phiAngle)
      .def("sumActivations", &Cluster::sumActivations);

  auto clusterContainer = py::bind_vector<ClusterContainer, py::smart_holder>(
      mex, "ClusterContainer");
  WhiteBoardRegistry::registerClass(clusterContainer);

  addEventDataMeasurements(mex);

  // MeasurementSimHitsMap and SimHitMeasurementsMap are the same C++ type
  // (flat_multimap<std::uint32_t, std::uint32_t>) because SimHitIndex == Index
  // == std::uint32_t. Bind once and alias the second name.
  auto simHitsMap =
      bindFlatMultimap<MeasurementSimHitsMap>(mex, "MeasurementSimHitsMap");
  simHitsMap.def("invert", [](const MeasurementSimHitsMap& self) {
    return invertIndexMultimap(self);
  });
  mex.attr("SimHitMeasurementsMap") = mex.attr("MeasurementSimHitsMap");

  bindIndexMultimapPair<SimBarcode>(mex, "MeasurementParticlesMap",
                                    "ParticleMeasurementsMap");

  py::enum_<SeedSpacePointSelection>(mex, "SeedSpacePointSelection")
      .value("FirstThree", SeedSpacePointSelection::FirstThree)
      .value("InnermostTriplet", SeedSpacePointSelection::InnermostTriplet)
      .value("SpreadTriplet", SeedSpacePointSelection::SpreadTriplet)
      .value("All", SeedSpacePointSelection::All);
}

}  // namespace ActsPython
