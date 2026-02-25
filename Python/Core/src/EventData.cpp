// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/GenericFreeTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/SeedContainer2.hpp"
#include "Acts/EventData/SeedProxy2.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/SpacePointColumns.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/EventData/SpacePointProxy2.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <array>
#include <memory>
#include <optional>
#include <span>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;

namespace ActsPython {

template <typename T>
auto spanToNumpy1d(std::span<T> s, const py::object& base) {
  using type = std::remove_cvref_t<T>;
  return py::array_t<type>(
      {static_cast<py::ssize_t>(s.size())},      // shape
      {static_cast<py::ssize_t>(sizeof(type))},  // strides (bytes)
      s.data(),                                  // data ptr
      base                                       // base/owner
  );
};

/// @brief This adds the classes from Core/EventData to the python module
/// @param m the pybind11 core module
void addEventData(py::module_& m) {
  // SpacePointColumns enum
  py::enum_<SpacePointColumns>(m, "SpacePointColumns")
      .value("None", SpacePointColumns::None)
      .value("SourceLinks", SpacePointColumns::SourceLinks)
      .value("X", SpacePointColumns::X)
      .value("Y", SpacePointColumns::Y)
      .value("Z", SpacePointColumns::Z)
      .value("R", SpacePointColumns::R)
      .value("Phi", SpacePointColumns::Phi)
      .value("Time", SpacePointColumns::Time)
      .value("VarianceZ", SpacePointColumns::VarianceZ)
      .value("VarianceR", SpacePointColumns::VarianceR)
      .value("TopStripVector", SpacePointColumns::TopStripVector)
      .value("BottomStripVector", SpacePointColumns::BottomStripVector)
      .value("StripCenterDistance", SpacePointColumns::StripCenterDistance)
      .value("TopStripCenter", SpacePointColumns::TopStripCenter)
      .value("CopyFromIndex", SpacePointColumns::CopyFromIndex)
      .value("PackedXY", SpacePointColumns::PackedXY)
      .value("PackedZR", SpacePointColumns::PackedZR)
      .value("PackedXYZ", SpacePointColumns::PackedXYZ)
      .value("PackedXYZR", SpacePointColumns::PackedXYZR)
      .value("PackedVarianceZR", SpacePointColumns::PackedVarianceZR)
      .value("Strip", SpacePointColumns::Strip)
      .value("All", SpacePointColumns::All)
      .def("__or__",
           [](SpacePointColumns a, SpacePointColumns b) {
             return static_cast<SpacePointColumns>(
                 static_cast<std::uint32_t>(a) | static_cast<std::uint32_t>(b));
           })
      .def("__and__", [](SpacePointColumns a, SpacePointColumns b) {
        return static_cast<SpacePointColumns>(static_cast<std::uint32_t>(a) &
                                              static_cast<std::uint32_t>(b));
      });

  // SpacePointContainer2
  auto spc2 =
      py::class_<SpacePointContainer2>(m, "SpacePointContainer2")
          .def(py::init<SpacePointColumns>(),
               py::arg("columns") = SpacePointColumns::None)
          .def_property_readonly("size", &SpacePointContainer2::size)
          .def_property_readonly("empty", &SpacePointContainer2::empty)
          .def_property_readonly("hasColumns",
                                 &SpacePointContainer2::hasColumns)
          .def("reserve", &SpacePointContainer2::reserve, py::arg("size"),
               py::arg("averageSourceLinks") = 1)
          .def("clear", &SpacePointContainer2::clear)
          .def("__len__", &SpacePointContainer2::size)
          .def("__getitem__",
               [](const SpacePointContainer2& self, SpacePointIndex2 idx) {
                 return ConstSpacePointProxy2(self, idx);
               })
          .def("__iter__", [](const SpacePointContainer2& self) {
            return py::make_iterator(self.begin(), self.end());
          });

  // ConstSpacePointProxy2
  py::class_<ConstSpacePointProxy2>(m, "ConstSpacePointProxy2")
      .def_property_readonly("index", &ConstSpacePointProxy2::index)
      .def_property_readonly("x", &ConstSpacePointProxy2::x)
      .def_property_readonly("y", &ConstSpacePointProxy2::y)
      .def_property_readonly("z", &ConstSpacePointProxy2::z)
      .def_property_readonly("r", &ConstSpacePointProxy2::r)
      .def_property_readonly("phi", &ConstSpacePointProxy2::phi)
      .def_property_readonly("time", &ConstSpacePointProxy2::time)
      .def_property_readonly("varianceZ", &ConstSpacePointProxy2::varianceZ)
      .def_property_readonly("varianceR", &ConstSpacePointProxy2::varianceR);

  // SeedContainer2
  py::class_<SeedContainer2>(m, "SeedContainer2")
      .def(py::init<>())
      .def_property_readonly("size", &SeedContainer2::size)
      .def_property_readonly("empty", &SeedContainer2::empty)
      .def("__len__", &SeedContainer2::size)
      .def("__getitem__",
           [](const SeedContainer2& self, SeedIndex2 idx) {
             return ConstSeedProxy2(self, idx);
           })
      .def("__iter__", [](const SeedContainer2& self) {
        return py::make_iterator(self.begin(), self.end());
      });

  // ConstSeedProxy2
  py::class_<ConstSeedProxy2>(m, "ConstSeedProxy2")
      .def_property_readonly("index", &ConstSeedProxy2::index)
      .def_property_readonly("size", &ConstSeedProxy2::size)
      .def_property_readonly("empty", &ConstSeedProxy2::empty)
      .def_property_readonly("quality", &ConstSeedProxy2::quality)
      .def_property_readonly("vertexZ", &ConstSeedProxy2::vertexZ)
      .def_property_readonly(
          "spacePointIndices", [](const ConstSeedProxy2& self) {
            return spanToNumpy1d(self.spacePointIndices(), py::cast(self));
          });

  // BoundTrackParameters (alias for
  // GenericBoundTrackParameters<ParticleHypothesis>)
  py::class_<BoundTrackParameters>(m, "BoundTrackParameters")
      .def_property_readonly("parameters",
                             [](const BoundTrackParameters& self) {
                               return BoundVector(self.parameters());
                             })
      .def_property_readonly(
          "covariance",
          [](const BoundTrackParameters& self) { return self.covariance(); })
      .def_property_readonly("referenceSurface",
                             [](const BoundTrackParameters& self)
                                 -> std::shared_ptr<const Surface> {
                               return self.referenceSurface().getSharedPtr();
                             })
      .def_property_readonly("particleHypothesis",
                             &BoundTrackParameters::particleHypothesis)
      .def_property_readonly("position", &BoundTrackParameters::position)
      .def_property_readonly("momentum", &BoundTrackParameters::momentum)
      .def_property_readonly("transverseMomentum",
                             &BoundTrackParameters::transverseMomentum);

  // SourceLink - type-erased wrapper (minimal binding for type identity)
  // Note: SourceLink has no default constructor; concrete types (e.g.
  // IndexSourceLink from Examples) are typically used to construct it.
  py::class_<SourceLink>(m, "SourceLink");
}

}  // namespace ActsPython
