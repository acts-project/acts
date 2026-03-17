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
#include "ActsPython/Utilities/WhiteBoardRegistry.hpp"

#include <array>
#include <memory>
#include <optional>
#include <span>
#include <string>
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
}

/// Zero-copy 2D view over a column of std::array<float, Cols>.
/// Returns shape (N, Cols), dtype float32, read-only.
/// Throws if the container does not have the required column.
template <std::size_t Cols>
using ArrayColumnGetter = ConstSpacePointColumnProxy<std::array<float, Cols>> (
    SpacePointContainer2::*)() const;

template <std::size_t Cols>
auto arrayColumn(ArrayColumnGetter<Cols> getColumn,
                 SpacePointColumns requiredColumn,
                 const std::string_view& columnName) {
  return [getColumn, requiredColumn,
          columnName](const SpacePointContainer2& self) {
    if (!self.hasColumns(requiredColumn)) {
      throw py::attribute_error(
          std::format("SpacePointContainer2 does not have "
                      "the {} column",
                      columnName));
    }
    const auto nRows = static_cast<py::ssize_t>(self.size());
    if (nRows == 0) {
      auto arr = py::array_t<float>(
          std::vector<py::ssize_t>{0, static_cast<py::ssize_t>(Cols)});
      arr.attr("flags").attr("writeable") = py::bool_(false);
      return arr;
    }
    const auto col = (self.*getColumn)();
    const auto& data = col.data();
    constexpr py::ssize_t rowStride =
        static_cast<py::ssize_t>(Cols * sizeof(float));
    constexpr py::ssize_t colStride = sizeof(float);
    auto arr = py::array_t<float>(
        {nRows, static_cast<py::ssize_t>(Cols)}, {rowStride, colStride},
        reinterpret_cast<const float*>(data.data()), py::cast(self));
    arr.attr("flags").attr("writeable") = py::bool_(false);
    return arr;
  };
}

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

  using FloatColumnGetter =
      ConstSpacePointColumnProxy<float> (SpacePointContainer2::*)() const;
  auto floatColumn = [](FloatColumnGetter column) {
    return [column](const SpacePointContainer2& self) {
      return spanToNumpy1d((self.*column)().data(), py::cast(self));
    };
  };

  // SpacePointContainer2
  auto spc2 =
      py::classh<SpacePointContainer2>(m, "SpacePointContainer2")
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
          .def("__iter__",
               [](const SpacePointContainer2& self) {
                 return py::make_iterator(self.begin(), self.end());
               })
          .def_property_readonly("x",
                                 floatColumn(&SpacePointContainer2::xColumn))
          .def_property_readonly("y",
                                 floatColumn(&SpacePointContainer2::yColumn))
          .def_property_readonly("z",
                                 floatColumn(&SpacePointContainer2::zColumn))
          .def_property_readonly("r",
                                 floatColumn(&SpacePointContainer2::rColumn))
          .def_property_readonly("phi",
                                 floatColumn(&SpacePointContainer2::phiColumn))
          .def_property_readonly("time",
                                 floatColumn(&SpacePointContainer2::timeColumn))
          .def_property_readonly(
              "varianceZ", floatColumn(&SpacePointContainer2::varianceZColumn))
          .def_property_readonly(
              "varianceR", floatColumn(&SpacePointContainer2::varianceRColumn))
          .def_property_readonly(
              "xyColumn", arrayColumn<2>(&SpacePointContainer2::xyColumn,
                                         SpacePointColumns::PackedXY, "xy"))
          .def_property_readonly(
              "zrColumn", arrayColumn<2>(&SpacePointContainer2::zrColumn,
                                         SpacePointColumns::PackedZR, "zr"))
          .def_property_readonly(
              "xy", arrayColumn<2>(&SpacePointContainer2::xyColumn,
                                   SpacePointColumns::PackedXY, "xy"))
          .def_property_readonly(
              "zr", arrayColumn<2>(&SpacePointContainer2::zrColumn,
                                   SpacePointColumns::PackedZR, "zr"))
          .def_property_readonly(
              "xyz", arrayColumn<3>(&SpacePointContainer2::xyzColumn,
                                    SpacePointColumns::PackedXYZ, "xyz"))
          .def_property_readonly(
              "xyzr", arrayColumn<4>(&SpacePointContainer2::xyzrColumn,
                                     SpacePointColumns::PackedXYZR, "xyzr"))
          .def_property_readonly(
              "varianceZRColumn",
              arrayColumn<2>(&SpacePointContainer2::varianceZRColumn,
                             SpacePointColumns::PackedVarianceZR, "varianceZR"))
          .def_property_readonly(
              "topStripVectorColumn",
              arrayColumn<3>(&SpacePointContainer2::topStripVectorColumn,
                             SpacePointColumns::TopStripVector,
                             "topStripVector"));

  WhiteBoardRegistry::registerClass(spc2);

  // ConstSpacePointProxy2
  // Note: SpacePointProxy2 has both mutable (requires(!ReadOnly)) and const
  // overloads of x/y/z etc. with different return types. GCC 13 includes the
  // constrained-out overloads in pointer-to-member resolution, so we must
  // static_cast to the exact return type to disambiguate.
  using FloatGetter = float (ConstSpacePointProxy2::*)() const noexcept;
  py::class_<ConstSpacePointProxy2>(m, "ConstSpacePointProxy2")
      .def_property_readonly("index", &ConstSpacePointProxy2::index)
      .def_property_readonly(
          "x", static_cast<FloatGetter>(&ConstSpacePointProxy2::x))
      .def_property_readonly(
          "y", static_cast<FloatGetter>(&ConstSpacePointProxy2::y))
      .def_property_readonly(
          "z", static_cast<FloatGetter>(&ConstSpacePointProxy2::z))
      .def_property_readonly(
          "r", static_cast<FloatGetter>(&ConstSpacePointProxy2::r))
      .def_property_readonly(
          "phi", static_cast<FloatGetter>(&ConstSpacePointProxy2::phi))
      .def_property_readonly(
          "time", static_cast<FloatGetter>(&ConstSpacePointProxy2::time))
      .def_property_readonly(
          "varianceZ",
          static_cast<FloatGetter>(&ConstSpacePointProxy2::varianceZ))
      .def_property_readonly(
          "varianceR",
          static_cast<FloatGetter>(&ConstSpacePointProxy2::varianceR));

  // SeedContainer2
  auto seedContainer2 =
      py::classh<SeedContainer2>(m, "SeedContainer2")
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

  WhiteBoardRegistry::registerClass(seedContainer2);

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
