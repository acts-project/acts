// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/SeedContainer2.hpp"
#include "Acts/EventData/SeedProxy2.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/SpacePointColumns.hpp"
#include "Acts/EventData/SpacePointContainer2.hpp"
#include "Acts/EventData/SpacePointProxy2.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsPython/Utilities/ProxyTether.hpp"
#include "ActsPython/Utilities/WhiteBoardRegistry.hpp"

#include <array>
#include <format>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
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

/// Build the "missing column" error message for a space point accessor.
inline std::string missingColumnMessage(const std::string_view& columnName) {
  return std::format(
      "SpacePoint has no '{}' column; construct the container with "
      "the corresponding SpacePointColumns flag",
      columnName);
}

/// Wrap a read accessor on a tethered space point proxy. First validates the
/// tether (raises py::value_error if the container was consumed/disowned), then
/// checks the column exists (raises py::attribute_error if missing), then calls
/// `access` with the inner proxy. `Tether` is given explicitly at the call
/// site; `access` is a callable taking `const Proxy&` (calling the accessor on
/// a const reference avoids the const/non-const overload ambiguity that would
/// otherwise force static_cast on plain member pointers).
template <typename Tether, typename Access>
auto guardedRead(SpacePointColumns requiredColumn,
                 const std::string_view& columnName, Access access) {
  return [requiredColumn, columnName, access](const Tether& self) {
    const auto& proxy = self.checked();
    if (!proxy.container().hasColumns(requiredColumn)) {
      throw py::attribute_error(missingColumnMessage(columnName));
    }
    return access(proxy);
  };
}

/// Wrap a write accessor on a tethered mutable space point proxy: validate the
/// tether, then check the column exists, then assign. `V` is the concrete value
/// type and is given explicitly at the call site: the returned setter takes
/// `const V&` (never auto) so pybind11 can derive the Python->C++ argument
/// conversion. `assign` is a callable taking `Proxy&` and `const V&`.
template <typename Tether, typename V, typename Assign>
auto guardedWrite(SpacePointColumns requiredColumn,
                  const std::string_view& columnName, Assign assign) {
  return [requiredColumn, columnName, assign](Tether& self, const V& value) {
    auto& proxy = self.checked();
    if (!proxy.container().hasColumns(requiredColumn)) {
      throw py::attribute_error(missingColumnMessage(columnName));
    }
    assign(proxy, value);
  };
}

/// @brief This adds the classes from Core/EventData to the python module
/// @param m the pybind11 core module
void addEventData(py::module_& m) {
  // SpacePointColumns enum
  py::enum_<SpacePointColumns>(m, "SpacePointColumns")
      .value("None", SpacePointColumns::None)
      .value("SourceLinks", SpacePointColumns::SourceLinks)
      .value("CopiedFromIndex", SpacePointColumns::CopiedFromIndex)
      .value("X", SpacePointColumns::X)
      .value("Y", SpacePointColumns::Y)
      .value("Z", SpacePointColumns::Z)
      .value("R", SpacePointColumns::R)
      .value("Phi", SpacePointColumns::Phi)
      .value("Time", SpacePointColumns::Time)
      .value("VarianceZ", SpacePointColumns::VarianceZ)
      .value("VarianceR", SpacePointColumns::VarianceR)
      .value("VarianceT", SpacePointColumns::VarianceT)
      .value("StripCalibrationDetails",
             SpacePointColumns::StripCalibrationDetails)
      .value("PackedXY", SpacePointColumns::PackedXY)
      .value("PackedZR", SpacePointColumns::PackedZR)
      .value("PackedXYZ", SpacePointColumns::PackedXYZ)
      .value("PackedXYZR", SpacePointColumns::PackedXYZR)
      .value("PackedVarianceZR", SpacePointColumns::PackedVarianceZR)
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

  // Proxies are bound as ProxyTether under the same Python names, so isinstance
  // is preserved; see ProxyTether.hpp for the disown/keep-alive rationale. The
  // aliases below keep the binding call sites readable.
  using ConstSpTether =
      ProxyTether<ConstSpacePointProxy2, SpacePointContainer2>;
  using MutSpTether =
      ProxyTether<MutableSpacePointProxy2, SpacePointContainer2>;
  using ConstSeedTether = ProxyTether<ConstSeedProxy2, SeedContainer2>;
  using MutSeedTether = ProxyTether<MutableSeedProxy2, SeedContainer2>;

  // Register iterator types before the __iter__ bindings that return them.
  bindIndexIteratorTether<SpacePointContainer2>(
      m, "_SpacePointContainer2Iterator");
  bindIndexIteratorTether<SeedContainer2>(m, "_SeedContainer2Iterator");

  using FloatColumnGetter =
      ConstSpacePointColumnProxy<float> (SpacePointContainer2::*)() const;
  auto floatColumn = [](FloatColumnGetter column,
                        SpacePointColumns requiredColumn,
                        const std::string_view& columnName) {
    return
        [column, requiredColumn, columnName](const SpacePointContainer2& self) {
          if (!self.hasColumns(requiredColumn)) {
            throw py::attribute_error(missingColumnMessage(columnName));
          }
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
          .def("createSpacePoint",
               [](py::object self) {
                 auto& c = self.cast<SpacePointContainer2&>();
                 return MutSpTether{self, c.createSpacePoint()};
               })
          .def("__len__", &SpacePointContainer2::size)
          .def("__getitem__",
               [](py::object self, SpacePointIndex2 idx) {
                 auto& c = self.cast<SpacePointContainer2&>();
                 return MutSpTether{self, MutableSpacePointProxy2(c, idx)};
               })
          .def("__iter__",
               [](py::object self) {
                 return IndexIteratorTether<SpacePointContainer2>{
                     self, 0,
                     [](const py::object& owner, SpacePointContainer2& c,
                        std::size_t i) {
                       return py::cast(MutSpTether{
                           owner, MutableSpacePointProxy2(
                                      c, static_cast<SpacePointIndex2>(i))});
                     }};
               })
          .def_property_readonly("x",
                                 floatColumn(&SpacePointContainer2::xColumn,
                                             SpacePointColumns::X, "x"))
          .def_property_readonly("y",
                                 floatColumn(&SpacePointContainer2::yColumn,
                                             SpacePointColumns::Y, "y"))
          .def_property_readonly("z",
                                 floatColumn(&SpacePointContainer2::zColumn,
                                             SpacePointColumns::Z, "z"))
          .def_property_readonly("r",
                                 floatColumn(&SpacePointContainer2::rColumn,
                                             SpacePointColumns::R, "r"))
          .def_property_readonly("phi",
                                 floatColumn(&SpacePointContainer2::phiColumn,
                                             SpacePointColumns::Phi, "phi"))
          .def_property_readonly("time",
                                 floatColumn(&SpacePointContainer2::timeColumn,
                                             SpacePointColumns::Time, "time"))
          .def_property_readonly(
              "varianceZ",
              floatColumn(&SpacePointContainer2::varianceZColumn,
                          SpacePointColumns::VarianceZ, "varianceZ"))
          .def_property_readonly(
              "varianceR",
              floatColumn(&SpacePointContainer2::varianceRColumn,
                          SpacePointColumns::VarianceR, "varianceR"))
          .def_property_readonly(
              "varianceT",
              floatColumn(&SpacePointContainer2::varianceTColumn,
                          SpacePointColumns::VarianceT, "varianceT"))
          .def_property_readonly(
              "xyColumn", arrayColumn<2>(&SpacePointContainer2::xyColumn,
                                         SpacePointColumns::PackedXY, "xy"))
          .def_property_readonly(
              "zrColumn", arrayColumn<2>(&SpacePointContainer2::zrColumn,
                                         SpacePointColumns::PackedZR, "zr"))
          .def_property_readonly(
              "xyzColumn", arrayColumn<3>(&SpacePointContainer2::xyzColumn,
                                          SpacePointColumns::PackedXYZ, "xyz"))
          .def_property_readonly(
              "xyzrColumn",
              arrayColumn<4>(&SpacePointContainer2::xyzrColumn,
                             SpacePointColumns::PackedXYZR, "xyzr"))
          .def_property_readonly(
              "varianceZRColumn",
              arrayColumn<2>(&SpacePointContainer2::varianceZRColumn,
                             SpacePointColumns::PackedVarianceZR,
                             "varianceZR"));

  WhiteBoardRegistry::registerClass(spc2);

  // ConstSpacePointProxy2
  // Note: every column is optional and only present when the corresponding
  // SpacePointColumns flag was requested. The proxy accessors guard existence
  // with assert only (no-op in release), so each getter is wrapped in
  // guardedRead to raise py::attribute_error instead of dereferencing a missing
  // column. The accessors are called on a const reference, which selects the
  // const overload unambiguously (no static_cast needed).
  py::class_<ConstSpTether>(m, "ConstSpacePointProxy2")
      .def_property_readonly("index", tetheredRead<ConstSpTether>(
                                          [](const ConstSpacePointProxy2& s) {
                                            return s.index();
                                          }))
      .def_property_readonly(
          "x", guardedRead<ConstSpTether>(
                   SpacePointColumns::X, "x",
                   [](const ConstSpacePointProxy2& s) { return s.x(); }))
      .def_property_readonly(
          "y", guardedRead<ConstSpTether>(
                   SpacePointColumns::Y, "y",
                   [](const ConstSpacePointProxy2& s) { return s.y(); }))
      .def_property_readonly(
          "z", guardedRead<ConstSpTether>(
                   SpacePointColumns::Z, "z",
                   [](const ConstSpacePointProxy2& s) { return s.z(); }))
      .def_property_readonly(
          "r", guardedRead<ConstSpTether>(
                   SpacePointColumns::R, "r",
                   [](const ConstSpacePointProxy2& s) { return s.r(); }))
      .def_property_readonly(
          "phi", guardedRead<ConstSpTether>(
                     SpacePointColumns::Phi, "phi",
                     [](const ConstSpacePointProxy2& s) { return s.phi(); }))
      .def_property_readonly(
          "time", guardedRead<ConstSpTether>(
                      SpacePointColumns::Time, "time",
                      [](const ConstSpacePointProxy2& s) { return s.time(); }))
      .def_property_readonly(
          "varianceZ",
          guardedRead<ConstSpTether>(
              SpacePointColumns::VarianceZ, "varianceZ",
              [](const ConstSpacePointProxy2& s) { return s.varianceZ(); }))
      .def_property_readonly(
          "varianceR",
          guardedRead<ConstSpTether>(
              SpacePointColumns::VarianceR, "varianceR",
              [](const ConstSpacePointProxy2& s) { return s.varianceR(); }))
      .def_property_readonly(
          "varianceT",
          guardedRead<ConstSpTether>(
              SpacePointColumns::VarianceT, "varianceT",
              [](const ConstSpacePointProxy2& s) { return s.varianceT(); }))
      .def_property_readonly(
          "sourceLinks",
          py::cpp_function(guardedRead<ConstSpTether>(
                               SpacePointColumns::SourceLinks, "sourceLinks",
                               [](const ConstSpacePointProxy2& self) {
                                 auto sls = self.sourceLinks();
                                 return py::make_iterator(sls.begin(),
                                                          sls.end());
                               }),
                           py::keep_alive<0, 1>()));

  // MutableSpacePointProxy2
  // Getters and setters are wrapped in guardedRead/guardedWrite so accessing or
  // writing a column that was not requested at construction raises
  // py::attribute_error instead of dereferencing a missing (disengaged)
  // optional column.
  using MutProxy = MutableSpacePointProxy2;
  py::class_<MutSpTether>(m, "MutableSpacePointProxy2")
      .def_property_readonly("index",
                             tetheredRead<MutSpTether>(
                                 [](const MutProxy& s) { return s.index(); }))
      .def_property(
          "x",
          guardedRead<MutSpTether>(SpacePointColumns::X, "x",
                                   [](const MutProxy& s) { return s.x(); }),
          guardedWrite<MutSpTether, float>(
              SpacePointColumns::X, "x",
              [](MutProxy& s, const float& v) { s.x() = v; }))
      .def_property(
          "y",
          guardedRead<MutSpTether>(SpacePointColumns::Y, "y",
                                   [](const MutProxy& s) { return s.y(); }),
          guardedWrite<MutSpTether, float>(
              SpacePointColumns::Y, "y",
              [](MutProxy& s, const float& v) { s.y() = v; }))
      .def_property(
          "z",
          guardedRead<MutSpTether>(SpacePointColumns::Z, "z",
                                   [](const MutProxy& s) { return s.z(); }),
          guardedWrite<MutSpTether, float>(
              SpacePointColumns::Z, "z",
              [](MutProxy& s, const float& v) { s.z() = v; }))
      .def_property(
          "r",
          guardedRead<MutSpTether>(SpacePointColumns::R, "r",
                                   [](const MutProxy& s) { return s.r(); }),
          guardedWrite<MutSpTether, float>(
              SpacePointColumns::R, "r",
              [](MutProxy& s, const float& v) { s.r() = v; }))
      .def_property(
          "phi",
          guardedRead<MutSpTether>(SpacePointColumns::Phi, "phi",
                                   [](const MutProxy& s) { return s.phi(); }),
          guardedWrite<MutSpTether, float>(
              SpacePointColumns::Phi, "phi",
              [](MutProxy& s, const float& v) { s.phi() = v; }))
      .def_property(
          "time",
          guardedRead<MutSpTether>(SpacePointColumns::Time, "time",
                                   [](const MutProxy& s) { return s.time(); }),
          guardedWrite<MutSpTether, float>(
              SpacePointColumns::Time, "time",
              [](MutProxy& s, const float& v) { s.time() = v; }))
      .def_property("varianceZ",
                    guardedRead<MutSpTether>(
                        SpacePointColumns::VarianceZ, "varianceZ",
                        [](const MutProxy& s) { return s.varianceZ(); }),
                    guardedWrite<MutSpTether, float>(
                        SpacePointColumns::VarianceZ, "varianceZ",
                        [](MutProxy& s, const float& v) { s.varianceZ() = v; }))
      .def_property("varianceR",
                    guardedRead<MutSpTether>(
                        SpacePointColumns::VarianceR, "varianceR",
                        [](const MutProxy& s) { return s.varianceR(); }),
                    guardedWrite<MutSpTether, float>(
                        SpacePointColumns::VarianceR, "varianceR",
                        [](MutProxy& s, const float& v) { s.varianceR() = v; }))
      .def_property("varianceT",
                    guardedRead<MutSpTether>(
                        SpacePointColumns::VarianceT, "varianceT",
                        [](const MutProxy& s) { return s.varianceT(); }),
                    guardedWrite<MutSpTether, float>(
                        SpacePointColumns::VarianceT, "varianceT",
                        [](MutProxy& s, const float& v) { s.varianceT() = v; }))
      .def_property("topStripVector",
                    guardedRead<MutSpTether>(
                        SpacePointColumns::TopStripVector, "topStripVector",
                        [](const MutProxy& s) { return s.topStripVector(); }),
                    guardedWrite<MutSpTether, std::array<float, 3>>(
                        SpacePointColumns::TopStripVector, "topStripVector",
                        [](MutProxy& s, const std::array<float, 3>& v) {
                          s.topStripVector() = v;
                        }))
      .def_property(
          "bottomStripVector",
          guardedRead<MutSpTether>(
              SpacePointColumns::BottomStripVector, "bottomStripVector",
              [](const MutProxy& s) { return s.bottomStripVector(); }),
          guardedWrite<MutSpTether, std::array<float, 3>>(
              SpacePointColumns::BottomStripVector, "bottomStripVector",
              [](MutProxy& s, const std::array<float, 3>& v) {
                s.bottomStripVector() = v;
              }))
      .def_property(
          "stripCenterDistance",
          guardedRead<MutSpTether>(
              SpacePointColumns::StripCenterDistance, "stripCenterDistance",
              [](const MutProxy& s) { return s.stripCenterDistance(); }),
          guardedWrite<MutSpTether, std::array<float, 3>>(
              SpacePointColumns::StripCenterDistance, "stripCenterDistance",
              [](MutProxy& s, const std::array<float, 3>& v) {
                s.stripCenterDistance() = v;
              }))
      .def_property("topStripCenter",
                    guardedRead<MutSpTether>(
                        SpacePointColumns::TopStripCenter, "topStripCenter",
                        [](const MutProxy& s) { return s.topStripCenter(); }),
                    guardedWrite<MutSpTether, std::array<float, 3>>(
                        SpacePointColumns::TopStripCenter, "topStripCenter",
                        [](MutProxy& s, const std::array<float, 3>& v) {
                          s.topStripCenter() = v;
                        }))
      .def_property("copyFromIndex",
                    guardedRead<MutSpTether>(
                        SpacePointColumns::CopyFromIndex, "copyFromIndex",
                        [](const MutProxy& s) { return s.copyFromIndex(); }),
                    guardedWrite<MutSpTether, SpacePointIndex2>(
                        SpacePointColumns::CopyFromIndex, "copyFromIndex",
                        [](MutProxy& s, const SpacePointIndex2& v) {
                          s.copyFromIndex() = v;
                        }))
      .def("assignSourceLinks",
           [](MutSpTether& self, const std::vector<SourceLink>& sourceLinks) {
             self.checked().assignSourceLinks(sourceLinks);
           })
      .def_property_readonly(
          "sourceLinks",
          py::cpp_function(guardedRead<MutSpTether>(
                               SpacePointColumns::SourceLinks, "sourceLinks",
                               [](const MutProxy& self) {
                                 auto sls = self.sourceLinks();
                                 return py::make_iterator(sls.begin(),
                                                          sls.end());
                               }),
                           py::keep_alive<0, 1>()));

  // SeedContainer2
  auto seedContainer2 =
      py::classh<SeedContainer2>(m, "SeedContainer2")
          .def(py::init<>())
          .def_property_readonly("size", &SeedContainer2::size)
          .def_property_readonly("empty", &SeedContainer2::empty)
          .def("__len__", &SeedContainer2::size)
          .def("__getitem__",
               [](py::object self, SeedIndex2 idx) {
                 const auto& c = self.cast<const SeedContainer2&>();
                 return ConstSeedTether{self, ConstSeedProxy2(c, idx)};
               })
          .def("__iter__",
               [](py::object self) {
                 return IndexIteratorTether<SeedContainer2>{
                     self, 0,
                     [](const py::object& owner, SeedContainer2& c,
                        std::size_t i) {
                       return py::cast(ConstSeedTether{
                           owner,
                           ConstSeedProxy2(c, static_cast<SeedIndex2>(i))});
                     }};
               })
          .def("createSeed",
               [](py::object self) {
                 auto& c = self.cast<SeedContainer2&>();
                 return MutSeedTether{self, c.createSeed()};
               })
          .def(
              "assignSpacePointContainer",
              [](SeedContainer2& self, const SpacePointContainer2& sp) {
                self.assignSpacePointContainer(sp);
              },
              py::keep_alive<1, 2>());

  WhiteBoardRegistry::registerClass(seedContainer2);

  // ConstSeedProxy2 (bound as a tether; seeds have no optional columns, so the
  // accessors only need the disown check).
  py::class_<ConstSeedTether>(m, "ConstSeedProxy2")
      .def_property_readonly(
          "index", tetheredRead<ConstSeedTether>(
                       [](const ConstSeedProxy2& s) { return s.index(); }))
      .def_property_readonly(
          "size", tetheredRead<ConstSeedTether>(
                      [](const ConstSeedProxy2& s) { return s.size(); }))
      .def_property_readonly(
          "empty", tetheredRead<ConstSeedTether>(
                       [](const ConstSeedProxy2& s) { return s.empty(); }))
      .def_property_readonly(
          "quality", tetheredRead<ConstSeedTether>(
                         [](const ConstSeedProxy2& s) { return s.quality(); }))
      .def_property_readonly(
          "vertexZ", tetheredRead<ConstSeedTether>(
                         [](const ConstSeedProxy2& s) { return s.vertexZ(); }))
      // spacePointIndices needs the owner (container) py::object as the numpy
      // base, so it is a custom tethered accessor rather than tetheredRead.
      .def_property_readonly(
          "spacePointIndices", [](const ConstSeedTether& self) {
            const auto& proxy = self.checked();
            return spanToNumpy1d(proxy.spacePointIndices(), self.owner);
          });

  // MutableSeedProxy2 (tether; disown check only).
  py::class_<MutSeedTether>(m, "MutableSeedProxy2")
      .def_property_readonly(
          "index", tetheredRead<MutSeedTether>(
                       [](const MutableSeedProxy2& s) { return s.index(); }))
      .def_property_readonly(
          "size", tetheredRead<MutSeedTether>(
                      [](const MutableSeedProxy2& s) { return s.size(); }))
      .def_property_readonly(
          "empty", tetheredRead<MutSeedTether>(
                       [](const MutableSeedProxy2& s) { return s.empty(); }))
      .def_property(
          "quality",
          tetheredRead<MutSeedTether>(
              [](const MutableSeedProxy2& s) { return s.quality(); }),
          tetheredWrite<MutSeedTether, float>(
              [](MutableSeedProxy2& s, const float& q) { s.quality() = q; }))
      .def_property(
          "vertexZ",
          tetheredRead<MutSeedTether>(
              [](const MutableSeedProxy2& s) { return s.vertexZ(); }),
          tetheredWrite<MutSeedTether, float>(
              [](MutableSeedProxy2& s, const float& z) { s.vertexZ() = z; }))
      .def("assignSpacePointIndices",
           [](MutSeedTether& self,
              const std::vector<SpacePointIndex2>& indices) {
             self.checked().assignSpacePointIndices(indices);
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

  // Test-only helpers: taking by unique_ptr triggers smart_holder disown,
  // replicating what the whiteboard does, so proxy tether failures can be
  // tested without acts.examples.
  auto mt = m.def_submodule("_testing");
  mt.def("consume_spacepoints", [](std::unique_ptr<SpacePointContainer2>) {});
  mt.def("consume_seeds", [](std::unique_ptr<SeedContainer2>) {});
}

}  // namespace ActsPython
