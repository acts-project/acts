// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/AnyTrackStateProxy.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "ActsExamples/EventData/Track.hpp"

#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "EventDataHelpers.hpp"

namespace py = pybind11;

using namespace ActsExamples;
using namespace ActsPython::EventDataHelpers;

namespace {

/// Setters only exist when `!Map::IsReadOnly`. isMeasurement is computed and
/// only gets a plain setIsMeasurement(), no getter/setter pair.
template <typename Map>
void addTrackStateTypeFlags(py::class_<Map>& cls) {
  cls.def_property_readonly("hasMaterial", &Map::hasMaterial)
      .def_property_readonly("hasParameters", &Map::hasParameters)
      .def_property_readonly("hasMeasurement", &Map::hasMeasurement)
      .def_property_readonly("isMeasurement", &Map::isMeasurement)
      .def_property_readonly("isOutlier", &Map::isOutlier)
      .def_property_readonly("isHole", &Map::isHole)
      .def_property_readonly("isSharedHit", &Map::isSharedHit)
      .def_property_readonly("isSplitHit", &Map::isSplitHit);

  if constexpr (!Map::IsReadOnly) {
    cls.def_property("hasMeasurement", &Map::hasMeasurement,
                     &Map::setHasMeasurement)
        .def("setIsMeasurement", &Map::setIsMeasurement)
        .def_property("isOutlier", &Map::isOutlier, &Map::setIsOutlier)
        .def_property("isHole", &Map::isHole, &Map::setIsHole)
        .def_property("isSharedHit", &Map::isSharedHit, &Map::setIsSharedHit)
        .def_property("isSplitHit", &Map::isSplitHit, &Map::setIsSplitHit);
  }
}

/// Copy a dynamic-size Eigen vector/matrix map into std::vector(s).
template <typename Map>
std::vector<double> toVector(const Map& map) {
  return std::vector<double>(map.data(), map.data() + map.size());
}
template <typename Map>
std::vector<std::vector<double>> toNestedVector(const Map& map) {
  std::vector<std::vector<double>> out(static_cast<std::size_t>(map.rows()));
  for (Eigen::Index i = 0; i < map.rows(); ++i) {
    out[static_cast<std::size_t>(i)].assign(map.row(i).begin(),
                                            map.row(i).end());
  }
  return out;
}

/// Assign a std::vector/nested std::vector into a dynamic-size Eigen map,
/// raising ValueError on a size mismatch.
template <typename Map>
void assignVector(Map map, const std::vector<double>& values) {
  if (static_cast<Eigen::Index>(values.size()) != map.size()) {
    throw std::invalid_argument("size mismatch");
  }
  std::copy(values.begin(), values.end(), map.begin());
}
template <typename Map>
void assignNestedVector(Map map,
                        const std::vector<std::vector<double>>& values) {
  if (static_cast<Eigen::Index>(values.size()) != map.rows()) {
    throw std::invalid_argument("size mismatch");
  }
  for (Eigen::Index i = 0; i < map.rows(); ++i) {
    const auto& row = values[static_cast<std::size_t>(i)];
    if (static_cast<Eigen::Index>(row.size()) != map.cols()) {
      throw std::invalid_argument("size mismatch");
    }
    std::copy(row.begin(), row.end(), map.row(i).begin());
  }
}

/// Bind the surface shared by the concrete and Any track state proxies (same
/// rationale as @ref bindTrackProxyCommon); `effectiveCalibrated`/
/// `effectiveCalibratedCovariance` stay read-only here, see below.
template <typename Cls>
void bindTrackStateProxyCommon(Cls& cls) {
  using Proxy = typename Cls::type;
  cls.def_property_readonly("index", &Proxy::index)
      .def_property_readonly("hasPrevious", &Proxy::hasPrevious)
      .def_property_readonly("referenceSurface", &Proxy::referenceSurface)
      .def_property_readonly("hasReferenceSurface", &Proxy::hasReferenceSurface)
      .def_property_readonly("mask", &Proxy::getMask)
      .def_property_readonly("hasPredicted", &Proxy::hasPredicted)
      .def_property_readonly("hasFiltered", &Proxy::hasFiltered)
      .def_property_readonly("hasSmoothed", &Proxy::hasSmoothed)
      .def_property_readonly("hasJacobian", &Proxy::hasJacobian)
      .def_property_readonly("hasProjector", &Proxy::hasProjector)
      .def_property_readonly("hasUncalibratedSourceLink",
                             &Proxy::hasUncalibratedSourceLink)
      .def_property_readonly("hasCalibrated", &Proxy::hasCalibrated)
      .def_property_readonly("uncalibratedSourceLink",
                             &Proxy::getUncalibratedSourceLink)
      .def_property_readonly("calibratedSize", &Proxy::calibratedSize)
      .def_property_readonly("projectorSubspaceIndices",
                             [](const Proxy& self) {
                               auto indices = self.projectorSubspaceIndices();
                               return std::vector<int>(
                                   indices.begin(),
                                   indices.begin() + self.calibratedSize());
                             })
      .def(
          "hasColumn",
          [](const Proxy& self, const std::string& key) {
            return self.hasColumn(Acts::hashStringDynamic(key));
          },
          py::arg("key"))
      // effectiveCalibrated(+Covariance) write access is Vector-only; see
      // bindTrackStateProxyVector.
      .def_property_readonly("effectiveCalibrated",
                             [](const Proxy& self) {
                               return toVector(self.effectiveCalibrated());
                             })
      .def_property_readonly(
          "effectiveCalibratedCovariance", [](const Proxy& self) {
            return toNestedVector(self.effectiveCalibratedCovariance());
          });

  // typeFlags returns genuinely different types depending on constness
  // (ConstTrackStateTypeMap vs. MutableTrackStateTypeMap), so it can't share
  // a getter+setter pair with the generic accessors below.
  cls.def_property_readonly(
      "typeFlags",
      py::cpp_function([](const Proxy& self) { return self.typeFlags(); },
                       py::keep_alive<0, 1>()));

  bindScalar(cls, "previous",
             [](auto& self) -> decltype(auto) { return self.previous(); });
  bindScalar(cls, "chi2",
             [](auto& self) -> decltype(auto) { return self.chi2(); });
  bindScalar(cls, "pathLength",
             [](auto& self) -> decltype(auto) { return self.pathLength(); });
  bindEigen(cls, "predicted",
            [](auto& self) -> decltype(auto) { return self.predicted(); });
  bindEigen(cls, "filtered",
            [](auto& self) -> decltype(auto) { return self.filtered(); });
  bindEigen(cls, "smoothed",
            [](auto& self) -> decltype(auto) { return self.smoothed(); });
  bindEigen(cls, "predictedCovariance", [](auto& self) -> decltype(auto) {
    return self.predictedCovariance();
  });
  bindEigen(cls, "filteredCovariance", [](auto& self) -> decltype(auto) {
    return self.filteredCovariance();
  });
  bindEigen(cls, "smoothedCovariance", [](auto& self) -> decltype(auto) {
    return self.smoothedCovariance();
  });
  bindEigen(cls, "jacobian",
            [](auto& self) -> decltype(auto) { return self.jacobian(); });
  // "best available" of predicted/filtered/smoothed; no mutable overload
  // exists at all, so always read-only regardless of Proxy.
  bindEigenReadonly(cls, "parameters",
                    [](const Proxy& self) { return self.parameters(); });
  bindEigenReadonly(cls, "covariance",
                    [](const Proxy& self) { return self.covariance(); });

  if constexpr (!Proxy::ReadOnly) {
    cls.def_property("referenceSurface", &Proxy::referenceSurface,
                     &Proxy::setReferenceSurface)
        // overrides the read-only common binding above: the mutable overload
        // returns MutableTrackStateTypeMap, not ConstTrackStateTypeMap.
        .def_property_readonly(
            "typeFlags",
            py::cpp_function([](Proxy& self) { return self.typeFlags(); },
                             py::keep_alive<0, 1>()))
        .def_property(
            "uncalibratedSourceLink", &Proxy::getUncalibratedSourceLink,
            [](Proxy& self, const Acts::SourceLink& sourceLink) {
              // Some concrete setUncalibratedSourceLink overloads take
              // SourceLink&&; pybind's argument caster only ever produces an
              // lvalue, so a direct member-pointer binding doesn't compile.
              self.setUncalibratedSourceLink(Acts::SourceLink{sourceLink});
            })
        .def(
            "setProjectorSubspaceIndices",
            [](Proxy& self, const std::vector<int>& indices) {
              self.setProjectorSubspaceIndices(indices);
            },
            py::arg("indices"))
        .def(
            "allocateCalibrated",
            [](Proxy& self, std::size_t measdim) {
              self.allocateCalibrated(measdim);
            },
            py::arg("measdim"))
        .def("unset", &Proxy::unset, py::arg("mask"))
        .def("addComponents", &Proxy::addComponents, py::arg("mask"));
  }
}

/// Bind what only concrete Vector-backed proxies have: `copyFrom`, and write
/// access to `effectiveCalibrated`/`effectiveCalibratedCovariance` (needs a
/// contiguous `Eigen::Map`, not guaranteed by `AnyTrackStateProxy`).
template <typename Cls>
void bindTrackStateProxyVector(Cls& cls) {
  using Proxy = typename Cls::type;
  if constexpr (!Proxy::ReadOnly) {
    cls.def("copyFrom", &Proxy::template copyFrom<ConstTrackStateProxy>,
            py::arg("other"), py::arg("mask") = Acts::TrackStatePropMask::All,
            py::arg("onlyAllocated") = true)
        .def("copyFrom", &Proxy::template copyFrom<TrackStateProxy>,
             py::arg("other"), py::arg("mask") = Acts::TrackStatePropMask::All,
             py::arg("onlyAllocated") = true)
        .def_property(
            "effectiveCalibrated",
            [](const Proxy& self) {
              return toVector(self.effectiveCalibrated());
            },
            [](Proxy& self, const std::vector<double>& v) {
              assignVector(self.effectiveCalibrated(), v);
            })
        .def_property(
            "effectiveCalibratedCovariance",
            [](const Proxy& self) {
              return toNestedVector(self.effectiveCalibratedCovariance());
            },
            [](Proxy& self, const std::vector<std::vector<double>>& v) {
              assignNestedVector(self.effectiveCalibratedCovariance(), v);
            });
  }
}

}  // namespace

namespace ActsPython {

void addEventDataTrackStates(py::module& mex) {
  py::class_<Acts::TrackStateType> trackStateTypeFlags(mex,
                                                       "TrackStateTypeFlags");
  addTrackStateTypeFlags(trackStateTypeFlags);

  py::class_<Acts::ConstTrackStateTypeMap> constTrackStateTypeFlags(
      mex, "ConstTrackStateTypeFlags");
  addTrackStateTypeFlags(constTrackStateTypeFlags);

  py::class_<Acts::MutableTrackStateTypeMap> mutableTrackStateTypeFlags(
      mex, "MutableTrackStateTypeFlags");
  addTrackStateTypeFlags(mutableTrackStateTypeFlags);

  py::class_<ConstTrackStateProxy> constTrackState(mex, "ConstTrackStateProxy");
  bindTrackStateProxyCommon(constTrackState);
  bindTrackStateProxyVector(constTrackState);

  py::class_<TrackStateProxy> trackState(mex, "TrackStateProxy");
  bindTrackStateProxyCommon(trackState);
  bindTrackStateProxyVector(trackState);

  py::classh<Acts::AnyConstTrackStateProxy> anyConstTrackState(
      mex, "AnyConstTrackStateProxy");
  bindTrackStateProxyCommon(anyConstTrackState);
  anyConstTrackState.def(py::init<ConstTrackStateProxy&>());
  py::implicitly_convertible<ConstTrackStateProxy,
                             Acts::AnyConstTrackStateProxy>();
  // Kept for backwards compatibility with the previous (misspelled) name.
  mex.attr("AnyConstTrackState") = mex.attr("AnyConstTrackStateProxy");

  py::classh<Acts::AnyMutableTrackStateProxy> anyMutableTrackState(
      mex, "AnyMutableTrackStateProxy");
  bindTrackStateProxyCommon(anyMutableTrackState);
  anyMutableTrackState.def(py::init<TrackStateProxy&>());
  py::implicitly_convertible<TrackStateProxy,
                             Acts::AnyMutableTrackStateProxy>();
}

}  // namespace ActsPython
