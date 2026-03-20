// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/SpacePointContainer2.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsPython/Utilities/WhiteBoardRegistry.hpp"

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

// Prevent stl.h's list_caster-based type_caster<std::vector<T>> from matching
// ProtoTrackContainer, which would break py::cast<std::unique_ptr<T>> needed
// by WhiteBoardRegistry. The full specialization takes priority over stl.h's
// partial specialization regardless of include order.
PYBIND11_MAKE_OPAQUE(ActsExamples::ProtoTrackContainer)

namespace py = pybind11;

using namespace ActsExamples;

namespace ActsPython {

void addEventData(py::module& mex) {
  py::class_<Acts::TrackStateType>(mex, "TrackStateTypeFlags")
      .def_property_readonly("isMeasurement",
                             &Acts::TrackStateType::isMeasurement)
      .def_property_readonly("isOutlier", &Acts::TrackStateType::isOutlier)
      .def_property_readonly("isHole", &Acts::TrackStateType::isHole)
      .def_property_readonly("hasMaterial", &Acts::TrackStateType::hasMaterial)
      .def_property_readonly("isSharedHit", &Acts::TrackStateType::isSharedHit)
      .def_property_readonly("hasParameters",
                             &Acts::TrackStateType::hasParameters);

  py::class_<ConstTrackStateProxy>(mex, "ConstTrackStateProxy")
      .def_property_readonly(
          "typeFlags",
          [](const ConstTrackStateProxy& self) {
            return Acts::TrackStateType{self.typeFlags().raw()};
          })
      .def_property_readonly("hasPredicted",
                             &ConstTrackStateProxy::hasPredicted)
      .def_property_readonly("hasFiltered", &ConstTrackStateProxy::hasFiltered)
      .def_property_readonly("hasSmoothed", &ConstTrackStateProxy::hasSmoothed)
      .def_property_readonly("predicted",
                             [](const ConstTrackStateProxy& self) {
                               return Acts::BoundVector{self.predicted()};
                             })
      .def_property_readonly("filtered",
                             [](const ConstTrackStateProxy& self) {
                               return Acts::BoundVector{self.filtered()};
                             })
      .def_property_readonly("smoothed",
                             [](const ConstTrackStateProxy& self) {
                               return Acts::BoundVector{self.smoothed()};
                             })
      .def_property_readonly(
          "predictedCovariance",
          [](const ConstTrackStateProxy& self) {
            return Acts::BoundMatrix{self.predictedCovariance()};
          })
      .def_property_readonly(
          "filteredCovariance",
          [](const ConstTrackStateProxy& self) {
            return Acts::BoundMatrix{self.filteredCovariance()};
          })
      .def_property_readonly(
          "smoothedCovariance",
          [](const ConstTrackStateProxy& self) {
            return Acts::BoundMatrix{self.smoothedCovariance()};
          })
      .def_property_readonly("pathLength", &ConstTrackStateProxy::pathLength);

  auto constTrackProxy =
      py::class_<ConstTrackProxy>(mex, "ConstTrackProxy")
          .def_property_readonly("index", &ConstTrackProxy::index)
          .def_property_readonly(
              "tipIndex",
              [](const ConstTrackProxy& self) { return self.tipIndex(); })
          .def_property_readonly(
              "stemIndex",
              [](const ConstTrackProxy& self) { return self.stemIndex(); })
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
                                 &ConstTrackProxy::particleHypothesis)
          .def_property_readonly(
              "nMeasurements",
              [](const ConstTrackProxy& self) { return self.nMeasurements(); })
          .def_property_readonly(
              "nHoles",
              [](const ConstTrackProxy& self) { return self.nHoles(); })
          .def_property_readonly("isForwardLinked",
                                 [](const ConstTrackProxy& self) {
                                   return self.isForwardLinked();
                                 })
          .def_property_readonly("trackStatesReversed",
                                 py::cpp_function(
                                     [](const ConstTrackProxy& self) {
                                       auto range = self.trackStatesReversed();
                                       return py::make_iterator(range.begin(),
                                                                range.end());
                                     },
                                     py::keep_alive<0, 1>()))
          .def_property_readonly("trackStates",
                                 py::cpp_function(
                                     [](const ConstTrackProxy& self) {
                                       auto range = self.trackStates();
                                       return py::make_iterator(range.begin(),
                                                                range.end());
                                     },
                                     py::keep_alive<0, 1>()));

  // Mark a numpy array as non-writeable and return it.
  const auto readOnly = [](auto arr) {
    arr.attr("flags").attr("writeable") = py::bool_(false);
    return arr;
  };

  // Factory for zero-copy 1-D views over a plain SoA column.
  // `accessor` is called with the backend and must return a const-ref to the
  // desired std::vector member.  dtype and stride are derived automatically.
  const auto col1D = [&readOnly](auto accessor) {
    return [accessor, &readOnly](const py::object& self_py) {
      const auto& backend =
          self_py.cast<const ConstTrackContainer&>().container();
      const auto N = static_cast<py::ssize_t>(backend.size_impl());
      const auto& vec = accessor(backend);
      using T = typename std::decay_t<decltype(vec)>::value_type;
      if (N == 0) {
        return readOnly(py::array_t<T>(py::ssize_t{0}));
      }
      return readOnly(py::array_t<T>({N}, {static_cast<py::ssize_t>(sizeof(T))},
                                     vec.data(), self_py));
    };
  };

  auto constTrackContainer =
      py::classh<ConstTrackContainer>(mex, "ConstTrackContainer")
          .def("__len__", &ConstTrackContainer::size)
          .def("__iter__",
               [](const ConstTrackContainer& self) {
                 return py::make_iterator(self.begin(), self.end());
               })
          .def("__getitem__", py::overload_cast<ConstTrackContainer::IndexType>(
                                  &ConstTrackContainer::getTrack, py::const_))

          // Zero-copy numpy array views of the underlying SoA columns.
          // The returned arrays are read-only and keep the container alive via
          // the numpy base-object mechanism as long as any array is alive.

          // shape (N, 6), float64 — bound track parameters [loc0, loc1, phi,
          // theta, q/p, t]. Strides account for potential Eigen alignment
          // padding between entries.
          .def_property_readonly(
              "parameters",
              [&readOnly](const py::object& self_py) -> py::array_t<double> {
                using CoeffsType = Acts::detail_tsp::FixedSizeTypes<
                    Acts::eBoundSize>::Coefficients;
                const auto& backend =
                    self_py.cast<const ConstTrackContainer&>().container();
                const auto N = static_cast<py::ssize_t>(backend.size_impl());
                if (N == 0) {
                  return readOnly(py::array_t<double>({N, py::ssize_t{6}}));
                }
                return readOnly(py::array_t<double>(
                    {N, py::ssize_t{6}},
                    {static_cast<py::ssize_t>(sizeof(CoeffsType)),
                     static_cast<py::ssize_t>(sizeof(double))},
                    backend.m_params[0].data(), self_py));
              })

          // shape (N, 6, 6), float64, column-major sub-matrices (Eigen
          // default). arr[k, i, j] is row i, column j of track k's covariance.
          .def_property_readonly(
              "covariance",
              [&readOnly](const py::object& self_py) -> py::array_t<double> {
                using CovType = Acts::detail_tsp::FixedSizeTypes<
                    Acts::eBoundSize>::Covariance;
                const auto& backend =
                    self_py.cast<const ConstTrackContainer&>().container();
                const auto N = static_cast<py::ssize_t>(backend.size_impl());
                if (N == 0) {
                  return readOnly(
                      py::array_t<double>({N, py::ssize_t{6}, py::ssize_t{6}}));
                }
                // Eigen column-major: row stride = 1 double, col stride = 6
                // doubles.
                constexpr py::ssize_t dbl = sizeof(double);
                return readOnly(py::array_t<double>(
                    {N, py::ssize_t{6}, py::ssize_t{6}},
                    {static_cast<py::ssize_t>(sizeof(CovType)), dbl, 6 * dbl},
                    backend.m_cov[0].data(), self_py));
              })

          .def_property_readonly(
              "tipIndex",
              col1D([](const auto& b) -> const auto& { return b.m_tipIndex; }))
          .def_property_readonly(
              "stemIndex",
              col1D([](const auto& b) -> const auto& { return b.m_stemIndex; }))
          .def_property_readonly("nMeasurements",
                                 col1D([](const auto& b) -> const auto& {
                                   return b.m_nMeasurements;
                                 }))
          .def_property_readonly(
              "nHoles",
              col1D([](const auto& b) -> const auto& { return b.m_nHoles; }))
          .def_property_readonly(
              "chi2",
              col1D([](const auto& b) -> const auto& { return b.m_chi2; }))
          .def_property_readonly("ndf", col1D([](const auto& b) -> const auto& {
                                   return b.m_ndf;
                                 }));

  WhiteBoardRegistry::registerClass(constTrackContainer);

  py::bind_vector<ProtoTrack>(mex, "ProtoTrack");

  auto protoTrackContainer =
      py::bind_vector<ProtoTrackContainer, py::smart_holder>(
          mex, "ProtoTrackContainer");
  WhiteBoardRegistry::registerClass(protoTrackContainer);

  mex.attr("kTrackIndexInvalid") = Acts::kTrackIndexInvalid;
}

}  // namespace ActsPython
