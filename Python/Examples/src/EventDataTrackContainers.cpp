// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/Track.hpp"
#include "ActsPython/Utilities/WhiteBoardRegistry.hpp"

#include <memory>
#include <string>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;

namespace ActsPython {

void addEventDataTrackContainers(py::module& mex) {
  // Mark a numpy array as non-writeable and return it.
  const auto readOnly = [](auto arr) {
    arr.attr("flags").attr("writeable") = py::bool_(false);
    return arr;
  };

  // Factory for zero-copy 1-D views over a plain SoA column.
  // `accessor` is called with the backend and must return a const-ref to the
  // desired std::vector member.  dtype and stride are derived automatically.
  const auto col1D = [readOnly](auto accessor) {
    return [accessor, readOnly](const py::object& self_py) {
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
          .def(
              "__iter__",
              [](const ConstTrackContainer& self) {
                return py::make_iterator(self.begin(), self.end());
              },
              py::keep_alive<0, 1>())
          .def("__getitem__",
               py::overload_cast<ConstTrackContainer::IndexType>(
                   &ConstTrackContainer::getTrack, py::const_),
               py::keep_alive<0, 1>())
          .def("getTrack",
               py::overload_cast<ConstTrackContainer::IndexType>(
                   &ConstTrackContainer::getTrack, py::const_),
               py::arg("index"), py::keep_alive<0, 1>())
          .def(
              "hasColumn",
              [](const ConstTrackContainer& self, const std::string& key) {
                return self.hasColumn(key);
              },
              py::arg("key"))
          // Build an independent, fully mutable copy; the source is left
          // untouched (unlike TrackContainer::makeConst, which moves out of
          // its backends) since a ConstTrackContainer may be shared.
          .def("makeMutable",
               [](const ConstTrackContainer& self) {
                 return TrackContainer{
                     std::make_shared<Acts::VectorTrackContainer>(
                         self.container()),
                     std::make_shared<Acts::VectorMultiTrajectory>(
                         self.trackStateContainer())};
               })

          // Zero-copy numpy array views of the underlying SoA columns.
          // The returned arrays are read-only and keep the container alive via
          // the numpy base-object mechanism as long as any array is alive.

          // shape (N, 6), float64 — bound track parameters [loc0, loc1, phi,
          // theta, q/p, t]. Strides account for potential Eigen alignment
          // padding between entries.
          .def_property_readonly(
              "parameters",
              [readOnly](const py::object& self_py) -> py::array_t<double> {
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
              [readOnly](const py::object& self_py) -> py::array_t<double> {
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
                                 }))
          .def_property_readonly(
              "nOutliers",
              col1D([](const auto& b) -> const auto& { return b.m_nOutliers; }))
          .def_property_readonly("nSharedHits",
                                 col1D([](const auto& b) -> const auto& {
                                   return b.m_nSharedHits;
                                 }))
          .def_property_readonly("nSplitHits",
                                 col1D([](const auto& b) -> const auto& {
                                   return b.m_nSplitHits;
                                 }));

  WhiteBoardRegistry::registerClass(constTrackContainer);

  py::class_<TrackContainer>(mex, "TrackContainer")
      .def(py::init([]() {
        return TrackContainer{std::make_shared<Acts::VectorTrackContainer>(),
                              std::make_shared<Acts::VectorMultiTrajectory>()};
      }))
      .def(
          py::init([](const ConstTrackContainer& other) {
            return TrackContainer{
                std::make_shared<Acts::VectorTrackContainer>(other.container()),
                std::make_shared<Acts::VectorMultiTrajectory>(
                    other.trackStateContainer())};
          }),
          py::arg("other"))
      .def("__len__", &TrackContainer::size)
      .def(
          "__iter__",
          [](TrackContainer& self) {
            return py::make_iterator(self.begin(), self.end());
          },
          py::keep_alive<0, 1>())
      .def("__getitem__",
           py::overload_cast<TrackContainer::IndexType>(
               &TrackContainer::getTrack),
           py::keep_alive<0, 1>())
      .def("getTrack",
           py::overload_cast<TrackContainer::IndexType>(
               &TrackContainer::getTrack),
           py::arg("index"), py::keep_alive<0, 1>())
      .def("makeTrack", &TrackContainer::makeTrack, py::keep_alive<0, 1>())
      .def("removeTrack", &TrackContainer::removeTrack, py::arg("index"))
      .def("clear", &TrackContainer::clear)
      .def(
          "hasColumn",
          [](const TrackContainer& self, const std::string& key) {
            return self.hasColumn(key);
          },
          py::arg("key"))
      .def(
          "ensureDynamicColumns",
          [](TrackContainer& self, const ConstTrackContainer& other) {
            self.ensureDynamicColumns(other);
          },
          py::arg("other"))
      .def(
          "ensureDynamicColumns",
          [](TrackContainer& self, const TrackContainer& other) {
            self.ensureDynamicColumns(other);
          },
          py::arg("other"))
      .def("makeConst", [](TrackContainer& self) {
        return ConstTrackContainer{
            std::make_shared<Acts::ConstVectorTrackContainer>(
                std::move(self.container())),
            std::make_shared<Acts::ConstVectorMultiTrajectory>(
                std::move(self.trackStateContainer()))};
      });
}

}  // namespace ActsPython
