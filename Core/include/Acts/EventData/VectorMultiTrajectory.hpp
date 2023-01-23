// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/detail/DynamicColumn.hpp"

#include <unordered_map>

#include <boost/histogram.hpp>

namespace Acts {

namespace detail_vmt {

using MultiTrajectoryTraits::IndexType;
constexpr auto kInvalid = MultiTrajectoryTraits::kInvalid;
constexpr auto MeasurementSizeMax = MultiTrajectoryTraits::MeasurementSizeMax;

class VectorMultiTrajectoryBase {
 public:
  struct Statistics {
    using axis_t = boost::histogram::axis::variant<
        boost::histogram::axis::category<std::string>,
        boost::histogram::axis::category<>>;

    using axes_t = std::vector<axis_t>;
    using hist_t = boost::histogram::histogram<axes_t>;

    hist_t hist;

    void toStream(std::ostream& os, size_t n = 1);
  };

  template <typename T>
  Statistics statistics(T& instance) const {
    using namespace boost::histogram;
    using cat = axis::category<std::string>;

    Statistics::axes_t axes;
    axes.emplace_back(cat({
        "count",
        "index",
        "parPred",
        "covPred",
        "parFilt",
        "covFilt",
        "parSmth",
        "covSmth",
        "meas",
        "measOffset",
        "measCov",
        "measCovOffset",
        "jac",
        "sourceLinks",
        "projectors",
    }));

    axes.emplace_back(axis::category<>({0, 1}));

    auto h = make_histogram(axes);

    for (IndexType i = 0; i < instance.size(); i++) {
      auto ts = instance.getTrackState(i);

      bool isMeas = ts.typeFlags().test(TrackStateFlag::MeasurementFlag);

      h("count", isMeas);

      h("index", isMeas, weight(sizeof(IndexData)));

      using scalar = typename decltype(ts.predicted())::Scalar;
      size_t par_size = eBoundSize * sizeof(scalar);
      size_t cov_size = eBoundSize * eBoundSize * sizeof(scalar);

      const IndexData& index = m_index[i];
      if (ts.hasPredicted() &&
          ACTS_CHECK_BIT(index.allocMask, TrackStatePropMask::Predicted)) {
        h("parPred", isMeas, weight(par_size));
        h("covPred", isMeas, weight(cov_size));
      }
      if (ts.hasFiltered() &&
          ACTS_CHECK_BIT(index.allocMask, TrackStatePropMask::Filtered)) {
        h("parFilt", isMeas, weight(par_size));
        h("covFilt", isMeas, weight(cov_size));
      }
      if (ts.hasSmoothed() &&
          ACTS_CHECK_BIT(index.allocMask, TrackStatePropMask::Smoothed)) {
        h("parSmth", isMeas, weight(par_size));
        h("covSmth", isMeas, weight(cov_size));
      }
      h("sourceLinks", isMeas, weight(sizeof(SourceLink)));
      h("measOffset", isMeas,
        weight(sizeof(decltype(m_measOffset)::value_type)));
      h("measCovOffset", isMeas,
        weight(sizeof(decltype(m_measCovOffset)::value_type)));
      if (ts.hasCalibrated() &&
          ACTS_CHECK_BIT(index.allocMask, TrackStatePropMask::Calibrated)) {
        size_t meas_size = ts.calibratedSize() * sizeof(scalar);
        size_t meas_cov_size =
            ts.calibratedSize() * ts.calibratedSize() * sizeof(scalar);

        h("meas", isMeas, weight(meas_size));
        h("measCov", isMeas, weight(meas_cov_size));
        h("sourceLinks", isMeas, weight(sizeof(const SourceLink)));
        h("projectors", isMeas, weight(sizeof(ProjectorBitset)));
      }

      if (ts.hasJacobian() &&
          ACTS_CHECK_BIT(index.allocMask, TrackStatePropMask::Jacobian)) {
        h("jac", isMeas, weight(cov_size));
      }
    }

    return Statistics{h};
  }

 protected:
  struct IndexData {
    IndexType iprevious = kInvalid;
    IndexType ipredicted = kInvalid;
    IndexType ifiltered = kInvalid;
    IndexType ismoothed = kInvalid;
    IndexType ijacobian = kInvalid;
    IndexType iprojector = kInvalid;

    double chi2 = 0;
    double pathLength = 0;
    TrackStateType typeFlags;

    IndexType iuncalibrated = kInvalid;
    IndexType icalibratedsourcelink = kInvalid;
    IndexType measdim = 0;

    TrackStatePropMask allocMask = TrackStatePropMask::None;
  };

  VectorMultiTrajectoryBase() = default;

  VectorMultiTrajectoryBase(const VectorMultiTrajectoryBase& other)
      : m_index{other.m_index},
        m_previous{other.m_previous},
        m_params{other.m_params},
        m_cov{other.m_cov},
        m_meas{other.m_meas},
        m_measCov{other.m_measCov},
        m_jac{other.m_jac},
        m_sourceLinks{other.m_sourceLinks},
        m_projectors{other.m_projectors},
        m_referenceSurfaces{other.m_referenceSurfaces} {
    for (const auto& [key, value] : other.m_dynamic) {
      m_dynamic.insert({key, value->clone()});
    }
  };

  VectorMultiTrajectoryBase(VectorMultiTrajectoryBase&& other) = default;

  // BEGIN INTERFACE HELPER
  template <typename T>
  static constexpr bool has_impl(T& instance, HashedString key,
                                 IndexType istate) {
    using namespace Acts::HashedStringLiteral;
    switch (key) {
      case "predicted"_hash:
        return instance.m_index[istate].ipredicted != kInvalid;
      case "filtered"_hash:
        return instance.m_index[istate].ifiltered != kInvalid;
      case "smoothed"_hash:
        return instance.m_index[istate].ismoothed != kInvalid;
      case "calibrated"_hash:
        return instance.m_measOffset[istate] != kInvalid;
      case "calibratedCov"_hash:
        return instance.m_measCovOffset[istate] != kInvalid;
      case "jacobian"_hash:
        return instance.m_index[istate].ijacobian != kInvalid;
      case "projector"_hash:
        return instance.m_index[istate].iprojector != kInvalid;
      case "uncalibratedSourceLink"_hash:
        return instance.m_sourceLinks[instance.m_index[istate].iuncalibrated]
            .has_value();
      case "previous"_hash:
      case "calibratedSourceLink"_hash:
      case "referenceSurface"_hash:
      case "measdim"_hash:
      case "chi2"_hash:
      case "pathLength"_hash:
      case "typeFlags"_hash:
        return true;
      default:
        return instance.m_dynamic.find(key) != instance.m_dynamic.end();
    }
  }

  template <bool EnsureConst, typename T>
  static std::any component_impl(T& instance, HashedString key,
                                 IndexType istate) {
    if constexpr (EnsureConst) {
      static_assert(std::is_const_v<std::remove_reference_t<T>>,
                    "Is not const");
    }
    using namespace Acts::HashedStringLiteral;
    switch (key) {
      case "previous"_hash:
        return &instance.m_index[istate].iprevious;
      case "predicted"_hash:
        return &instance.m_index[istate].ipredicted;
      case "filtered"_hash:
        return &instance.m_index[istate].ifiltered;
      case "smoothed"_hash:
        return &instance.m_index[istate].ismoothed;
      case "calibrated"_hash:
        return &instance.m_measOffset[istate];
      case "calibratedCov"_hash:
        return &instance.m_measCovOffset[istate];
      case "jacobian"_hash:
        return &instance.m_index[istate].ijacobian;
      case "projector"_hash:
        return &instance.m_projectors[instance.m_index[istate].iprojector];
      case "uncalibratedSourceLink"_hash:
        return &instance.m_sourceLinks[instance.m_index[istate].iuncalibrated];
      case "calibratedSourceLink"_hash:
        return &instance.m_sourceLinks[instance.m_index[istate]
                                           .icalibratedsourcelink];
      case "referenceSurface"_hash:
        return &instance.m_referenceSurfaces[istate];
      case "measdim"_hash:
        return &instance.m_index[istate].measdim;
      case "chi2"_hash:
        return &instance.m_index[istate].chi2;
      case "pathLength"_hash:
        return &instance.m_index[istate].pathLength;
      case "typeFlags"_hash:
        return &instance.m_index[istate].typeFlags;
      default:
        auto it = instance.m_dynamic.find(key);
        if (it == instance.m_dynamic.end()) {
          throw std::runtime_error("Unable to handle this component");
        }
        std::conditional_t<EnsureConst, const detail::DynamicColumnBase*,
                           detail::DynamicColumnBase*>
            col = it->second.get();
        assert(col && "Dynamic column is null");
        return col->get(istate);
    }
  }

  template <typename T>
  static constexpr bool hasColumn_impl(T& instance, HashedString key) {
    using namespace Acts::HashedStringLiteral;
    switch (key) {
      case "predicted"_hash:
      case "filtered"_hash:
      case "smoothed"_hash:
      case "calibrated"_hash:
      case "calibratedCov"_hash:
      case "jacobian"_hash:
      case "projector"_hash:
      case "previous"_hash:
      case "uncalibratedSourceLink"_hash:
      case "calibratedSourceLink"_hash:
      case "referenceSurface"_hash:
      case "measdim"_hash:
      case "chi2"_hash:
      case "pathLength"_hash:
      case "typeFlags"_hash:
        return true;
      default:
        return instance.m_dynamic.find(key) != instance.m_dynamic.end();
    }
  }

  IndexType calibratedSize_impl(IndexType istate) const {
    return m_index[istate].measdim;
  }

  // END INTERFACE HELPER

  /// index to map track states to the corresponding
  std::vector<IndexData> m_index;
  std::vector<IndexType> m_previous;
  std::vector<typename detail_lt::Types<eBoundSize>::Coefficients> m_params;
  std::vector<typename detail_lt::Types<eBoundSize>::Covariance> m_cov;

  std::vector<double> m_meas;
  std::vector<MultiTrajectoryTraits::IndexType> m_measOffset;
  std::vector<double> m_measCov;
  std::vector<MultiTrajectoryTraits::IndexType> m_measCovOffset;

  std::vector<typename detail_lt::Types<eBoundSize>::Covariance> m_jac;
  std::vector<std::optional<SourceLink>> m_sourceLinks;
  std::vector<ProjectorBitset> m_projectors;

  // owning vector of shared pointers to surfaces
  //
  // This might be problematic when appending a large number of surfaces
  // trackstates, because vector has to reallocated and thus copy. This might
  // be handled in a smart way by moving but not sure.
  std::vector<std::shared_ptr<const Surface>> m_referenceSurfaces;

  std::unordered_map<HashedString, std::unique_ptr<detail::DynamicColumnBase>>
      m_dynamic;
};

}  // namespace detail_vmt

class VectorMultiTrajectory;
template <>
struct IsReadOnlyMultiTrajectory<VectorMultiTrajectory> : std::false_type {};

class VectorMultiTrajectory final
    : public detail_vmt::VectorMultiTrajectoryBase,
      public MultiTrajectory<VectorMultiTrajectory> {
#ifndef DOXYGEN
  friend MultiTrajectory<VectorMultiTrajectory>;
#endif

 public:
  VectorMultiTrajectory() = default;
  VectorMultiTrajectory(const VectorMultiTrajectory& other)
      : VectorMultiTrajectoryBase{other} {}

  VectorMultiTrajectory(VectorMultiTrajectory&& other)
      : VectorMultiTrajectoryBase{std::move(other)} {}

  Statistics statistics() const {
    return detail_vmt::VectorMultiTrajectoryBase::statistics(*this);
  }

 private:
  // BEGIN INTERFACE
  TrackStateProxy::Parameters parameters_impl(IndexType parIdx) {
    return TrackStateProxy::Parameters{m_params[parIdx].data()};
  }

  ConstTrackStateProxy::Parameters parameters_impl(IndexType parIdx) const {
    return ConstTrackStateProxy::Parameters{m_params[parIdx].data()};
  }

  TrackStateProxy::Covariance covariance_impl(IndexType parIdx) {
    return TrackStateProxy::Covariance{m_cov[parIdx].data()};
  }

  ConstTrackStateProxy::Covariance covariance_impl(IndexType parIdx) const {
    return ConstTrackStateProxy::Covariance{m_cov[parIdx].data()};
  }

  TrackStateProxy::Covariance jacobian_impl(IndexType parIdx) {
    return TrackStateProxy::Covariance{m_jac[parIdx].data()};
  }

  ConstTrackStateProxy::Covariance jacobian_impl(IndexType parIdx) const {
    return ConstTrackStateProxy::Covariance{m_jac[parIdx].data()};
  }

  template <size_t measdim>
  TrackStateProxy::Measurement<measdim> measurement_impl(IndexType offset) {
    return TrackStateProxy::Measurement<measdim>{&m_meas[offset]};
  }

  template <size_t measdim>
  ConstTrackStateProxy::Measurement<measdim> measurement_impl(
      IndexType offset) const {
    return ConstTrackStateProxy::Measurement<measdim>{&m_meas[offset]};
  }

  template <size_t measdim>
  TrackStateProxy::MeasurementCovariance<measdim> measurementCovariance_impl(
      IndexType offset) {
    return TrackStateProxy::MeasurementCovariance<measdim>{&m_measCov[offset]};
  }

  template <size_t measdim>
  ConstTrackStateProxy::MeasurementCovariance<measdim>
  measurementCovariance_impl(IndexType offset) const {
    return ConstTrackStateProxy::MeasurementCovariance<measdim>{
        &m_measCov[offset]};
  }

  IndexType addTrackState_impl(
      TrackStatePropMask mask = TrackStatePropMask::All,
      IndexType iprevious = kInvalid);

  void shareFrom_impl(IndexType iself, IndexType iother,
                      TrackStatePropMask shareSource,
                      TrackStatePropMask shareTarget);

  void unset_impl(TrackStatePropMask target, IndexType istate);

  constexpr bool has_impl(HashedString key, IndexType istate) const {
    return detail_vmt::VectorMultiTrajectoryBase::has_impl(*this, key, istate);
  }

  IndexType size_impl() const { return m_index.size(); }

  void clear_impl();

  std::any component_impl(HashedString key, IndexType istate) {
    return detail_vmt::VectorMultiTrajectoryBase::component_impl<false>(
        *this, key, istate);
  }

  std::any component_impl(HashedString key, IndexType istate) const {
    return detail_vmt::VectorMultiTrajectoryBase::component_impl<true>(
        *this, key, istate);
  }

  template <typename T>
  constexpr void addColumn_impl(const std::string& key) {
    m_dynamic.insert(
        {hashString(key), std::make_unique<detail::DynamicColumn<T>>()});
  }

  constexpr bool hasColumn_impl(HashedString key) const {
    return detail_vmt::VectorMultiTrajectoryBase::hasColumn_impl(*this, key);
  }

  void allocateCalibrated_impl(IndexType istate, size_t measdim) {
    if (m_measOffset[istate] != kInvalid &&
        m_measCovOffset[istate] != kInvalid &&
        m_index[istate].measdim == measdim) {
      return;
    }

    m_index[istate].measdim = measdim;

    m_measOffset[istate] = static_cast<IndexType>(m_meas.size());
    m_meas.resize(m_meas.size() + measdim);

    m_measCovOffset[istate] = static_cast<IndexType>(m_measCov.size());
    m_measCov.resize(m_measCov.size() + measdim * measdim);
  }

  // END INTERFACE
};

class ConstVectorMultiTrajectory;
template <>
struct IsReadOnlyMultiTrajectory<ConstVectorMultiTrajectory> : std::true_type {
};

class ConstVectorMultiTrajectory final
    : public detail_vmt::VectorMultiTrajectoryBase,
      public MultiTrajectory<ConstVectorMultiTrajectory> {
#ifndef DOXYGEN
  friend MultiTrajectory<ConstVectorMultiTrajectory>;
#endif

 public:
  ConstVectorMultiTrajectory() = default;

  ConstVectorMultiTrajectory(const ConstVectorMultiTrajectory& other)
      : VectorMultiTrajectoryBase{other} {}

  ConstVectorMultiTrajectory(const VectorMultiTrajectory& other)
      : VectorMultiTrajectoryBase{other} {}

  ConstVectorMultiTrajectory(VectorMultiTrajectory&& other)
      : VectorMultiTrajectoryBase{std::move(other)} {}

  ConstVectorMultiTrajectory(ConstVectorMultiTrajectory&&) = default;

  Statistics statistics() const {
    return detail_vmt::VectorMultiTrajectoryBase::statistics(*this);
  }

 private:
  // BEGIN INTERFACE

  ConstTrackStateProxy::Parameters parameters_impl(IndexType parIdx) const {
    return ConstTrackStateProxy::Parameters{m_params[parIdx].data()};
  }

  ConstTrackStateProxy::Covariance covariance_impl(IndexType parIdx) const {
    return ConstTrackStateProxy::Covariance{m_cov[parIdx].data()};
  }

  ConstTrackStateProxy::Covariance jacobian_impl(IndexType parIdx) const {
    return ConstTrackStateProxy::Covariance{m_jac[parIdx].data()};
  }

  template <size_t measdim>
  ConstTrackStateProxy::Measurement<measdim> measurement_impl(
      IndexType offset) const {
    return ConstTrackStateProxy::Measurement<measdim>{&m_meas[offset]};
  }

  template <size_t measdim>
  ConstTrackStateProxy::MeasurementCovariance<measdim>
  measurementCovariance_impl(IndexType offset) const {
    return ConstTrackStateProxy::MeasurementCovariance<measdim>{
        &m_measCov[offset]};
  }

  constexpr bool has_impl(HashedString key, IndexType istate) const {
    return detail_vmt::VectorMultiTrajectoryBase::has_impl(*this, key, istate);
  }

  IndexType size_impl() const { return m_index.size(); }

  std::any component_impl(HashedString key, IndexType istate) const {
    return detail_vmt::VectorMultiTrajectoryBase::component_impl<true>(
        *this, key, istate);
  }

  constexpr bool hasColumn_impl(HashedString key) const {
    return detail_vmt::VectorMultiTrajectoryBase::hasColumn_impl(*this, key);
  }

  // END INTERFACE
};

}  // namespace Acts
