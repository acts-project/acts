// This file is part of the Acts project.
//
// Copyright (C) 2020-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/SubspaceHelpers.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/EventData/detail/CalculateResiduals.hpp"
#include "Acts/EventData/detail/ParameterTraits.hpp"
#include "Acts/EventData/detail/PrintParameters.hpp"

#include <array>
#include <cstddef>
#include <iosfwd>
#include <type_traits>
#include <variant>
#include <vector>

#include <boost/container/static_vector.hpp>

namespace ActsExamples {

template <typename Derived, std::size_t FullSize, bool ReadOnly>
class MeasurementProxyBase;
template <std::size_t FullSize, std::size_t Size, bool ReadOnly>
class FixedMeasurementProxy;
template <std::size_t FullSize, bool ReadOnly>
class VariableMeasurementProxy;

template <std::size_t Size>
using FixedBoundMeasurementProxy =
    FixedMeasurementProxy<Acts::eBoundSize, Size, false>;
template <std::size_t Size>
using ConstFixedBoundMeasurementProxy =
    FixedMeasurementProxy<Acts::eBoundSize, Size, true>;
using VariableBoundMeasurementProxy =
    VariableMeasurementProxy<Acts::eBoundSize, false>;
using ConstVariableBoundMeasurementProxy =
    VariableMeasurementProxy<Acts::eBoundSize, true>;

class MeasurementContainer {
 public:
  template <std::size_t Size>
  using FixedProxy = FixedMeasurementProxy<Acts::eBoundSize, Size, false>;
  template <std::size_t Size>
  using ConstFixedProxy = FixedMeasurementProxy<Acts::eBoundSize, Size, true>;
  using VariableProxy = VariableMeasurementProxy<Acts::eBoundSize, false>;
  using ConstVariableProxy = VariableMeasurementProxy<Acts::eBoundSize, true>;

  MeasurementContainer();

  std::size_t size() const;

  void reserve(std::size_t size);

  std::size_t addMeasurement(std::uint8_t size);

  VariableProxy getMeasurement(std::size_t index);
  ConstVariableProxy getMeasurement(std::size_t index) const;

  template <std::size_t Size>
  FixedProxy<Size> getMeasurement(std::size_t index) {
    return FixedProxy<Size>{*this, index};
  }
  template <std::size_t Size>
  ConstFixedProxy<Size> getMeasurement(std::size_t index) const {
    return ConstFixedProxy<Size>{*this, index};
  }

  VariableProxy makeMeasurement(std::uint8_t size);
  template <std::size_t Size>
  FixedProxy<Size> makeMeasurement() {
    return getMeasurement<Size>(addMeasurement(Size));
  }

  template <bool Const>
  class IteratorImpl {
   public:
    using value_type =
        std::conditional_t<Const, ConstVariableProxy, VariableProxy>;
    using reference = value_type;
    using pointer = value_type*;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::forward_iterator_tag;

    using Container = std::conditional_t<Const, const MeasurementContainer,
                                         MeasurementContainer>;

    IteratorImpl(Container& container, std::size_t index)
        : m_container(container), m_index(index) {}

    reference operator*() const {
      if constexpr (Const) {
        return m_container.getMeasurement(m_index);
      } else {
        return m_container.getMeasurement(m_index);
      }
    }

    pointer operator->() const { return &operator*(); }

    IteratorImpl& operator++() {
      ++m_index;
      return *this;
    }

    IteratorImpl operator++(int) {
      auto copy = *this;
      ++*this;
      return copy;
    }

    bool operator==(const IteratorImpl& other) const {
      return m_index == other.m_index;
    }

    bool operator!=(const IteratorImpl& other) const {
      return !(*this == other);
    }

   private:
    Container& m_container;
    std::size_t m_index;
  };

  using Iterator = IteratorImpl<false>;
  using ConstIterator = IteratorImpl<true>;

  Iterator begin();
  Iterator end();
  ConstIterator begin() const;
  ConstIterator end() const;

 public:
  struct MeasurementEntry {
    std::size_t subspaceIndexOffset{};
    std::size_t parameterOffset{};
    std::size_t covarianceOffset{};
    std::uint8_t size{};
  };

  std::vector<MeasurementEntry> m_entries;

  std::vector<Acts::SourceLink> m_sourceLinks;
  std::vector<std::uint8_t> m_subspaceIndices;
  std::vector<double> m_parameters;
  std::vector<double> m_covariances;
};

template <typename Derived, std::size_t FullSize, bool ReadOnly>
class MeasurementProxyBase {
 public:
  using Index = std::uint8_t;
  using Scalar = double;

  using FullVector = Acts::ActsVector<FullSize>;
  using FullSquareMatrix = Acts::ActsSquareMatrix<FullSize>;

  using Container = std::conditional_t<ReadOnly, const MeasurementContainer,
                                       MeasurementContainer>;

  MeasurementProxyBase(Container& container_, std::size_t index_)
      : m_container(&container_), m_index(index_) {}
  template <typename OtherDerived, bool OtherReadOnly>
  MeasurementProxyBase(
      const MeasurementProxyBase<OtherDerived, FullSize, OtherReadOnly>& other)
    requires(ReadOnly == OtherReadOnly || ReadOnly)
      : m_container(&other.container()), m_index(other.index()) {}

  Container& container() const { return *m_container; }
  std::size_t index() const { return m_index; }

  std::size_t size() const { return container().m_entries.at(m_index).size; }

  template <typename indices_t>
  bool contains(indices_t i) const {
    return self().subspaceHelper().contains(i);
  }

  template <typename indices_t>
  std::size_t indexOf(indices_t i) const {
    return self().subspaceHelper().indexOf(i);
  }

  Acts::SourceLink& sourceLink() {
    return container().m_sourceLinks.at(m_index);
  }

  const Acts::SourceLink& sourceLink() const {
    return container().m_sourceLinks.at(m_index);
  }

  template <typename IndexContainer>
  void setSubspaceIndices(const IndexContainer& indices)
    requires(!ReadOnly)
  {
    assert(checkSubspaceIndices(indices, FullSize, size()) &&
           "Invalid indices");
    std::transform(indices.begin(), indices.end(),
                   self().subspaceIndexVector().begin(),
                   [](auto index) { return static_cast<Index>(index); });
  }

  FullVector fullParameters() const {
    return self().subspaceHelper().expandVector(self().parameters());
  }

  FullSquareMatrix fullCovariance() const {
    return self().subspaceHelper().expandMatrix(self().covariance());
  }

  template <typename OtherDerived>
  void copyFrom(const OtherDerived& other)
    requires(!ReadOnly)
  {
    assert(size() == other.size() && "Size mismatch");
    sourceLink() = other.sourceLink();
    self().subspaceIndexVector() = other.subspaceIndexVector();
    self().parameters() = other.parameters();
    self().covariance() = other.covariance();
  }

 protected:
  Derived& self() { return static_cast<Derived&>(*this); }
  const Derived& self() const { return static_cast<const Derived&>(*this); }

  Container* m_container;
  std::size_t m_index;
};

template <std::size_t FullSize, std::size_t Size, bool ReadOnly>
class FixedMeasurementProxy
    : public MeasurementProxyBase<
          FixedMeasurementProxy<FullSize, Size, ReadOnly>, FullSize, ReadOnly> {
 public:
  using Base =
      MeasurementProxyBase<FixedMeasurementProxy<FullSize, Size, ReadOnly>,
                           FullSize, ReadOnly>;
  using Index = Base::Index;
  using Scalar = Base::Scalar;
  using Container = Base::Container;

  using SubspaceHelper = Acts::FixedSubspaceHelper<FullSize, Size>;

  using SubspaceVector = Eigen::Matrix<Index, Size, 1>;
  using SubspaceVectorMap =
      std::conditional_t<ReadOnly, Eigen::Map<const SubspaceVector>,
                         Eigen::Map<SubspaceVector>>;

  using ParametersVector = Eigen::Matrix<Scalar, Size, 1>;
  using ParametersVectorMap =
      std::conditional_t<ReadOnly, Eigen::Map<const ParametersVector>,
                         Eigen::Map<ParametersVector>>;

  using CovarianceMatrix = Eigen::Matrix<Scalar, Size, Size>;
  using CovarianceMatrixMap =
      std::conditional_t<ReadOnly, Eigen::Map<const CovarianceMatrix>,
                         Eigen::Map<CovarianceMatrix>>;

  FixedMeasurementProxy(Container& container_, std::size_t index_)
      : Base(container_, index_) {
    assert(container().m_entries.at(index()).size == Size && "Size mismatch");
  }
  template <typename OtherDerived, bool OtherReadOnly>
  FixedMeasurementProxy(
      const MeasurementProxyBase<OtherDerived, FullSize, OtherReadOnly>& other)
    requires(ReadOnly == OtherReadOnly || ReadOnly)
      : Base(other) {
    assert(container().m_entries.at(index()).size == Size && "Size mismatch");
  }

  using Base::container;
  using Base::index;

  static constexpr std::size_t size() { return Size; }

  SubspaceHelper subspaceHelper() const {
    return SubspaceHelper{subspaceIndexVector()};
  }

  Acts::SubspaceIndices<Size> subspaceIndices() const {
    return subspaceHelper().indices();
  }

  SubspaceVectorMap subspaceIndexVector() const {
    return SubspaceVectorMap{
        container().m_subspaceIndices.data() +
        container().m_entries.at(index()).subspaceIndexOffset};
  }

  ParametersVectorMap parameters() const {
    return ParametersVectorMap{
        container().m_parameters.data() +
        container().m_entries.at(index()).parameterOffset};
  }

  CovarianceMatrixMap covariance() const {
    return CovarianceMatrixMap{
        container().m_covariances.data() +
        container().m_entries.at(index()).covarianceOffset};
  }
};

template <std::size_t FullSize, bool ReadOnly>
class VariableMeasurementProxy
    : public MeasurementProxyBase<VariableMeasurementProxy<FullSize, ReadOnly>,
                                  FullSize, ReadOnly> {
 public:
  using Base =
      MeasurementProxyBase<VariableMeasurementProxy<FullSize, ReadOnly>,
                           FullSize, ReadOnly>;
  using Index = Base::Index;
  using Scalar = Base::Scalar;
  using Container = Base::Container;

  using SubspaceHelper = Acts::VariableSubspaceHelper<FullSize>;

  using SubspaceVector = Eigen::Matrix<Index, Eigen::Dynamic, 1>;
  using SubspaceVectorMap =
      std::conditional_t<ReadOnly, Eigen::Map<const SubspaceVector>,
                         Eigen::Map<SubspaceVector>>;

  using ParametersVector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
  using ParametersVectorMap =
      std::conditional_t<ReadOnly, Eigen::Map<const ParametersVector>,
                         Eigen::Map<ParametersVector>>;

  using CovarianceMatrix =
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  using CovarianceMatrixMap =
      std::conditional_t<ReadOnly, Eigen::Map<const CovarianceMatrix>,
                         Eigen::Map<CovarianceMatrix>>;

  VariableMeasurementProxy(Container& container_, std::size_t index_)
      : Base(container_, index_) {}
  template <typename OtherDerived, bool OtherReadOnly>
  VariableMeasurementProxy(
      const MeasurementProxyBase<OtherDerived, FullSize, OtherReadOnly>& other)
    requires(ReadOnly == OtherReadOnly || ReadOnly)
      : Base(other) {}

  using Base::container;
  using Base::index;

  SubspaceHelper subspaceHelper() const {
    return SubspaceHelper{subspaceIndexVector()};
  }

  SubspaceVectorMap subspaceIndexVector() const {
    return SubspaceVectorMap{
        container().m_subspaceIndices.data() +
            container().m_entries.at(index()).subspaceIndexOffset,
        static_cast<Eigen::Index>(this->size())};
  }

  ParametersVectorMap parameters() const {
    return ParametersVectorMap{
        container().m_parameters.data() +
            container().m_entries.at(index()).parameterOffset,
        static_cast<Eigen::Index>(this->size())};
  }

  CovarianceMatrixMap covariance() const {
    const auto size = this->size();
    return CovarianceMatrixMap{
        container().m_covariances.data() +
            container().m_entries.at(index()).covarianceOffset,
        static_cast<Eigen::Index>(size), static_cast<Eigen::Index>(size)};
  }
};

}  // namespace ActsExamples
