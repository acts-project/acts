// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/IAxis.hpp"
#include "Acts/Utilities/detail/grid_helper.hpp"

#include <iosfwd>

#include <boost/container/small_vector.hpp>

namespace Acts {

/// Common base class for all MultiAxis instances. This allows generice handling
/// such as for inspection.
class IMultiAxis {
 public:
  template <typename T>
  using SmallVector = boost::container::small_vector<T, 3>;

  using FlatIndex = std::size_t;
  using AnyMultiIndex = SmallVector<std::size_t>;
  using AnyPoint = SmallVector<double>;
  using AnyAxesVector = SmallVector<const IAxis*>;

  virtual ~IMultiAxis() = default;

  virtual std::size_t getNAxes() const = 0;

  virtual const IAxis& getAxis(std::size_t i) const = 0;

  virtual AnyMultiIndex getNBinsAny() const {
    AnyMultiIndex result;
    result.reserve(getNAxes());
    for (const IAxis& axis : *this) {
      result.push_back(axis.getNBins());
    }
    return result;
  }

  virtual std::size_t getNTotalBins(bool fullCounter = true) const {
    std::size_t result = 1;
    for (const IAxis& axis : *this) {
      result *= axis.getNBins() + (fullCounter ? 2 : 0);
    }
    return result;
  }

  virtual AnyAxesVector getAnyAxesVector() const {
    AnyAxesVector result;
    std::ranges::transform(*this, std::back_inserter(result),
                           [](const IAxis& axis) { return &axis; });
    return result;
  }

  virtual AnyPoint getMinPointAny() const {
    AnyPoint result;
    result.reserve(getNAxes());
    for (const IAxis& axis : *this) {
      result.push_back(axis.getMin());
    }
    return result;
  }

  virtual AnyPoint getMaxPointAny() const {
    AnyPoint result;
    result.reserve(getNAxes());
    for (const IAxis& axis : *this) {
      result.push_back(axis.getMax());
    }
    return result;
  }

  virtual bool isInsideAny(const AnyPoint& point) const {
    if (point.size() != getNAxes()) {
      throw std::invalid_argument("Invalid number of coordinates");
    }
    for (std::size_t i = 0; i < point.size(); ++i) {
      const IAxis& axis = getAxis(i);
      if (!axis.isInside(point[i])) {
        return false;
      }
    }
    return true;
  }

  virtual AnyPoint getLowerLeftBinCornerAny(
      const AnyMultiIndex& indices) const {
    AnyPoint result;
    result.reserve(getNAxes());
    for (std::size_t i = 0; i < indices.size(); ++i) {
      const IAxis& axis = getAxis(i);
      result.push_back(axis.getBinLowerBound(indices[i]));
    }
    return result;
  }

  virtual AnyPoint getUpperRightBinCornerAny(
      const AnyMultiIndex& indices) const {
    AnyPoint result;
    result.reserve(getNAxes());
    for (std::size_t i = 0; i < indices.size(); ++i) {
      const IAxis& axis = getAxis(i);
      result.push_back(axis.getBinUpperBound(indices[i]));
    }
    return result;
  }

  virtual AnyPoint getBinCenterAny(const AnyMultiIndex& indices) const {
    AnyPoint result;
    result.reserve(getNAxes());
    for (std::size_t i = 0; i < indices.size(); ++i) {
      const IAxis& axis = getAxis(i);
      result.push_back(axis.getBinCenter(indices[i]));
    }
    return result;
  }

  class iterator {
   public:
    using value_type = const IAxis;
    using difference_type = std::ptrdiff_t;
    using pointer = const IAxis*;
    using reference = const IAxis&;

    using iterator_category = std::random_access_iterator_tag;
    using iterator_concept = std::random_access_iterator_tag;

    constexpr iterator() noexcept = default;
    constexpr iterator(const IMultiAxis& multiAxis, std::size_t index) noexcept
        : m_multiAxis(&multiAxis), m_index(index) {}

    constexpr reference operator*() const {
      return m_multiAxis->getAxis(m_index);
    }
    constexpr iterator& operator++() noexcept {
      ++m_index;
      return *this;
    }
    constexpr iterator operator++(int) noexcept {
      auto tmp = *this;
      ++(*this);
      return tmp;
    }
    constexpr iterator& operator--() noexcept {
      --m_index;
      return *this;
    }
    constexpr iterator operator--(int) noexcept {
      auto tmp = *this;
      --(*this);
      return tmp;
    }
    constexpr iterator& operator+=(difference_type n) noexcept {
      m_index += n;
      return *this;
    }
    constexpr iterator& operator-=(difference_type n) noexcept {
      m_index -= n;
      return *this;
    }

   private:
    const IMultiAxis* m_multiAxis{};
    std::size_t m_index{};

    friend constexpr iterator operator+(iterator it,
                                        difference_type n) noexcept {
      return it += n;
    }

    friend constexpr iterator operator+(difference_type n,
                                        iterator it) noexcept {
      return it += n;
    }

    friend constexpr iterator operator-(iterator it,
                                        difference_type n) noexcept {
      return it -= n;
    }

    friend constexpr difference_type operator-(const iterator& lhs,
                                               const iterator& rhs) noexcept {
      return lhs.m_index - rhs.m_index;
    }

    friend constexpr auto operator<=>(const iterator& a,
                                      const iterator& b) noexcept {
      return a.m_index <=> b.m_index;
    }

    friend constexpr bool operator==(const iterator& a,
                                     const iterator& b) noexcept {
      return a.m_index == b.m_index;
    }
  };

  iterator begin() const { return iterator(*this, 0); }

  iterator end() const { return iterator(*this, getNAxes()); }

 protected:
  virtual void toStream(std::ostream& os) const {
    for (std::size_t i = 0; i < getNAxes(); ++i) {
      os << getAxis(i);
      if (i < getNAxes() - 1) {
        os << ", ";
      }
    }
  }

 private:
  friend bool operator==(const IMultiAxis& lhs, const IMultiAxis& rhs) {
    if (lhs.getNAxes() != rhs.getNAxes()) {
      return false;
    }
    return std::ranges::equal(lhs, rhs);
  }

  friend std::ostream& operator<<(std::ostream& os,
                                  const IMultiAxis& multiAxis) {
    multiAxis.toStream(os);
    return os;
  }
};

template <std::size_t _DIM>
class IMultiAxisND : public IMultiAxis {
 public:
  static constexpr std::size_t DIM = _DIM;

  using MultiIndex = std::array<std::size_t, DIM>;
  using Point = std::array<double, DIM>;
  using AnyAxesArray = std::array<const IAxis*, DIM>;
  using AnyAxesTuple = decltype(std::apply(
      [](auto&&... xs) { return std::tie(*xs...); }, AnyAxesArray{}));

  std::size_t getNAxes() const override { return DIM; }

  virtual AnyAxesArray getAnyAxesArray() const {
    AnyAxesArray result{};
    std::ranges::transform(*this, result.begin(),
                           [](const IAxis& axis) { return &axis; });
    return result;
  }

  virtual AnyAxesTuple getAnyAxesTuple() const {
    return std::apply([](auto&&... xs) { return std::tie(*xs...); },
                      getAnyAxesArray());
  }

  virtual MultiIndex getNBins() const {
    MultiIndex result{};
    for (std::size_t i = 0; i < DIM; ++i) {
      result[i] = getAxis(i).getNBins();
    }
    return result;
  }

  virtual Point getMinPoint() const {
    Point result{};
    for (std::size_t i = 0; i < DIM; ++i) {
      result[i] = getAxis(i).getMin();
    }
    return result;
  }

  virtual Point getMaxPoint() const {
    Point result{};
    for (std::size_t i = 0; i < DIM; ++i) {
      result[i] = getAxis(i).getMax();
    }
    return result;
  }

  virtual bool isInside(const Point& point) const {
    for (std::size_t i = 0; i < DIM; ++i) {
      if (!getAxis(i).isInside(point[i])) {
        return false;
      }
    }
    return true;
  }

  virtual Point getLowerLeftBinCorner(const MultiIndex& multiIndex) const {
    Point result{};
    for (std::size_t i = 0; i < DIM; ++i) {
      result[i] = getAxis(i).getBinLowerBound(multiIndex[i]);
    }
    return result;
  }

  virtual Point getUpperRightBinCorner(const MultiIndex& multiIndex) const {
    Point result{};
    for (std::size_t i = 0; i < DIM; ++i) {
      result[i] = getAxis(i).getBinUpperBound(multiIndex[i]);
    }
    return result;
  }

  virtual Point getBinCenter(const MultiIndex& multiIndex) const {
    Point result{};
    for (std::size_t i = 0; i < DIM; ++i) {
      result[i] = getAxis(i).getBinCenter(multiIndex[i]);
    }
    return result;
  }

  virtual FlatIndex getFlatIndexFromPosition(const Point& point) const {
    return getFlatIndexFromMultiIndex(getMultiIndexFromPosition(point));
  }

  virtual FlatIndex getFlatIndexFromMultiIndex(
      const MultiIndex& multiIndex) const {
    return detail::grid_helper::getGlobalBin(multiIndex, getAnyAxesTuple());
  }

  virtual MultiIndex getMultiIndexFromPosition(const Point& point) const {
    return detail::grid_helper::getLocalBinIndices(point, getAnyAxesTuple());
  }

  virtual MultiIndex getMultiIndexFromFlatIndex(FlatIndex flatIndex) const {
    return detail::grid_helper::getLocalBinIndices(flatIndex,
                                                   getAnyAxesTuple());
  }

  virtual detail::GlobalNeighborHoodIndices<DIM> getNeighborHoodIndices(
      const MultiIndex& multiIndex, std::size_t size = 1u) const = 0;

  virtual detail::GlobalNeighborHoodIndices<DIM> getNeighborHoodIndices(
      const MultiIndex& multiIndex,
      std::array<std::pair<int, int>, DIM>& sizePerAxis) const = 0;

  virtual detail::GlobalNeighborHoodIndices<DIM> getClosestPointsIndices(
      const MultiIndex& multiIndex) const = 0;

  virtual detail::GlobalNeighborHoodIndices<DIM> getClosestPointsIndices(
      const Point& position) const {
    return getClosestPointsIndices(getMultiIndexFromPosition(position));
  }
};

}  // namespace Acts
