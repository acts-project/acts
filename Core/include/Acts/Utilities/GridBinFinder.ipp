// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <std::size_t DIM>
template <typename... args>
Acts::GridBinFinder<DIM>::GridBinFinder(args&&... vals) {
  static_assert(sizeof...(args) == DIM);
  static_assert(
      std::conjunction<std::disjunction<
          std::is_same<int, typename std::decay<args>::type>,
          std::is_same<std::pair<int, int>, typename std::decay<args>::type>,
          std::is_same<std::vector<std::pair<int, int>>,
                       typename std::decay<args>::type>>...>::value);
  storeValue(std::forward<args>(vals)...);
}

template <std::size_t DIM>
template <typename first_value_t, typename... vals>
void Acts::GridBinFinder<DIM>::storeValue(first_value_t&& fv,
                                          vals&&... others) {
  constexpr std::size_t N = sizeof...(vals);
  static_assert(N < DIM);
  /// Check the fist value is reasonable
  using decayed_value_t = typename std::decay<first_value_t>::type;
  if constexpr (std::is_same<int, decayed_value_t>::value) {
    /// if int -> value is positive
    assert(fv >= 0);
    m_values[DIM - N - 1ul] = fv;
  } else if constexpr (std::is_same<std::pair<int, int>,
                                    decayed_value_t>::value) {
    m_values[DIM - N - 1ul] = fv;
  } else {
    /// If vector of pairs -> also allow for empty vectors
    if (!fv.empty()) {
      m_values[DIM - N - 1ul] = std::forward<first_value_t>(fv);
    } else {
      m_values[DIM - N - 1ul] = 1;
    }
  }
  if constexpr (N != 0ul) {
    storeValue(std::forward<vals>(others)...);
  }
}

template <std::size_t DIM>
std::array<std::pair<int, int>, DIM> Acts::GridBinFinder<DIM>::getSizePerAxis(
    const std::array<std::size_t, DIM>& locPosition) const {
  std::array<std::pair<int, int>, DIM> output;
  for (std::size_t i(0ul); i < DIM; ++i) {
    output[i] = std::visit(
        [&locPosition, i](const auto& val) -> std::pair<int, int> {
          using value_t = typename std::decay<decltype(val)>::type;
          if constexpr (std::is_same<int, value_t>::value) {
            assert(val >= 0);
            return std::make_pair(-val, val);
          } else if constexpr (std::is_same<std::pair<int, int>,
                                            value_t>::value) {
            return std::make_pair(-val.first, val.second);
          } else {
            assert(locPosition.size() > i);
            assert(locPosition[i] > 0ul);
            assert(val.size() >= locPosition[i]);
            return val[locPosition[i] - 1ul];
          }
        },
        m_values[i]);
  }
  return output;
}

template <std::size_t DIM>
template <typename stored_t, class... Axes>
boost::container::small_vector<std::size_t, Acts::detail::ipow(3, DIM)>
Acts::GridBinFinder<DIM>::findBins(
    const std::array<std::size_t, DIM>& locPosition,
    const Acts::Grid<stored_t, Axes...>& grid) const {
  static_assert(sizeof...(Axes) == DIM);
  assert(isGridCompatible(grid));
  std::array<std::pair<int, int>, DIM> sizePerAxis =
      getSizePerAxis(locPosition);
  return grid.neighborHoodIndices(locPosition, sizePerAxis).collect();
}

template <std::size_t DIM>
template <typename stored_t, class... Axes>
bool Acts::GridBinFinder<DIM>::isGridCompatible(
    const Acts::Grid<stored_t, Axes...>& grid) const {
  const std::array<std::size_t, DIM> nLocBins = grid.numLocalBins();
  for (std::size_t i(0ul); i < DIM; ++i) {
    std::size_t nBins = nLocBins[i];
    bool isCompabile = std::visit(
        [nBins](const auto& val) -> bool {
          using value_t = typename std::decay<decltype(val)>::type;
          if constexpr (std::is_same<int, value_t>::value ||
                        std::is_same<std::pair<int, int>, value_t>::value) {
            return true;
          } else {
            return val.size() == nBins;
          }
        },
        m_values[i]);
    if (!isCompabile) {
      return false;
    }
  }
  return true;
}
