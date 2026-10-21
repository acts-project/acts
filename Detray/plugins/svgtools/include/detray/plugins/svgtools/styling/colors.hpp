// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Actsvg include(s)
#include "actsvg/core/style.hpp"

// System include(s)
#include <array>
#include <random>
#include <stdexcept>
#include <vector>

namespace detray::svgtools::styling::colors {

/// @brief Picks a random element in the container.
template <typename container_t>
inline auto pick_random(container_t container) {
  typename container_t::value_type c{};
  std::sample(container.cbegin(), container.cend(), &c, 1,
              std::mt19937{std::random_device{}()});
  return c;
}

// Black
inline constexpr std::array black{0, 0, 0};
inline constexpr std::array dark_grey{171, 171, 171};
inline constexpr std::array mortar{89, 89, 89};
inline constexpr std::array suva_grey{137, 137, 137};
inline constexpr std::array very_light_grey{207, 207, 207};

// Red tones
inline constexpr std::array red{255, 0, 0};
inline constexpr std::array cardinal{167, 51, 63};
inline constexpr std::array madder{165, 28, 48};
inline constexpr std::array auburn{167, 51, 63};
inline constexpr std::array burgundy{116, 18, 29};
inline constexpr std::array chocolate_cosmos{88, 12, 31};
inline constexpr std::array macaroni_and_cheese{255, 188, 121};
inline constexpr std::array pumpkin{255, 128, 14};
inline constexpr std::array tenne{200, 82, 0};

// Blue tones
inline constexpr std::array blue{0, 0, 255};
inline constexpr std::array celestial_blue{62, 146, 204};
inline constexpr std::array cerulean{0, 107, 164};
inline constexpr std::array lapis_lazulli{42, 98, 143};
inline constexpr std::array picton_blue{95, 158, 209};
inline constexpr std::array prussian_blue1{19, 41, 61};
inline constexpr std::array prussian_blue2{22, 50, 79};
inline constexpr std::array indigo_dye{24, 67, 90};
inline constexpr std::array sail{162, 200, 236};

// Green tones
inline constexpr std::array green{0, 255, 0};
inline constexpr std::array celadon{190, 230, 206};
inline constexpr std::array aquamarine1{188, 255, 219};
inline constexpr std::array aquamarine2{141, 255, 205};
inline constexpr std::array emerald{104, 216, 155};
inline constexpr std::array shamrock_green{79, 157, 105};

// Yellow/orange tones
inline constexpr std::array yellow{0, 255, 0};

namespace gradient {

inline std::vector<actsvg::style::color> rainbow_scale{
    actsvg::style::color{{255, 0, 0}, 1.f},
    actsvg::style::color{{255, 51, 0}, 1.f},
    actsvg::style::color{{255, 102, 0}, 1.f},
    actsvg::style::color{{255, 153, 0}, 1.f},
    actsvg::style::color{{255, 204, 0}, 1.f},
    actsvg::style::color{{255, 255, 0}, 1.f},
    actsvg::style::color{{204, 255, 0}, 1.f},
    actsvg::style::color{{153, 255, 0}, 1.f},
    actsvg::style::color{{102, 255, 0}, 1.f},
    actsvg::style::color{{51, 255, 0}, 1.f},
    actsvg::style::color{{0, 255, 0}, 1.f},
    actsvg::style::color{{0, 255, 51}, 1.f},
    actsvg::style::color{{0, 255, 102}, 1.f},
    actsvg::style::color{{0, 255, 153}, 1.f},
    actsvg::style::color{{0, 255, 204}, 1.f},
    actsvg::style::color{{0, 255, 255}, 1.f},
    actsvg::style::color{{0, 204, 255}, 1.f},
    actsvg::style::color{{0, 153, 255}, 1.f},
    actsvg::style::color{{0, 102, 255}, 1.f},
    actsvg::style::color{{0, 51, 255}, 1.f}};

inline std::vector<actsvg::style::color> viridis_scale{
    actsvg::style::color{{253, 231, 37}, 1.f},
    actsvg::style::color{{221, 227, 24}, 1.f},
    actsvg::style::color{{186, 222, 40}, 1.f},
    actsvg::style::color{{149, 216, 64}, 1.f},
    actsvg::style::color{{117, 208, 84}, 1.f},
    actsvg::style::color{{86, 198, 103}, 1.f},
    actsvg::style::color{{61, 188, 116}, 1.f},
    actsvg::style::color{{41, 175, 127}, 1.f},
    actsvg::style::color{{32, 163, 134}, 1.f},
    actsvg::style::color{{31, 150, 139}, 1.f},
    actsvg::style::color{{35, 138, 141}, 1.f},
    actsvg::style::color{{40, 125, 142}, 1.f},
    actsvg::style::color{{45, 113, 142}, 1.f},
    actsvg::style::color{{51, 99, 141}, 1.f},
    actsvg::style::color{{57, 85, 140}, 1.f},
    actsvg::style::color{{64, 70, 136}, 1.f},
    actsvg::style::color{{69, 55, 129}, 1.f},
    actsvg::style::color{{72, 37, 118}, 1.f},
    actsvg::style::color{{72, 20, 103}, 1.f},
    actsvg::style::color{{68, 1, 84}, 1.f}};

inline std::vector<actsvg::style::color> plasma_scale{
    actsvg::style::color{{240, 249, 33}, 1.f},
    actsvg::style::color{{247, 226, 37}, 1.f},
    actsvg::style::color{{252, 205, 37}, 1.f},
    actsvg::style::color{{254, 183, 45}, 1.f},
    actsvg::style::color{{252, 163, 56}, 1.f},
    actsvg::style::color{{247, 144, 68}, 1.f},
    actsvg::style::color{{240, 127, 79}, 1.f},
    actsvg::style::color{{231, 110, 91}, 1.f},
    actsvg::style::color{{221, 94, 102}, 1.f},
    actsvg::style::color{{209, 78, 114}, 1.f},
    actsvg::style::color{{197, 64, 126}, 1.f},
    actsvg::style::color{{182, 48, 139}, 1.f},
    actsvg::style::color{{167, 33, 151}, 1.f},
    actsvg::style::color{{149, 17, 161}, 1.f},
    actsvg::style::color{{131, 5, 167}, 1.f},
    actsvg::style::color{{110, 0, 168}, 1.f},
    actsvg::style::color{{89, 1, 165}, 1.f},
    actsvg::style::color{{67, 3, 158}, 1.f},
    actsvg::style::color{{44, 5, 148}, 1.f},
    actsvg::style::color{{13, 8, 135}, 1.f},
};

/// Generate stops for actsvg gradients
inline std::vector<actsvg::style::gradient::stop> generate_stops(
    const std::vector<actsvg::style::color> &scale, unsigned int n_stops) {
  if (n_stops > scale.size()) {
    throw std::invalid_argument(
        "Too many gradient stops for given color scale! Color scale has "
        "only " +
        std::to_string(scale.size()) + " entries.");
  }

  std::vector<actsvg::style::gradient::stop> stops{};
  stops.reserve(n_stops);

  // Choose a color from the scale
  std::size_t color_step{static_cast<std::size_t>(scale.size() / n_stops)};
  std::size_t i_color{0u};

  // Find the gradient percentage for the color
  float grad_step{1.f / static_cast<float>(n_stops - 1u)};
  float grad{0.f};

  for (std::size_t i = 0u; i < n_stops; ++i) {
    stops.push_back(actsvg::style::gradient::stop{grad, scale.at(i_color)});

    i_color += color_step;
    grad += grad_step;
  }

  return stops;
}

}  // namespace gradient

inline std::vector<actsvg::style::color> black_theme(
    const actsvg::scalar opacity) {
  return {{black, opacity}};
}

inline std::vector<actsvg::style::color> red_theme(
    const actsvg::scalar opacity) {
  return {{cardinal, opacity},
          {madder, opacity},
          {auburn, opacity},
          {burgundy, opacity},
          {chocolate_cosmos, opacity}};
}

inline std::vector<actsvg::style::color> blue_theme(
    const actsvg::scalar opacity) {
  return {{celestial_blue, opacity},
          {lapis_lazulli, opacity},
          {prussian_blue1, opacity},
          {prussian_blue2, opacity},
          {indigo_dye, opacity}};
}

inline std::vector<actsvg::style::color> green_theme(
    const actsvg::scalar opacity) {
  return {{emerald, opacity},
          {shamrock_green, opacity},
          {celadon, opacity},
          {aquamarine1, opacity},
          {aquamarine2, opacity}};
}

// Same color circle that is used in matplot plugin
struct tableau_colorblind10 {
  static std::vector<actsvg::style::color> grey_tones(
      const actsvg::scalar opacity) {
    return {{dark_grey, opacity},
            {mortar, opacity},
            {suva_grey, opacity},
            {very_light_grey, opacity}};
  }
  static std::vector<actsvg::style::color> blue_tones(
      const actsvg::scalar opacity) {
    return {{cerulean, opacity}, {picton_blue, opacity}, {sail, opacity}};
  }
  static std::vector<actsvg::style::color> red_tones(
      const actsvg::scalar opacity) {
    return {
        {tenne, opacity}, {pumpkin, opacity}, {macaroni_and_cheese, opacity}};
  }
};

}  // namespace detray::svgtools::styling::colors
