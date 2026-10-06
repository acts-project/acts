// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/ButcherTableau.hpp"

#include <cmath>
#include <numeric>
#include <stdexcept>

namespace Acts {

ButcherTableau::ButcherTableau(std::string name, unsigned order,
                               unsigned embeddedOrder, std::vector<double> c,
                               std::vector<std::vector<double>> a,
                               std::vector<double> b,
                               std::vector<double> bEmbedded)
    : m_name(std::move(name)),
      m_order(order),
      m_embeddedOrder(embeddedOrder),
      m_c(std::move(c)),
      m_b(std::move(b)),
      m_bEmbedded(std::move(bEmbedded)) {
  const std::size_t s = m_c.size();
  if (s == 0 || m_b.size() != s || a.size() != s ||
      (!m_bEmbedded.empty() && m_bEmbedded.size() != s)) {
    throw std::invalid_argument("ButcherTableau " + m_name +
                                ": the sizes do not match");
  }
  if (m_c[0] != 0.) {
    throw std::invalid_argument("ButcherTableau " + m_name +
                                ": the first node must be zero");
  }

  m_a.assign(s * s, 0.);
  for (std::size_t i = 0; i < s; ++i) {
    if (a[i].size() != i) {
      throw std::invalid_argument("ButcherTableau " + m_name + ": row " +
                                  std::to_string(i) +
                                  " must have one coefficient per previous "
                                  "stage");
    }
    const double rowSum = std::accumulate(a[i].begin(), a[i].end(), 0.);
    if (std::abs(rowSum - m_c[i]) > 1e-14) {
      throw std::invalid_argument("ButcherTableau " + m_name + ": row " +
                                  std::to_string(i) +
                                  " does not sum to its node");
    }
    for (std::size_t j = 0; j < i; ++j) {
      m_a[i * s + j] = a[i][j];
    }
  }
}

std::shared_ptr<const ButcherTableau> ButcherTableau::classicalRk4() {
  static const auto tableau = std::make_shared<const ButcherTableau>(
      "ClassicalRk4", 4, 0, std::vector<double>{0., 1. / 2, 1. / 2, 1.},
      std::vector<std::vector<double>>{
          {}, {1. / 2}, {0., 1. / 2}, {0., 0., 1.}},
      std::vector<double>{1. / 6, 1. / 3, 1. / 3, 1. / 6},
      std::vector<double>{});
  return tableau;
}

std::shared_ptr<const ButcherTableau> ButcherTableau::dormandPrince54() {
  static const auto tableau = std::make_shared<const ButcherTableau>(
      "DormandPrince54", 5, 4,
      std::vector<double>{0., 1. / 5, 3. / 10, 4. / 5, 8. / 9, 1., 1.},
      std::vector<std::vector<double>>{
          {},
          {1. / 5},
          {3. / 40, 9. / 40},
          {44. / 45, -56. / 15, 32. / 9},
          {19372. / 6561, -25360. / 2187, 64448. / 6561, -212. / 729},
          {9017. / 3168, -355. / 33, 46732. / 5247, 49. / 176, -5103. / 18656},
          {35. / 384, 0., 500. / 1113, 125. / 192, -2187. / 6784, 11. / 84}},
      std::vector<double>{35. / 384, 0., 500. / 1113, 125. / 192, -2187. / 6784,
                          11. / 84, 0.},
      std::vector<double>{5179. / 57600, 0., 7571. / 16695, 393. / 640,
                          -92097. / 339200, 187. / 2100, 1. / 40});
  return tableau;
}

std::shared_ptr<const ButcherTableau> ButcherTableau::verner98() {
  // The exact coefficients are rationals and rationals times sqrt(6).
  const double sqrt6 = std::sqrt(6.);
  static const auto tableau = std::make_shared<const ButcherTableau>(
      "Verner98", 9, 8,
      std::vector<double>{
          0., 1731. / 50000, 7630049. / 53810000 - 983539. / 53810000 * sqrt6,
          22890147. / 107620000 - 2950617. / 107620000 * sqrt6, 561. / 1000,
          387. / 1000 - 129. / 2000 * sqrt6, 387. / 1000 + 129. / 2000 * sqrt6,
          129. / 200, 387. / 800, 6757. / 100000, 1. / 4, 0.6590650618730999,
          4103. / 5000, 2253. / 2500, 1., 1.},
      std::vector<std::vector<double>>{
          {},
          {1731. / 50000},
          {-177968356965557. / 1002427673820000 +
               14180534491313. / 250606918455000 * sqrt6,
           64021741529527. / 200485534764000 -
               7504450763411. / 100242767382000 * sqrt6},
          {22890147. / 430480000 - 2950617. / 430480000 * sqrt6, 0.,
           68670441. / 430480000 - 8851851. / 430480000 * sqrt6},
          {1.1541094857182945 + 0.35585143038209954 * sqrt6, 0.,
           -4.244643541842229 - 1.3853417041862812 * sqrt6,
           3.651534056123934 + 1.0294902738041818 * sqrt6},
          {11380823631. / 157617812000 - 339148869. / 39404453000 * sqrt6, 0.,
           0., 0.27509409584108324 - 0.04001311565984247 * sqrt6,
           165912282616977. / 4179075230308000 -
               33181894472511. / 2089537615154000 * sqrt6},
          {26523528363. / 231790900000 + 863255358. / 123138915625 * sqrt6, 0.,
           0., -0.04535066278469676 - 0.10221596690492023 * sqrt6,
           0.3181366802236169 - 0.09402893012545875 * sqrt6,
           -362925891. / 1690350537500 + 857800423623. / 3380701075000 * sqrt6},
          {43. / 600, 0., 0., 0., 0., 43. / 150 + 43. / 2400 * sqrt6,
           43. / 150 - 43. / 2400 * sqrt6},
          {7353. / 102400, 0., 0., 0., 0.,
           22833. / 102400 + 8901. / 204800 * sqrt6,
           22833. / 102400 - 8901. / 204800 * sqrt6, -3483. / 102400},
          {0.04836757646340647, 0., 0., 0., 0.,
           0.07238199692289805 - 0.01350979230006338 * sqrt6,
           0.07238199692289805 + 0.01350979230006338 * sqrt6,
           -0.021438652846483126, -0.10412291746271944},
          {-426968570497. / 54394415898750 -
               92754382349. / 12087647977500 * sqrt6,
           0., 0., 0., 0., 1. / 30,
           -0.024935138638394777 - 0.05640851783761993 * sqrt6,
           4389715333607. / 309890657317500 +
               92754382349. / 11477431752500 * sqrt6,
           4990058173976. / 83757096376875 +
               371017529396. / 9306344041875 * sqrt6,
           0.1757081936006537 + 0.016133382194788114 * sqrt6},
          {0.0013154210043370917 + 0.014523298816645871 * sqrt6, 0., 0., 0., 0.,
           -0.1470277986219964 + 0.00020806004636011267 * sqrt6,
           -0.03674538349491771 + 0.10655409400308863 * sqrt6,
           0.06041011360833149 - 0.0152954534911017 * sqrt6,
           0.18124121813885546 - 0.0754549897479561 * sqrt6,
           0.16148742604165614 - 0.030535009627036813 * sqrt6,
           0.43838406519683376},
          {-0.9529699587733558 + 0.19039423971301084 * sqrt6, 0., 0., 0., 0.,
           -5.729962318480114 - 0.23459593309005522 * sqrt6,
           -4.284206994781277 + 1.6342021387456511 * sqrt6,
           -2.187855250643843 - 0.200516857448822 * sqrt6,
           2.9418084517982352 - 0.9891826634558618 * sqrt6,
           2.345886196104345 - 0.40030092446392296 * sqrt6, 5.8850910885039465,
           2.8028087862720628},
          {1.0566810851455466 - 0.26052133561707513 * sqrt6, 0., 0., 0., 0.,
           6.4842136850429215 + 0.09811590253505763 * sqrt6,
           4.505949487169823 - 2.013233075118432 * sqrt6,
           2.671206837567015 + 0.27437237384412294 * sqrt6,
           -2.698361693987128 + 1.3535240826677508 * sqrt6,
           -2.271654661242113 + 0.547742051688576 * sqrt6, -6.099948804751011,
           -3.002206187889399, 0.2553202529443446},
          {-2.215207438880872 + 0.5861765116544452 * sqrt6, 0., 0., 0., 0.,
           -13.845681302507453 - 0.03742054273563755 * sqrt6,
           -9.394560230624968 + 4.346459957781801 * sqrt6,
           -13.17932865245738 - 0.6173415340948993 * sqrt6,
           6.965087119532673 - 3.045447403911066 * sqrt6,
           5.261792176682924 - 1.2324269886946428 * sqrt6, 13.367893803828643,
           14.396650486650687, -0.79758133317768, 0.4409353709534278},
          {3.7340160511235347 - 0.6842097292280843 * sqrt6, 0., 0., 0., 0.,
           20.39157340157072 + 0.8027648746808281 * sqrt6,
           15.196038809510034 - 5.832455817227349 * sqrt6,
           34.126030893706144 + 0.7205868462592434 * sqrt6,
           -12.149896310859955 + 3.554773523297931 * sqrt6,
           -8.389171072898177 + 1.4385403022174306 * sqrt6, -18.909803813543427,
           -34.26354448030452, 1.2647565216956427, 0., 0.}},
      std::vector<double>{0.014611976858423152, 0., 0., 0., 0., 0., 0.,
                          -0.3915211862331339, 0.23109325002895065,
                          0.12747667699928525, 0.2246434176204158,
                          0.5684352689748513, 0.058258715572158275,
                          0.13643174034822156, 0.030570139830827976, 0.},
      std::vector<double>{
          0.01996996514886773, 0., 0., 0., 0., 0., 0., 2.19149930494933,
          0.08857071848208438, 0.11405602348659656, 0.2533163805345107,
          -2.056564386240941, 0.340809679901312, 0., 0., 0.048342313738239585});
  return tableau;
}

}  // namespace Acts
