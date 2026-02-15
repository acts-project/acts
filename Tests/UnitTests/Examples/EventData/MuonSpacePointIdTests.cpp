// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"

using Identifier = ActsExamples::MuonSpacePoint::MuonId;

using TechField = Identifier::TechField;
using StationName = Identifier::StationName;
using DetSide = Identifier::DetSide;
using namespace Acts;
namespace ActsExamples {

std::ostream& operator<<(std::ostream& ostr, const TechField& field) {
  ostr << Identifier::toString(field);
  return ostr;
}
std::ostream& operator<<(std::ostream& ostr, const StationName& field) {
  ostr << Identifier::toString(field);
  return ostr;
}
std::ostream& operator<<(std::ostream& ostr, const DetSide& field) {
  ostr << Identifier::toString(field);
  return ostr;
}

namespace Test {

BOOST_AUTO_TEST_SUITE(MuonSpacePointIdTests)

BOOST_AUTO_TEST_CASE(MainTest) {
  constexpr std::uint16_t sectorMax = 63;
  constexpr std::uint8_t layerMax = 15;
  constexpr std::uint16_t channelMax = 2;
  for (const auto tech : {TechField::Mm, TechField::sTgc, TechField::Mdt,
                          TechField::Rpc, TechField::Tgc}) {
    for (const auto side : {DetSide::A, DetSide::C}) {
      for (std::int8_t st = toUnderlying(StationName::UnDef) + 1;
           st < toUnderlying(StationName::MaxVal); ++st) {
        const auto stName = static_cast<StationName>(st);
        for (std::uint16_t sector = 1; sector <= sectorMax; ++sector) {
          /// Create a reference spacepoint object
          Identifier refId{};
          refId.setChamber(stName, side, sector, tech);
          BOOST_CHECK_EQUAL(refId.msStation(), stName);
          BOOST_CHECK_EQUAL(refId.side(), side);
          BOOST_CHECK_EQUAL(refId.sector(), sector);
          BOOST_CHECK_EQUAL(refId.technology(), tech);
          for (const auto& [measEta, measPhi, measTime] :
               std::vector<std::tuple<bool, bool, bool>>{{false, false, false},
                                                         {true, false, false},
                                                         {true, true, false},
                                                         {true, false, true},
                                                         {false, true, false},
                                                         {false, true, true},
                                                         {false, false, true},
                                                         {true, true, true}}) {
            refId.setCoordFlags(measEta, measPhi, measTime);
            BOOST_CHECK_EQUAL(refId.measuresEta(), measEta);
            BOOST_CHECK_EQUAL(refId.measuresPhi(), measPhi);
            BOOST_CHECK_EQUAL(refId.measuresTime(), measTime);

            for (std::uint8_t layer = 1; layer <= layerMax; ++layer) {
              for (std::uint16_t ch = 1; ch <= channelMax; ++ch) {
                refId.setLayAndCh(layer, ch);
                BOOST_CHECK_EQUAL(refId.channel(), ch);
                BOOST_CHECK_EQUAL(refId.detLayer(), layer);

                Identifier testId{refId.toInt()};
                BOOST_CHECK_EQUAL(refId.toInt(), testId.toInt());
                BOOST_CHECK_EQUAL(testId.msStation(), refId.msStation());
                BOOST_CHECK_EQUAL(testId.side(), refId.side());
                BOOST_CHECK_EQUAL(testId.technology(), refId.technology());
                BOOST_CHECK_EQUAL(testId.sector(), refId.sector());

                BOOST_CHECK_EQUAL(refId.measuresEta(), testId.measuresEta());
                BOOST_CHECK_EQUAL(refId.measuresPhi(), testId.measuresPhi());
                BOOST_CHECK_EQUAL(refId.measuresTime(), testId.measuresTime());
                BOOST_CHECK_EQUAL(testId.detLayer(), refId.detLayer());
                BOOST_CHECK_EQUAL(refId.sameStation(testId), true);
                BOOST_CHECK_EQUAL(testId.channel(), refId.channel());
              }
            }  /// layer flag
          }  // coord flag loop
        }  // sector loop
      }  // side loop
    }  // station index loop
  }  // tech index loop
}

BOOST_AUTO_TEST_SUITE_END()
}  // namespace Test
}  // namespace ActsExamples
