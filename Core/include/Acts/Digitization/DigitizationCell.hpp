// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
namespace Acts {

/// @brief pair of ints for definition of a cell
struct DigitizationCell final {
  // identification and data
  std::size_t channel0 = 0;
  std::size_t channel1 = 1;
  float data = 0.;

  // construct them
  DigitizationCell(std::size_t ch0, std::size_t ch1, float d = 0.)
      : channel0(ch0), channel1(ch1), data(d) {}

  /// To merge cells in case they are at the same position
  /// @param dc the cell to be added to the current cell
  /// @param analogueReadout flag indicating if we have analogue readout
  /// @note this function is needed because possible derived classes may
  /// calculate the energy deposit differently. Furthermore this allows to apply
  /// an energy cut, because the energy deposit can also be stored for digital
  /// readout.
  void addCell(const DigitizationCell& dc, bool analogueReadout) {
    if (analogueReadout) {
      data += dc.data;
    }
  }

  /// the deposited energy
  /// @note this function is needed because possible derived classes may
  /// calculate the energy deposit differently. Furthermore this allows to apply
  /// an energy cut, because the energy deposit can also be stored for digital
  /// readout.
  double depositedEnergy() const { return data; }
};

/// @brief DigitizationStep for further handling
struct DigitizationStep final {
  double stepLength{0.};   /// this is the path length within the cell
  double driftLength{0.};  /// this is the path length of the setp center to the
                           /// readout surface
  DigitizationCell stepCell;     /// this is the cell identifier of the segment
  Vector3 stepEntry;             /// this is the Entry point into the segment
  Vector3 stepExit;              /// this is the Exit point from the segment
  Vector2 stepReadoutProjected;  /// this is the projected position at the
                                 /// readout surface
  Vector2 stepCellCenter;        /// this is the cell position

  /// Standard constructor
  DigitizationStep()
      : stepCell(0, 0),
        stepEntry(0., 0., 0.),
        stepExit(0., 0., 0.),
        stepReadoutProjected(0., 0.),
        stepCellCenter(0., 0.) {}

  /// Constructor with arguments
  ///
  /// @param sl step length of this step
  /// @param dl drift length of this step
  /// @param dc is the digitization zell (with indices)
  /// @param entryP is the entry position into the cell
  /// @param exitP is the exit position from the cell
  /// @param projectedPosition is the position on the readout surface
  /// @param cellPosition is the nominal position of the cell
  DigitizationStep(double sl, double dl, const DigitizationCell& dc,
                   const Vector3& entryP, const Vector3& exitP,
                   const Vector2& projectedPosition,
                   const Vector2& cellPosition)
      : stepLength(sl),
        driftLength(dl),
        stepCell(dc),
        stepEntry(entryP),
        stepExit(exitP),
        stepReadoutProjected(projectedPosition),
        stepCellCenter(cellPosition) {}
};

}  // namespace Acts
