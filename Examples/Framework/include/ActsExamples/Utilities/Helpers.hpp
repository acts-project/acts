// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <string>
#include <utility>

#include "TEfficiency.h"
#include "TProfile.h"

class TEfficiency;
class TH1D;
class TH1F;
class TH2F;
class TProfile;

namespace ActsExamples::PlotHelpers {


/// @brief book a 1D histogram
/// @param histName the name of histogram
/// @param histTitle the title of histogram
/// @param varBinning the binning info of variable
/// @return histogram pointer
TH1F* bookHisto(const std::string& histName, const std::string& histTitle,
                const Binning& varBinning);

/// @brief book a 2D histogram
/// @param histName the name of histogram
/// @param histTitle the title of histogram
/// @param varXBinning the binning info of variable at x axis
/// @param varYBinning the binning info of variable at y axis
/// @return histogram pointer
TH2F* bookHisto(const std::string& histName, const std::string& histTitle,
                const Binning& varXBinning, const Binning& varYBinning);

/// @brief fill a 1D histogram
/// @param hist histogram to fill
/// @param value value to fill
/// @param weight weight to fill
void fillHisto(TH1F* hist, float value, float weight = 1.0);

/// @brief fill a 2D histogram
/// @param hist histogram to fill
/// @param xValue x value to fill
/// @param yValue y value to fill
/// @param weight weight to fill
void fillHisto(TH2F* hist, float xValue, float yValue, float weight = 1.0);

/// @brief extract details, i.e. mean and width of a 1D histogram and fill
/// them into histograms
/// @param inputHist histogram to investigate
/// @param j  which bin number of meanHist and widthHist to fill
/// @param meanHist histogram to fill the mean value of inputHist
/// @param widthHist  histogram to fill the width value of inputHist
///
/// @todo  write specialized helper class to extract details of hists
void anaHisto(TH1D* inputHist, int j, TH1F* meanHist, TH1F* widthHist);

/// @brief book a 1D efficiency plot
/// @param effName the name of plot
/// @param effTitle the title of plot
/// @param varBinning the binning info of variable
/// @return TEfficiency pointer
TEfficiency* bookEff(const std::string& effName, const std::string& effTitle,
                     const Binning& varBinning);

/// @brief book a 2D efficiency plot
/// @param effName the name of plot
/// @param effTitle the title of plot
/// @param varXBinning the binning info of variable at x axis
/// @param varYBinning the binning info of variable at y axis
/// @return TEfficiency pointer
TEfficiency* bookEff(const std::string& effName, const std::string& effTitle,
                     const Binning& varXBinning, const Binning& varYBinning);

/// @brief fill a 1D efficiency plot
/// @param efficiency plot to fill
/// @param value value to fill
/// @param status bool to denote passed or not
void fillEff(TEfficiency* efficiency, float value, bool status);

/// @brief fill a 2D efficiency plot
/// @param efficiency plot to fill
/// @param xValue x value to fill
/// @param yValue y value to fill
/// @param status bool to denote passed or not
void fillEff(TEfficiency* efficiency, float xValue, float yValue, bool status);

/// @brief book a TProfile plot
/// @param profName the name of plot
/// @param profTitle the title of plot
/// @param varXBinning the binning info of variable at x axis
/// @param varYBinning the binning info of variable at y axis
/// @return TProfile pointer
TProfile* bookProf(const std::string& profName, const std::string& profTitle,
                   const Binning& varXBinning, const Binning& varYBinning);

/// @brief fill a TProfile plot
/// @param profile plot to fill
/// @param xValue  xvalue to fill
/// @param yValue  yvalue to fill
/// @param weight weight to fill
void fillProf(TProfile* profile, float xValue, float yValue,
              float weight = 1.0);

}  // namespace ActsExamples::PlotHelpers
