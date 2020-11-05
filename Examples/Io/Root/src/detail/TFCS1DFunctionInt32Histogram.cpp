// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TFCS1DFunctionInt32Histogram.hpp"
#include "TMath.h"
#include "TH1.h"

#include <iostream>

TFCS1DFunctionInt32Histogram::TFCS1DFunctionInt32Histogram(const TH1* hist, std::string name, Acts::Logging::Level lvl) : m_logger(Acts::getDefaultLogger(name, lvl))
{  Initialize(hist);
}

void TFCS1DFunctionInt32Histogram::Initialize(const TH1* hist)
{
	// Retrieve the number of bins & borders
  const int nBins = hist->GetNbinsX();
  m_HistoBorders.resize(nBins+1);
  m_HistoContents.resize(nBins);
  
  // The integral of the original histogram
  const double histIntegral = hist->Integral();
  const double invHistIntegral = 1. / histIntegral;
  
  // Fill the cumulative histogram
  float integral = 0.;
  std::vector<double> temp_HistoContents(nBins);  
  int iBin;
  for (iBin = 0; iBin < nBins; iBin++){
    float binval = hist->GetBinContent(iBin + 1);
    // Avoid negative bin values
    if(binval < 0) {
		ACTS_DEBUG("Bin content is negative in histogram "<<hist->GetName()<<" : "<<hist->GetTitle()<<" binval="<<binval
			<<" "<<binval * invHistIntegral * 100<<"% of integral="<<histIntegral<<". Forcing bin to 0.");
      binval = 0.;
    }
    // Store the value
    integral += binval;
    temp_HistoContents[iBin] = integral;
  }
  
  // Ensure that content is available
  if(integral == 0.) {
	  ACTS_DEBUG("Histogram "<<hist->GetName()<<" : "<<hist->GetTitle()<<" integral="<<integral<<" is 0");
    m_HistoBorders.clear();
    m_HistoContents.clear();
    return;
  }

  // Set the bin borders
  for (iBin = 1; iBin <= nBins; iBin++) 
	m_HistoBorders[iBin - 1]=hist->GetXaxis()->GetBinLowEdge(iBin);
  m_HistoBorders[nBins]=hist->GetXaxis()->GetXmax();
  
  // Set the bin content
  const float invIntegral = 1. / integral;
  for(iBin = 0; iBin < nBins; ++iBin) {
    m_HistoContents[iBin] = s_MaxValue * (temp_HistoContents[iBin] * invIntegral);
  }
}

double TFCS1DFunctionInt32Histogram::rndToFunction(double rnd) const
{
	// Fast exit
  if(m_HistoContents.empty()) {
    return 0;
  }
  
  // Find the bin
  const HistoContent_t int_rnd = s_MaxValue * rnd;
  const auto it = std::upper_bound(m_HistoContents.begin(), m_HistoContents.end(), int_rnd);
  size_t iBin = std::min((size_t) std::distance(m_HistoContents.begin(), it), m_HistoContents.size() - 1);
  
  // Interpolate between neighbouring bins and return a diced intermediate value
  const HistoContent_t basecont = (iBin > 0 ? m_HistoContents[iBin - 1] : 0);
  const HistoContent_t dcont = m_HistoContents[iBin] - basecont;  
  return m_HistoBorders[iBin] + (m_HistoBorders[iBin + 1] - m_HistoBorders[iBin]) * (dcont > 0 ? (int_rnd-basecont) / dcont : 0.5);
}