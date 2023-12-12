// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// @file HoughTransformSeeder.hpp
// @author Max Goblirsch, most of the work by Riley Xu and Jahred Adelman
// @brief Implements helpers to support Hough transforms for a variety of detectors

#include "HoughVectors.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include <set>
#include <unordered_set>

#pragma once

namespace Acts {
    namespace HoughTransformUtils{

        /// this type is responsible for counting the density in each bin of our hough hist 
        using yieldType = int;
        /// this type is responsible for encoding the coordinates of our hough hist
        using coordType = double; 

        // this type is used as an identifier proxy 
        using idType = long long int;
        
        /// @brief this function represents a mapping of a coordinate point in detector space to a line in 
        /// hough space. Given the value of the first hough coordinate, it will return the 
        /// corresponding second coordinate according to the line parametrisation. 
        /// @tparam PointType: The class representing a point in detector space (can differ between implementations)
        template<class PointType> using lineParametrisation = std::function<coordType(coordType, const PointType&)>; 

        /// @brief this function represents a way of obtaining an identifier for an input object to the hough.
        /// Will allow us to return a list of "hits on trajectory" for all of our hough candidates
        template<class PointType> using idGetter = std::function<idType(const PointType&)>; 

        /// @brief this function represents a way of obtaining a layer index for an input object to the hough.
        template<class PointType> using layerGetter = std::function<unsigned(const PointType&)>; 

        struct HoughPlaneConfig{
            unsigned nLayers = 10; 
            int nBinsX=0;
            int nBinsY=0; 
            coordType xMin = 0;
            coordType xMax = 0;
            coordType yMin = 0;
            coordType yMax = 0;
            double threshold = 4; 
            double fracForWidth = 0.5;

            int localMaxWindowSize = 2;  // Only create candidates from a local maximum
        };
        struct HoughMaximum{
            coordType x;    // x value of the maximum   
            coordType y;    // y value of the maximum
            coordType wx;   // x width of the maximum
            coordType wy;   // y width of the maximum
            std::set<idType> hitIdentifiers; 
            yieldType nHits; 
        };
        using HoughHist = vector2D<std::pair<yieldType, std::unordered_set<idType>>>;
        class HoughPlane{
            public:
                /// @brief instantiate the (empty) hough plane 
                /// @param cfg: configuration 
                HoughPlane(HoughPlaneConfig & cfg);
                /// @brief add one measurement to the hough plane 
                /// @tparam PointType: Type of the measurement
                /// @param measurement: The measurement to add
                /// @param linePar: The function y(x) parametrising the hough space line for a given measurement 
                /// @param widthPar: The function dy(x) parametrising the width of the y(x) curve 
                ///                   for a given measurement 
                /// @param idGetter: A function that maps our measurement to a (unique) identifier. 
                /// @param layerGetter: A function that maps our measurement to a layer index
                template <class PointType> void fill(const PointType & measurement, 
                                                lineParametrisation<PointType> linePar, 
                                                lineParametrisation<PointType> widthPar, 
                                                idGetter<PointType> idGetter,
                                                layerGetter<PointType> layerGetter); 
                /// performs a search for the maxima and returns a list of the found candidates
                std::vector<HoughMaximum> getMaxima(); 
                /// resets the contents of the grid, starting from scratch
                void reset();
                yieldType operator()(int x, int y){
                    mergeLayers();
                    return m_totalHoughHist(x,y).first; 
                }   
            private: 
                static inline int quant(double min, double max, unsigned nSteps, double val) {
                    return static_cast<int>((val - min) / (max - min) * nSteps);
                }

                // Returns the lower bound of the bin specified by step
                static inline double unquant(double min, double max, unsigned nSteps,
                                            int step) {
                    return min + (max - min) * step / nSteps;
                }

                bool passThreshold(HoughHist const& houghHist, unsigned x,
                     unsigned y) const;  // did we pass extensions?
                void drawHoughHist(HoughHist const& houghHist,
                                    std::string const& name);  // for making pretty plots
                std::vector<std::vector<int>> getComboIndices(std::vector<size_t>& sizes)
                    const;  // useful to find all candidates from given bins that pass
                            // (looping over hit combinatorics)
                void mergeLayers(); 
                HoughPlaneConfig m_cfg; 
                HoughHist m_totalHoughHist;
                std::vector<HoughHist> m_layerHists; 
                std::vector<coordType> m_binCentersX; 
                std::vector<coordType> m_binCentersY; 
                bool m_mergedLayers;
        };
    };
};

#include "HoughTransformUtils.ipp"
