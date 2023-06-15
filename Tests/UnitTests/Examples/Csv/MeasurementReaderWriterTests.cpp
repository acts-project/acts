// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "ActsExamples/Digitization/Smearers.hpp"
#include "ActsExamples/Io/Csv/CsvMeasurementReader.hpp"
#include "ActsExamples/Io/Csv/CsvMeasurementWriter.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Cluster.hpp"

#include <fstream>
#include <iostream>

using namespace ActsExamples;

const CsvMeasurementWriter::Config writerConfig;
const CsvMeasurementReader::Config readerConfig{
    writerConfig.outputDir,
    writerConfig.inputMeasurements,
    writerConfig.inputMeasurementSimHitsMap,
    "sourcelinks",
    writerConfig.inputClusters
};

// void write(const MeasurementContainer &m, const ClusterContainer &b;) {
//     ActsExamples::WhiteBoard board;
//     board.add
//
//     auto cfg = CsvMeasurementWriter::Config();
//
//
// }

BOOST_AUTO_TEST_CASE(RoundTrip) {


    // A: Write
    {

    }
    // B: Read
    {

    }

}
