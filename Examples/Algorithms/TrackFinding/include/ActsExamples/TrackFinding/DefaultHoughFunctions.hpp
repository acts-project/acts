// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

double fieldCorrectionDefault(unsigned region, double y, double r) {
   return 0.0;
}

double findLayerIDSPDefault(double r) {
   if (r < 50) return 0;
   else if (r < 100) return 1;
   else if (r < 150) return 2;
   else if (r < 200) return 3;
   return 9999; /// for now pixel only, so these won't be used
}

double findLayerIDMeasurementDefault(double r) {
   if (r < 200) return 9999; // this is a pixel, ignore beacuse it will be a SP, this won't be used
   else if (r < 300) return 4;
   else if (r < 400) return 5;
   else if (r < 550) return 6;
   else if (r < 700) return 7;
   else if (r < 900) return 8;
   else if (r < 1100) return 9;
   return 9999; /// shouldn't be here, this won't be used
}

