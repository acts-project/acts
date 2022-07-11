// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

double fieldCorrectionDefault(unsigned region, double y, double r) {
   if (region == 999) return y+r; // 
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

// default with two slices, one for negative and one for positive z, counting some small overlaps, and -1 means "just take everything"
bool inSliceSPDefault(double z, unsigned layer, int slice) {
   if (slice == -1) return true;

   double absz = abs(z);
   if (slice == 0 && z > 50) return false;
   else if (slice == 1 && z < -50) return false;
   else {
      if (layer <= 3) {
         if (absz < 200) return true;
         else return false;
      }
      else if (layer == 4) {
         if (absz < 300) return true;
         else return false;
      }
      else if (layer == 5) {
         if (absz < 400) return true;
         else return false;
      }
      else if (layer == 6) {
         if (absz < 600) return true;
         else return false;
      }
      else if (layer == 7) {
         if (absz < 700) return true;
         else return false;
      }
      else if (layer == 8) {
         if (absz < 800) return true;
         else return false;
      }
      else if (layer == 9) {
         if (absz < 1100) return true;
         else return false;
      }
      else return false;
   }
}

// default with two slices, one for negative and one for positive z, counting some small overlaps, and -1 means "just take everything"
bool inSliceMeasurementDefault(double z, unsigned layer, int slice) {
   if (slice == -1) return true;

   double absz = abs(z);
   if (slice == 0 && z > 50) return false;
   else if (slice == 1 && z < -50) return false;
   else {
      if (layer <= 3) {
         if (absz < 200) return true;
         else return false;
      }
      else if (layer == 4) {
         if (absz < 300) return true;
         else return false;
      }
      else if (layer == 5) {
         if (absz < 400) return true;
         else return false;
      }
      else if (layer == 6) {
         if (absz < 600) return true;
         else return false;
      }
      else if (layer == 7) {
         if (absz < 700) return true;
         else return false;
      }
      else if (layer == 8) {
         if (absz < 800) return true;
         else return false;
      }
      else if (layer == 9) {
         if (absz < 1100) return true;
         else return false;
      }
      else return false;
   }
}

