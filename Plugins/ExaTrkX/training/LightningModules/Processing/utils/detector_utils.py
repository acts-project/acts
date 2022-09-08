import os
import time
import pickle
import logging
import argparse
import numpy as np
import pandas as pd
import functools

import trackml.dataset


#############################################
#               DETECTOR UTILS              #
#############################################
def load_detector(detector_path):
    detector_orig = pd.read_csv(detector_path)
    detector_pfx = detector_path.split(".")[0]
    detector_preproc = detector_pfx + ".pickle"
    try:
        print("Loading detector...")
        with open(detector_preproc, "rb") as f:
            detector = pickle.load(f)
        print("Detector loaded.")
    except:
        print("Failed to load preprocessed detector. Building...")
        detector = pd.read_csv(detector_path)
        detector = preprocess_detector(detector)
        with open(detector_preproc, "xb") as f:
            pickle.dump(detector, f)
        print("Detector preprocessed and saved.")
    return detector_orig, detector


def preprocess_detector(detector):
    thicknesses = Detector_Thicknesses(detector).get_thicknesses()
    rotations = Detector_Rotations(detector).get_rotations()
    pixel_size = Detector_Pixel_Size(detector).get_pixel_size()
    det = dict(thicknesses=thicknesses, rotations=rotations, pixel_size=pixel_size)
    return det


def determine_array_size(detector):
    max_v, max_l, max_m = (0, 0, 0)
    unique_vols = detector.volume_id.unique()
    max_v = max(unique_vols) + 1
    for v in unique_vols:
        vol = detector.loc[detector["volume_id"] == v]
        unique_layers = vol.layer_id.unique()
        max_l = max(max_l, max(unique_layers) + 1)
        for l in unique_layers:
            lay = vol.loc[vol["layer_id"] == l]
            unique_modules = lay.module_id.unique()
            max_m = max(max_m, max(unique_modules) + 1)
    return max_v, max_l, max_m


class Detector_Rotations(object):
    def __init__(self, detector):
        self.detector = detector
        self.max_v, self.max_l, self.max_m = determine_array_size(detector)

    def get_rotations(self):
        print("  Extracting rotations...")
        self._init_rotation_array()
        self._extract_all_rotations()
        print("  Done.")
        return self.rot

    def _init_rotation_array(self):
        self.rot = np.zeros((self.max_v, self.max_l, self.max_m, 3, 3))

    def _extract_all_rotations(self):
        for i, r in self.detector.iterrows():
            v, l, m = tuple(map(int, (r.volume_id, r.layer_id, r.module_id)))
            rot = self._extract_rotation_matrix(r)
            self.rot[v, l, m] = rot

    def _extract_rotation_matrix(self, mod):
        """
        Extract the rotation matrix from module dataframe
        """
        r = np.matrix(
            [
                [mod.rot_xu.item(), mod.rot_xv.item(), mod.rot_xw.item()],
                [mod.rot_yu.item(), mod.rot_yv.item(), mod.rot_yw.item()],
                [mod.rot_zu.item(), mod.rot_zv.item(), mod.rot_zw.item()],
            ]
        )
        return r


class Detector_Thicknesses(object):
    def __init__(self, detector):
        self.detector = detector
        self.max_v, self.max_l, self.max_m = determine_array_size(detector)

    def get_thicknesses(self):
        print("  Extracting thicknesses...")
        self._init_thickness_array()
        self._extract_all_thicknesses()
        print("  Done.")
        return self.all_t

    def _init_thickness_array(self):
        self.all_t = np.zeros((self.max_v, self.max_l, self.max_m))

    def _extract_all_thicknesses(self):
        for i, r in self.detector.iterrows():
            v, l, m = tuple(map(int, (r.volume_id, r.layer_id, r.module_id)))
            self.all_t[v, l, m] = r.module_t


class Detector_Pixel_Size(object):
    def __init__(self, detector):
        print(detector.keys())
        self.detector = detector
        self.max_v, self.max_l, self.max_m = determine_array_size(detector)

    def get_pixel_size(self):
        print("  Extracting thicknesses...")
        self._init_size_array()
        self._extract_all_size()
        print("  Done.")
        return self.all_s

    def _init_size_array(self):
        self.all_s = np.zeros((self.max_v, self.max_l, self.max_m, 2))

    def _extract_all_size(self):
        for i, r in self.detector.iterrows():
            v, l, m = tuple(map(int, (r.volume_id, r.layer_id, r.module_id)))
            self.all_s[v, l, m, 0] = r.pitch_u
            self.all_s[v, l, m, 1] = r.pitch_v
