# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import copy
import json
from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.parametrize("legacy_map", [False, True])
@pytest.mark.parametrize("legacy_config", [False, True])
def test_material_mapping_configuration(tmp_path, legacy_map, legacy_config):
    scripts = Path(__file__).resolve().parents[3] / "Examples/Scripts/MaterialMapping"

    def material(legacy, bins):
        result = {"type": "proto", "mappingType": "Default", "mapMaterial": False}
        if legacy:
            result["binUtility"] = {
                "binningdata": [
                    {"value": "AxisPhi", "bins": bins[1], "min": -1, "max": 1},
                    {"value": "AxisR", "bins": bins[0], "min": 5, "max": 10},
                ]
            }
        else:
            result["axis_specs"] = [
                {"direction": "AxisR", "type": "equidistant", "bins": bins[0]},
                {"direction": "AxisPhi", "type": "equidistant", "bins": bins[1]},
            ]
        return result

    entries = []
    for identifiers in (
        {"layer": 2},
        {"boundary": 1},
        {"layer": 2, "approach": 1},
        {"layer": 2, "sensitive": 1},
    ):
        entries.append(
            {
                "volume": 1,
                **identifiers,
                "value": {
                    "bounds": {"type": "RadialBounds"},
                    "material": material(legacy_map, (1, 1)),
                },
            }
        )
    geometry = {"Surfaces": {"entries": entries}, "Volumes": {"entries": []}}
    map_path = tmp_path / "geometry.json"
    config_path = tmp_path / "config.json"
    map_path.write_text(json.dumps(geometry))
    subprocess.run(
        [sys.executable, scripts / "writeMapConfig.py", map_path, config_path],
        check=True,
    )
    config = json.loads(config_path.read_text())
    assert len(config["Surfaces"]["1"]) == 4
    for entry in config["Surfaces"]["1"]:
        assert entry["value"]["material"] == material(legacy_map, (1, 1))
        entry["value"]["material"] = material(legacy_config, (7, 11))
        entry["value"]["material"]["mapMaterial"] = True
    config_path.write_text(json.dumps(config))
    subprocess.run(
        [sys.executable, scripts / "configureMap.py", map_path, config_path],
        check=True,
    )
    expected = copy.deepcopy(geometry)
    for entry in expected["Surfaces"]["entries"]:
        entry["value"]["material"] = material(legacy_map, (7, 11))
        entry["value"]["material"]["mapMaterial"] = True
    assert json.loads(map_path.read_text()) == expected
