# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import json
from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.parametrize("legacy", [False, True])
@pytest.mark.parametrize("legacy_config", [False, True])
def test_material_map_configuration(tmp_path, legacy, legacy_config):
    scripts = Path(__file__).resolve().parents[3] / "Examples/Scripts/MaterialMapping"
    entries = []
    for identifiers, surface_type, directions in (
        ({"layer": 2}, "CylinderSurface", ("RPhi", "Z")),
        ({"layer": 2, "approach": 1}, "CylinderSurface", ("RPhi", "Z")),
        ({"boundary": 1}, "DiscSurface", ("R", "Phi")),
        ({"layer": 2, "sensitive": 1}, "PlaneSurface", ("X", "Y")),
    ):
        material = {
            "type": "proto-grid",
            "axis_specs": [
                {"type": "equidistant", "bins": 1, "direction": "Axis" + direction}
                for direction in directions
            ],
            "mapMaterial": True,
            "mappingType": "Default",
        }
        if legacy:
            material["type"] = "proto"
            material["binUtility"] = {
                "binningdata": [
                    {"bins": 1, "value": "bin" + direction.replace("RPhi", "Phi")}
                    for direction in directions
                ]
            }
            del material["axis_specs"]
        entries.append(
            {
                "volume": 1,
                **identifiers,
                "value": {
                    "type": surface_type,
                    "bounds": {"type": surface_type + "Bounds"},
                    "material": material,
                },
            }
        )
    geometry = tmp_path / "geometry.json"
    config_file = tmp_path / "config.json"
    document = {"Surfaces": {"entries": entries}, "Volumes": {"entries": []}}
    geometry.write_text(json.dumps(document))
    subprocess.run(
        [sys.executable, scripts / "writeMapConfig.py", geometry, config_file],
        check=True,
    )
    config = json.loads(config_file.read_text())
    assert len(config["Surfaces"]["1"]) == 4
    for entry in config["Surfaces"]["1"]:
        material = entry["value"]["material"]
        assert material["mapMaterial"] is False
        material["mapMaterial"] = True
        material["mappingType"] = "PostMapping"
        axes = material.get(
            "axis_specs", material.get("binUtility", {}).get("binningdata")
        )
        for i, axis in enumerate(axes):
            axis["bins"] = (3, 5)[i]
        # Exercise both representations and reordered configs independently
        # of the geometry's format, including old cylinder phi names.
        if legacy_config:
            material.pop("axis_specs", None)
            material["binUtility"] = {
                "binningdata": [
                    {
                        "bins": axis["bins"],
                        "value": axis.get("value", axis.get("direction", ""))
                        .replace("Axis", "bin")
                        .replace("RPhi", "Phi"),
                    }
                    for axis in reversed(axes)
                ]
            }
        else:
            material.pop("binUtility", None)
            material["axis_specs"] = [
                {
                    "type": "equidistant",
                    "bins": axis["bins"],
                    "direction": axis.get("direction", axis.get("value", "")).replace(
                        "bin", "Axis"
                    ),
                }
                for axis in reversed(axes)
            ]
    config_file.write_text(json.dumps(config))
    subprocess.run(
        [sys.executable, scripts / "configureMap.py", geometry, config_file], check=True
    )
    configured = json.loads(geometry.read_text())
    for entry in configured["Surfaces"]["entries"]:
        material = entry["value"]["material"]
        assert material["mapMaterial"] is True
        assert material["mappingType"] == "PostMapping"
        axes = material.get(
            "axis_specs", material.get("binUtility", {}).get("binningdata")
        )
        assert [axis["bins"] for axis in axes] == [3, 5]
