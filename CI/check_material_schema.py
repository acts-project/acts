#!/usr/bin/env python3
"""Offline schema/fixture QA, not a production material reader.

Requires jsonschema>=4.18. Run from any directory. No network references are used.
Direct execution validates the schema and examples. Run the self-tests with
``python -m pytest CI/check_material_schema.py`` (also requires pytest).
"""

import copy
import json
from pathlib import Path

from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = ROOT / "Plugins/Json/schema/material-map-v1.schema.json"
EXAMPLES = ROOT / "docs/examples/material-map-v1"
ID_COMPONENTS = ("volume", "portal", "layer", "passive", "sensitive", "extra")


def load(path):
    return json.loads(path.read_text())


def changed(document, path, value):
    result = copy.deepcopy(document)
    cursor = result
    for key in path[:-1]:
        cursor = cursor[key]
    cursor[path[-1]] = value
    return result


def validation_inputs():
    schema = load(SCHEMA)
    Draft202012Validator.check_schema(schema)
    validator = Draft202012Validator(schema)
    examples = {p.stem: load(p) for p in sorted(EXAMPLES.glob("*.json"))}
    return validator, examples


def main():
    validator, examples = validation_inputs()
    for name, document in examples.items():
        validator.validate(document)
        print(f"Validated {name}.json (schema)")


def test_schema_and_examples():
    main()


def test_valid_geometry_ids():
    validator, examples = validation_inputs()
    m = examples["minimal"]
    # All-zero, omitted-zero and full-width component IDs are valid objects.
    maxima = dict(zip(ID_COMPONENTS, (255, 255, 4095, 255, 1048575, 255)))
    for identifier in ({}, {"volume": 0}, maxima):
        document = changed(m, ["surfaces", 0, "target", "geometry_id"], identifier)
        validator.validate(document)


def test_invalid_structures():
    validator, examples = validation_inputs()
    m, s = (examples[n] for n in ("minimal", "surfaces"))
    maxima = dict(zip(ID_COMPONENTS, (255, 255, 4095, 255, 1048575, 255)))
    structural = [
        changed(m, ["header", "version"], 2),
        changed(m, ["header", "units"], {}),
        changed(m, ["header", "units", "length"], "GeV"),
        changed(m, ["header", "units", "material_amount"], "mol/mm^3"),
        changed(m, ["header"], {}),
        changed(m, ["header"], None),
        changed(m, ["header", "description"], 42),
        changed(m, ["description"], "misplaced"),
        changed(m, ["surfaces", 0, "target", "geometry_id"], 72057594037927936),
        changed(m, ["surfaces", 0, "target", "geometry_id"], "72057594037927936"),
        changed(m, ["surfaces", 0, "target", "geometry_id"], {"volume": -1}),
        changed(m, ["surfaces", 0, "target", "geometry_id"], {"volume": 1.5}),
        changed(m, ["surfaces", 0, "target", "geometry_id"], {"approach": 1}),
        changed(m, ["surfaces", 0, "target", "geometry_id"], {"boundary": 1}),
        changed(m, ["surfaces", 0, "material", "kind"], "unknown"),
        changed(m, ["surfaces", 0, "material", "slab", "thickness"], -1),
        changed(
            s,
            ["surfaces", 3, "material", "axes"],
            [{"kind": "equidistant", "bins": 2}] * 2,
        ),
        changed(s, ["surfaces", 5, "material"], None),
        changed(m, ["volumes"], []),
        changed(m, ["volumes"], [{"geometry_id": {"volume": 1}, "material": None}]),
        changed(m, ["typo"], True),
    ]
    for name, maximum in maxima.items():
        structural.append(
            changed(m, ["surfaces", 0, "target", "geometry_id"], {name: maximum + 1})
        )
    for document in structural:
        assert not validator.is_valid(document), document


def test_unit_symbols():
    validator, examples = validation_inputs()
    for dimension, definition in validator.schema["properties"]["header"]["properties"][
        "units"
    ]["properties"].items():
        for symbol in definition["enum"]:
            validator.validate(
                changed(examples["minimal"], ["header", "units", dimension], symbol)
            )


if __name__ == "__main__":
    main()
