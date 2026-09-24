#!/usr/bin/env python3
"""Offline draft-schema/fixture QA, not a production material reader.

Requires jsonschema>=4.18. Run from any directory. No network references are used.
Direct execution validates the schema and examples. Run the self-tests with
``python -m pytest CI/check_material_schema.py`` (also requires pytest).
"""

import copy
import json
import math
from pathlib import Path

from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = ROOT / "Plugins/Json/schema/material-map-v1.schema.json"
EXAMPLES = ROOT / "docs/examples/material-map-v1"
ID_COMPONENTS = ("volume", "boundary", "layer", "approach", "sensitive", "extra")


def geometry_identity(components):
    """Normalize omitted zero components before checking exact identity."""
    return tuple(components.get(name, 0) for name in ID_COMPONENTS)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f"duplicate JSON key: {key}")
        result[key] = value
    return result


def load(path):
    return json.loads(path.read_text(), object_pairs_hook=unique_object)


def axis_count(axis):
    if axis["kind"] == "subdivided":
        base = axis_count(axis["base"])
        sub = axis_count(axis["subdivision"])
        return base + sub - 1 if axis["mode"] == "replace" else base * sub
    if axis["kind"] == "equidistant":
        if "range" in axis:
            require(axis["range"][0] < axis["range"][1], "axis range")
        return axis["bins"]
    edges = axis.get("edges", axis.get("normalized_edges"))
    require(all(a < b for a, b in zip(edges, edges[1:])), "axis edge order")
    if axis["kind"] == "deferred-variable":
        require(edges[0] == 0 and edges[-1] == 1, "normalized endpoints")
    return len(edges) - 1


def check_tree(value):
    if isinstance(value, float):
        require(math.isfinite(value), "nonfinite number")
    elif isinstance(value, list):
        for child in value:
            check_tree(child)
    elif isinstance(value, dict):
        for child in value.values():
            check_tree(child)


def check_examples(document):
    """Selected cross-field checks; see the documentation for the full contract."""
    check_tree(document)
    stores = document.get("slab_stores", {})
    identities = set()
    for assignment in document["surfaces"]:
        target = assignment["target"]
        identity = (
            target["kind"],
            (
                target["key"]
                if target["kind"] == "stable-key"
                else geometry_identity(target["geometry_id"])
            ),
        )
        require(identity not in identities, "duplicate surface identity")
        identities.add(identity)
        payload = assignment["material"]
        if payload is None:
            continue
        if target["kind"] == "stable-key" and "material_key" in payload:
            require(payload["material_key"] == target["key"], "proto key mismatch")
        kind = payload["kind"]
        if kind in ("proto", "binned"):
            counts = [axis_count(a) for a in payload["binning"]["axes"]]
            if kind == "binned":
                require(counts, "binned material requires axes")
                require(len(payload["values"]) == math.prod(counts), "binned size")
        elif kind in ("proto-grid", "grid"):
            counts = [axis_count(a) for a in payload["axes"]]
            if kind == "proto-grid":
                continue
            storage = payload["storage"]
            values = storage.get("values", storage.get("indices"))
            require(len(values) == math.prod(n + 2 for n in counts), "grid size")
            if storage["kind"] == "direct":
                continue
            if storage["kind"] == "indexed":
                slabs = storage["slabs"]
            else:
                require(storage["store"] in stores, "missing slab store")
                slabs = stores[storage["store"]]
            require(all(i < len(slabs) for i in values), "slab index out of range")


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
    require(len(examples) == 3, "expected three fixture documents")
    return validator, examples


def main():
    validator, examples = validation_inputs()
    for name, document in examples.items():
        validator.validate(document)
        check_examples(document)
        print(f"Validated {name}.json (schema + selected fixture invariants)")


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
        check_examples(document)
    assert geometry_identity({}) == geometry_identity({"volume": 0})


def test_invalid_structures():
    validator, examples = validation_inputs()
    m, s = (examples[n] for n in ("minimal", "surfaces"))
    maxima = dict(zip(ID_COMPONENTS, (255, 255, 4095, 255, 1048575, 255)))
    structural = [
        changed(m, ["header", "version"], 2),
        changed(m, ["header"], {}),
        changed(m, ["header"], None),
        changed(m, ["header", "description"], 42),
        changed(m, ["description"], "misplaced"),
        changed(m, ["surfaces", 0, "target", "geometry_id"], 72057594037927936),
        changed(m, ["surfaces", 0, "target", "geometry_id"], "72057594037927936"),
        changed(m, ["surfaces", 0, "target", "geometry_id"], {"volume": -1}),
        changed(m, ["surfaces", 0, "target", "geometry_id"], {"volume": 1.5}),
        changed(m, ["surfaces", 0, "target", "geometry_id"], {"passive": 1}),
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


def test_invalid_semantics():
    import pytest

    validator, examples = validation_inputs()
    m, s, t = (examples[n] for n in ("minimal", "surfaces", "templates"))
    # These schema-valid documents need cross-field validation.
    duplicate = changed(
        m, ["surfaces"], [copy.deepcopy(m["surfaces"][0]) for _ in range(2)]
    )
    duplicate["surfaces"][1]["target"]["geometry_id"]["extra"] = 0
    semantic = [
        duplicate,
        changed(m, ["surfaces"], m["surfaces"] * 2),
        changed(s, ["surfaces", 5, "material", "storage", "store"], "missing"),
        changed(s, ["surfaces", 4, "material", "storage", "indices", 4], 2),
        changed(
            s,
            ["surfaces", 3, "material", "storage", "values"],
            s["surfaces"][3]["material"]["storage"]["values"][:-1],
        ),
        changed(t, ["surfaces", 1, "material", "material_key"], "wrong-key"),
        changed(
            t, ["surfaces", 1, "material", "axes", 1, "normalized_edges"], [0.1, 1]
        ),
    ]
    for document in semantic:
        validator.validate(document)
        with pytest.raises(ValueError):
            check_examples(document)


if __name__ == "__main__":
    main()
