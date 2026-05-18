#!/usr/bin/env python3

import pytest
import tempfile
import pathlib
import os

# Test imports to verify availability - make optional dependencies optional
try:
    pytest.importorskip("acts.ActsPythonBindingsDD4hep")
    HAS_DD4HEP = True
except pytest.skip.Exception:
    HAS_DD4HEP = False

try:
    pytest.importorskip("acts.ActsPythonBindingsGeoModel")
    HAS_GEOMODEL = True
except pytest.skip.Exception:
    HAS_GEOMODEL = False

import acts
import acts.examples

# Optional imports
if HAS_DD4HEP:
    try:
        import acts.examples.dd4hep
        import acts.examples.geant4
        import acts.examples.geant4.dd4hep
    except ImportError:
        pass

if HAS_GEOMODEL:
    try:
        import acts.examples.geomodel
    except ImportError:
        pass

u = acts.UnitConstants


def test_toroidal_field_basic():
    """Test basic ToroidField functionality."""

    # Test default configuration
    config = acts.ToroidField.Config()
    field = acts.ToroidField(config)
    assert field is not None

    # Test field at origin
    ctx = acts.MagneticFieldContext()
    cache = field.makeCache(ctx)

    # Test field calculation at a point in barrel region
    position = acts.Vector3(6000.0, 0.0, 0.0)  # 6m radius
    field_value = field.getField(position, cache)

    # Should have non-zero field components
    assert (
        abs(field_value[0]) > 0.0
        or abs(field_value[1]) > 0.0
        or abs(field_value[2]) > 0.0
    )


def test_toroidal_field_custom():
    """Test ToroidField with custom parameters."""

    config = acts.ToroidField.Config()

    # Customize barrel configuration
    config.barrel.R_in = 5.0
    config.barrel.R_out = 9.0
    config.barrel.I = 15000.0

    field = acts.ToroidField(config)
    assert field is not None

    # Test field calculation
    ctx = acts.MagneticFieldContext()
    cache = field.makeCache(ctx)

    position = acts.Vector3(7000.0, 0.0, 0.0)
    field_value = field.getField(position, cache)

    # Verify field is calculated (should be a Vector3 with components)
    assert hasattr(field_value, "__getitem__")  # Can access components
    assert field_value[0] is not None
    assert field_value[1] is not None
    assert field_value[2] is not None


def test_toroidal_field_symmetry():
    """Test that the field has expected symmetries."""

    config = acts.ToroidField.Config()
    field = acts.ToroidField(config)
    ctx = acts.MagneticFieldContext()
    cache = field.makeCache(ctx)

    # Test azimuthal symmetry - field should be similar at same radius
    radius = 7000.0
    z = 1000.0

    pos1 = acts.Vector3(radius, 0.0, z)
    pos2 = acts.Vector3(0.0, radius, z)

    field1 = field.getField(pos1, cache)
    field2 = field.getField(pos2, cache)

    # Due to toroidal symmetry, field magnitudes should be similar
    mag1 = (field1[0] ** 2 + field1[1] ** 2 + field1[2] ** 2) ** 0.5
    mag2 = (field2[0] ** 2 + field2[1] ** 2 + field2[2] ** 2) ** 0.5

    assert abs(mag1 - mag2) / max(mag1, mag2, 1e-10) < 0.1  # Within 10%


def test_toroidal_field_regions():
    """Test field behavior in different regions (barrel vs endcap)."""

    config = acts.ToroidField.Config()
    field = acts.ToroidField(config)
    ctx = acts.MagneticFieldContext()
    cache = field.makeCache(ctx)

    # Test in barrel region
    pos_barrel = acts.Vector3(7000.0, 0.0, 1000.0)
    field_barrel = field.getField(pos_barrel, cache)
    mag_barrel = (
        field_barrel[0] ** 2 + field_barrel[1] ** 2 + field_barrel[2] ** 2
    ) ** 0.5

    # Test in endcap region
    pos_endcap = acts.Vector3(2000.0, 0.0, 15000.0)
    field_endcap = field.getField(pos_endcap, cache)
    mag_endcap = (
        field_endcap[0] ** 2 + field_endcap[1] ** 2 + field_endcap[2] ** 2
    ) ** 0.5

    # Both should have some field
    assert mag_barrel > 0.0
    assert mag_endcap > 0.0


def test_toroidal_field_configuration():
    """Test configuration classes."""

    # Test BarrelConfig
    barrel_config = acts.ToroidField.BarrelConfig()
    assert barrel_config.R_in > 0
    assert barrel_config.R_out > barrel_config.R_in
    assert barrel_config.I > 0

    # Test EctConfig
    ect_config = acts.ToroidField.EctConfig()
    assert ect_config.R_in > 0
    assert ect_config.R_out > ect_config.R_in
    assert ect_config.I > 0

    # Test LayoutConfig
    layout_config = acts.ToroidField.LayoutConfig()
    assert layout_config.nCoils > 0
    assert layout_config.nArc > 0
    assert layout_config.nStraight > 0


if __name__ == "__main__":
    # Run basic tests if called directly
    test_toroidal_field_basic()
    test_toroidal_field_custom()
    test_toroidal_field_symmetry()
    test_toroidal_field_regions()
    test_toroidal_field_configuration()
    print("All ToroidField tests passed!")
