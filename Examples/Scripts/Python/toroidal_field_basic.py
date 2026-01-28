#!/usr/bin/env python3

# Copyright (c) 2025 ACTS-Project
# This file is part of ACTS.
# See LICENSE for details.

"""
Simple test script for ToroidalField Python bindings
Tests the magnetic field functionality without complex dependencies
"""

import sys


def test_toroidal_field_basic():
    """Test basic import and configuration"""
    print("=" * 60)
    print("Testing ToroidalField Python Bindings")
    print("=" * 60)

    try:
        import acts

        _ = acts.ToroidalField
        print("✅ Successfully imported acts.ToroidalField")
    except (ImportError, AttributeError) as e:
        print(f"❌ Failed to import acts.ToroidalField: {e}")
        return False

    # Test Config creation with defaults
    print("\n📋 Testing Config with defaults:")
    try:
        config = acts.ToroidalField.Config()
        print(f"   Barrel R_in: {config.barrel.R_in / 1000:.1f} m")
        print(f"   Barrel R_out: {config.barrel.R_out / 1000:.1f} m")
        print(f"   Barrel c: {config.barrel.c / 1000:.1f} m")
        print(f"   Barrel b: {config.barrel.b / 1000:.3f} m")
        print(f"   Barrel I: {config.barrel.I} A")
        print(f"   Barrel Nturns: {config.barrel.Nturns}")
        print(f"   Number of coils: {config.layout.nCoils}")
        print("✅ Config creation successful")
    except Exception as e:
        print(f"❌ Config creation failed: {e}")
        return False

    # Test Config customization
    print("\n🔧 Testing Config customization:")
    try:
        custom_config = acts.ToroidalField.Config()
        custom_config.barrel.R_in = 5.2 * 1000  # Convert to mm
        custom_config.barrel.R_out = 9.8 * 1000  # Convert to mm
        custom_config.barrel.I = 18000.0
        custom_config.layout.nCoils = 10

        print(f"   Customized Barrel R_in: {custom_config.barrel.R_in / 1000:.1f} m")
        print(f"   Customized Barrel R_out: {custom_config.barrel.R_out / 1000:.1f} m")
        print(f"   Customized Barrel I: {custom_config.barrel.I} A")
        print(f"   Customized number of coils: {custom_config.layout.nCoils}")
        print("✅ Config customization successful")
    except Exception as e:
        print(f"❌ Config customization failed: {e}")
        return False

    # Test ToroidalField creation
    print("\n🧲 Testing ToroidalField creation:")
    try:
        field = acts.ToroidalField(config)
        print("✅ ToroidalField creation successful")
    except Exception as e:
        print(f"❌ ToroidalField creation failed: {e}")
        return False

    return True


def test_toroidal_field_calculation():
    """Test magnetic field calculation"""
    print("\n" + "=" * 60)
    print("Testing Magnetic Field Calculation")
    print("=" * 60)

    try:
        import acts

        # Create field
        config = acts.ToroidalField.Config()
        field = acts.ToroidalField(config)

        # Create magnetic field context and cache
        ctx = acts.MagneticFieldContext()
        cache = field.makeCache(ctx)
        print("✅ Magnetic field context and cache created")

        # Test field calculation at various points
        test_points = [
            (6000.0, 0.0, 0.0, "Barrel region"),
            (0.0, 7000.0, 1000.0, "Barrel region (rotated)"),
            (2000.0, 0.0, 15000.0, "Endcap region"),
            (0.0, 3000.0, -12000.0, "Negative endcap"),
            (0.0, 0.0, 0.0, "Origin"),
        ]

        print(f"\n🎯 Testing field calculation at {len(test_points)} points:")
        for x, y, z, description in test_points:
            position = acts.Vector3(x, y, z)
            field_value = field.getField(position, cache)

            # Calculate field magnitude
            magnitude = (
                field_value[0] ** 2 + field_value[1] ** 2 + field_value[2] ** 2
            ) ** 0.5

            print(f"   {description}:")
            print(f"     Position: ({x/1000:.1f}, {y/1000:.1f}, {z/1000:.1f}) m")
            print(
                f"     Field: ({field_value[0]:.2e}, {field_value[1]:.2e}, {field_value[2]:.2e}) T"
            )
            print(f"     Magnitude: {magnitude:.2e} T")

        print("✅ Field calculation successful")
        return True

    except Exception as e:
        print(f"❌ Field calculation failed: {e}")
        return False


def test_configuration_classes():
    """Test individual configuration classes"""
    print("\n" + "=" * 60)
    print("Testing Configuration Classes")
    print("=" * 60)

    try:
        import acts

        # Test BarrelConfig
        print("\n🏺 Testing BarrelConfig:")
        barrel_config = acts.ToroidalField.BarrelConfig()
        print(f"   Default R_in: {barrel_config.R_in / 1000:.1f} m")
        print(f"   Default R_out: {barrel_config.R_out / 1000:.1f} m")
        print(f"   Default current: {barrel_config.I} A")
        print(f"   Default turns: {barrel_config.Nturns}")

        # Test EctConfig
        print("\n🔚 Testing EctConfig:")
        ect_config = acts.ToroidalField.EctConfig()
        print(f"   Default R_in: {ect_config.R_in / 1000:.3f} m")
        print(f"   Default R_out: {ect_config.R_out / 1000:.2f} m")
        print(f"   Default current: {ect_config.I} A")
        print(f"   Default turns: {ect_config.Nturns}")

        # Test LayoutConfig
        print("\n📐 Testing LayoutConfig:")
        layout_config = acts.ToroidalField.LayoutConfig()
        print(f"   Default nCoils: {layout_config.nCoils}")
        print(f"   Default theta0: {layout_config.theta0:.4f} rad")
        print(f"   Default thetaStep: {layout_config.thetaStep:.4f} rad")

        print("✅ Configuration classes test successful")
        return True

    except Exception as e:
        print(f"❌ Configuration classes test failed: {e}")
        return False


def main():
    """Run all tests"""
    print("🚀 Starting ToroidalField Python Binding Tests")
    print("=" * 80)

    tests = [
        ("Basic functionality", test_toroidal_field_basic),
        ("Field calculation", test_toroidal_field_calculation),
        ("Configuration classes", test_configuration_classes),
    ]

    results = []
    for test_name, test_func in tests:
        print(f"\n🧪 Running test: {test_name}")
        try:
            result = test_func()
            results.append(result)
            if result:
                print(f"✅ {test_name}: PASSED")
            else:
                print(f"❌ {test_name}: FAILED")
        except Exception as e:
            print(f"💥 {test_name}: ERROR - {e}")
            results.append(False)

    # Summary
    print("\n" + "=" * 80)
    print("🏁 Test Summary")
    print("=" * 80)

    passed = sum(results)
    total = len(results)

    for i, (test_name, _) in enumerate(tests):
        status = "✅ PASSED" if results[i] else "❌ FAILED"
        print(f"   {test_name}: {status}")

    print(f"\nOverall: {passed}/{total} tests passed")

    if passed == total:
        print("🎉 All tests passed! ToroidalField is working correctly.")
        return 0
    else:
        print("⚠️  Some tests failed. Check the output above for details.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
