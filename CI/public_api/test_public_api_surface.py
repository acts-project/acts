#!/usr/bin/env python3
"""Self-tests for the public API surface tooling.

Guards two things against silent breakage:
  * the pure classification/extraction logic (fast, no external deps), and
  * the end-to-end pipeline over real Doxygen output, using committed fixture
    headers with a known set of API differences -- this also catches a Doxygen
    upgrade quietly changing its XML.

Run directly (``python3 CI/public_api/test_public_api_surface.py``): all ``test_*``
functions run; a non-zero exit means a failure. The e2e test is skipped when
``doxygen`` is not on PATH. Also collectable by pytest.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import public_api_surface as pas  # noqa: E402
import public_api_diff as pad  # noqa: E402

# --- pure unit tests -------------------------------------------------------


def test_callable_forms_expands_defaulted_args():
    md = ET.fromstring(
        '<memberdef kind="function" const="no">'
        "<type>void</type><name>foo</name>"
        "<param><type>int</type><declname>x</declname></param>"
        "<param><type>double</type><declname>y</declname><defval>1.0</defval></param>"
        "</memberdef>"
    )
    forms = pas.callable_forms(md, "Acts::foo")
    assert "Acts::foo(int)" in forms, forms  # shorter call still valid
    assert "Acts::foo(int, double)" in forms, forms


def test_callable_forms_const_qualifier():
    md = ET.fromstring(
        '<memberdef kind="function" const="yes">'
        "<type>int</type><name>size</name></memberdef>"
    )
    assert "Acts::C::size() const" in pas.callable_forms(md, "Acts::C::size")


def test_module_of_components():
    assert pas.module_of("/w/Core/include/Acts/Surfaces/Plane.hpp") == "Core/Surfaces"
    assert (
        pas.module_of("/w/Plugins/Json/include/ActsPlugins/Json/X.hpp") == "Plugin:Json"
    )
    assert pas.module_of("/w/Fatras/include/ActsFatras/Y.hpp") == "Fatras"
    assert pas.module_of("/w/Alignment/include/ActsAlignment/Z.hpp") == "Alignment"


def test_sanitize_strips_bidi_and_zero_width_chars():
    # U+202E right-to-left override, U+200B zero-width space: neither is
    # visible whitespace, so a crafted C++20 Unicode identifier could smuggle
    # either into a report without this stripping them.
    assert pas._sanitize("Acts::Foo‮​") == "Acts::Foo"


def test_name_and_norm_type_sanitize_extracted_text():
    md = ET.fromstring(
        "<memberdef><name>size‮</name>" "<type>Acts::Bar​*</type></memberdef>"
    )
    assert pas._name(md, "name") == "size"
    assert pas._norm_type(md.find("type")) == "Acts::Bar*"


def test_classify_sanitizes_untrusted_json_input():
    # public_api_diff.py is the trusted (base-branch) script, but its input
    # JSON comes from an unprivileged job a PR fully controls -- it could
    # hand-craft this JSON directly, bidi/zero-width characters included.
    base = {"symbols": [], "callables": {}, "fields": {}}
    head = {"symbols": ["type Acts::New‮​"], "callables": {}, "fields": {}}
    c = pad.classify(base, head)
    assert c["added_names"] == ["type Acts::New"]


def test_classify_additions_and_breaking():
    base = {
        "symbols": ["type Acts::Old", "type Acts::Keep"],
        "callables": {
            "Acts::foo(int)": "void",
            "Acts::gone(double)": "void",
            "Acts::baz()": "int",
        },
        "fields": {"Acts::S::a": "double", "Acts::S::gone": "int"},
    }
    head = {
        "symbols": ["type Acts::Keep", "type Acts::New"],
        "callables": {
            "Acts::foo(int)": "void",
            "Acts::foo(int, double)": "void",
            "Acts::bar(int, double)": "void",
            "Acts::baz()": "long",
        },
        "fields": {"Acts::S::a": "float", "Acts::S::added": "bool"},
    }
    c = pad.classify(base, head)
    # additions
    assert "type Acts::New" in c["added_names"]
    assert "Acts::foo(int, double)" in c["added_signatures"]  # defaulted-arg add
    assert "Acts::S::added" in c["added_fields"]
    # breaking
    assert "type Acts::Old" in c["removed_names"]
    assert "Acts::gone(double)" in c["removed_signatures"]
    assert any("Acts::baz()" in s for s in c["return_type_changes"])
    assert "Acts::S::gone" in c["removed_fields"]
    assert any("Acts::S::a" in s for s in c["field_type_changes"])
    assert c["has_additions"] and c["has_breaking"]


def test_classify_defaulted_arg_is_not_breaking():
    base = {"symbols": [], "callables": {"Acts::f(int)": "void"}, "fields": {}}
    head = {
        "symbols": [],
        "callables": {"Acts::f(int)": "void", "Acts::f(int, double)": "void"},
        "fields": {},
    }
    c = pad.classify(base, head)
    assert c["has_additions"] and not c["has_breaking"], c


def test_classify_no_change():
    snap = {
        "symbols": ["type Acts::A"],
        "callables": {"Acts::f()": "void"},
        "fields": {"Acts::A::x": "int"},
    }
    c = pad.classify(snap, snap)
    assert not c["has_additions"] and not c["has_breaking"]


# --- end-to-end fixture test ----------------------------------------------


def _measure(input_dir: Path) -> dict:
    out = Path(tempfile.mkdtemp(prefix="api-selftest-"))
    pas.run_doxygen(repo=HERE.parents[1], out=out, input_dirs=[input_dir])
    return pas.parse_xml(out / "xml")


def test_end_to_end_fixture():
    if not shutil.which("doxygen"):
        print("  (skipped: doxygen not on PATH)")
        return
    data = HERE / "testdata"
    base = _measure(data / "base" / "Acts")
    head = _measure(data / "head" / "Acts")
    c = pad.classify(base, head)

    # additions
    assert "type Acts::NewThing" in c["added_names"], c
    assert "Acts::Demo::flag" in c["added_fields"], c
    assert any("Acts::Demo::doThing(int, int)" == s for s in c["added_signatures"]), c
    assert any(
        s.startswith("Acts::compute(int, double)") for s in c["added_signatures"]
    ), c
    # breaking
    assert any(s == "Acts::Demo::doThing(int)" for s in c["removed_signatures"]), c
    assert any(
        s.startswith("Acts::Demo::oldFn(double)") for s in c["removed_signatures"]
    ), c
    assert any("Acts::Demo::tolerance" in s for s in c["field_type_changes"]), c
    # compute(int) must survive (defaulted arg) -> not a removal
    assert not any(s == "Acts::compute(int)" for s in c["removed_signatures"]), c


def test_detray_roots_and_gpu_api():
    with tempfile.TemporaryDirectory() as tmp:
        repo = Path(tmp)
        public = repo / "Detray/core/include"
        public.mkdir(parents=True)
        for rel in (
            "Detray/tests/include",
            "Detray/tutorials/include",
            "Detray/extern/vecmem/include",
            "Detray/plugins/algebra/array/include",
        ):
            (repo / rel).mkdir(parents=True)
        roots = pas.standard_roots(repo)
        assert set(roots) == {public, repo / "Detray/plugins/algebra/array/include"}
        assert pas.module_of(str(public / "detray/Thing.hpp")) == "Detray/core"
        if not shutil.which("doxygen"):
            return
        header = public / "Thing.hpp"
        header.write_text("""
namespace detray {
/// Public device-compatible type.
struct Thing {
  /// A shared operation.
  DETRAY_HOST_DEVICE int run(int value) const;
  /// A device operation.
  DETRAY_DEVICE int device(int value);
  /// A host operation.
  DETRAY_HOST int host(int value);
};
namespace detail {
/// Internal type.
struct Hidden {};
}
}
""")
        base = _measure(public)
        assert "type detray::Thing" in base["symbols"]
        assert "type detray::detail::Hidden" not in base["symbols"]
        header.write_text(
            header.read_text().replace("run(int value)", "run(double value)")
        )
        head = _measure(public)
        changes = pad.classify(base, head)
        assert "detray::Thing::run(int) const" in changes["removed_signatures"], changes
        assert (
            "detray::Thing::run(double) const" in changes["added_signatures"]
        ), changes


def test_detray_validation_roots():
    with tempfile.TemporaryDirectory() as tmp:
        repo = Path(tmp)
        included = {
            "Detray/tests/include/detray/test/common",
            "Detray/tests/include/detray/test/utils",
            "Detray/tests/include/detray/test/validation",
            "Detray/tests/tools/include",
        }
        excluded = {
            "Detray/tests/include/detray/test/framework",
            "Detray/tests/include/detray/test/cpu",
            "Detray/tests/include/detray/test/device",
            "Detray/tests/unit_tests",
            "Detray/tests/integration_tests",
            "Detray/tests/benchmarks/include",
            "Detray/tests/tools/src",
        }
        for rel in included | excluded:
            (repo / rel).mkdir(parents=True)
        config = (
            repo / "Detray/tests/include/detray/test/framework/test_configuration.hpp"
        )
        config.touch()
        roots = pas.standard_roots(repo)
        assert set(roots) == {repo / rel for rel in included} | {config}
        assert pas.module_of(str(config)) == "Detray/tests"
        assert (
            pas.module_of(str(repo / "Detray/tests/tools/include/detray/options/x.hpp"))
            == "Detray/tests/tools"
        )


def test_detray_generated_declarations_are_types():
    if not shutil.which("doxygen"):
        return
    with tempfile.TemporaryDirectory() as tmp:
        headers = Path(tmp)
        (headers / "declarations.hpp").write_text(
            "#define DETRAY_DECLARE_MATERIAL(NAME, VALUE) struct NAME {};\n"
            "#define DETRAY_DECLARE_MATERIAL_WITH_DED(NAME, VALUE) struct NAME {};\n"
            "#define DETRAY_DECLARE_PARTICLE(NAME, PDG, MASS, CHARGE) struct NAME {};\n"
        )
        (headers / "materials.hpp").write_text(
            '#include "declarations.hpp"\n'
            "namespace detray {\n"
            "/// A predefined material.\n"
            "DETRAY_DECLARE_MATERIAL(silicon, 14);\n"
            "/// A material with density-effect data.\n"
            "DETRAY_DECLARE_MATERIAL_WITH_DED(copper, 29);\n"
            "/// A particle hypothesis.\n"
            "DETRAY_DECLARE_PARTICLE(muon, 13, 105.66, -1);\n"
            "}\n"
        )
        surface = _measure(headers)
        assert "type detray::silicon" in surface["symbols"]
        assert "type detray::copper" in surface["symbols"]
        assert "type detray::muon" in surface["symbols"]
        assert not any("DETRAY_DECLARE" in s for s in surface["callables"])


def main() -> int:
    tests = [
        v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)
    ]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS {t.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL {t.__name__}: {e}")
        except Exception as e:  # noqa: BLE001
            failed += 1
            print(f"ERROR {t.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
