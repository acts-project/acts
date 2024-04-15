import pathlib, acts, acts.examples
from acts import covfie_conversion as cc


def test_constant_field_conversion():
    v = acts.Vector3(1, 2, 3)
    af = acts.ConstantBField(v)
    cf = cc.covfieField(af)
    view = cc.newView(cf)
    for x, y, z in [(0, 0, 1), (1, 1, 1), (1, 0, 2)]:
        assert view.at(x, y, z) == [1, 2, 3]


def test_root_field_conversion():
    current_file_path = pathlib.Path(__file__).resolve().parent
    p = (
        current_file_path.parent.parent.parent
        / "thirdparty"
        / "OpenDataDetector"
        / "data"
        / "odd-bfield.root"
    )

    af = acts.examples.MagneticFieldMapXyz(str(p))
    cf = cc.covfieField(af)
    view = cc.newView(cf)
    points = [
        (9300.0, 4700.0, 11200.0),
        (10000.0, 10000.0, 14300.0),
        (-2900.0, -4700.0, 5200.0),
        (-2900.0, -4800.0, 9100.0),
        (-2900.0, -5200.0, -8800.0),
        (-4400.0, 4800.0, -12700.0),
        (-6600.0, 1900.0, 7700.0),
        (-9700.0, -900.0, 12700.0),
        (-10000.0, -10000.0, -13000.0),
        (10000.0, 10000.0, 14300.0),
    ]

    error_margin_half_width = 0.0001
    for x, y, z in points:
        val = af.getFieldUnchecked(acts.Vector3(x, y, z))
        Bx1, By1, Bz1 = val[0], val[1], val[2]

        Bx2, By2, Bz2 = tuple(view.at(x, y, z))

        assert (
            abs(Bx1 - Bx2) < error_margin_half_width
            and abs(By1 - By2) < error_margin_half_width
            and abs(Bz1 - Bz2) < error_margin_half_width
        )
