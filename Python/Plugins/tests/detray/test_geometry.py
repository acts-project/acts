from acts import detray


def test_barcode():
    barcode = detray.geometry.barcode()
    assert barcode.is_invalid
    barcode.volume = 1
    barcode.id = detray.surface_id.e_sensitive
    barcode.index = 3
    barcode.transform = 4
    barcode.extra = 5
    assert not barcode.is_invalid
    assert barcode.volume == 1
    assert barcode.id == detray.surface_id.e_sensitive
    assert barcode.index == 3
    assert barcode.transform == 4
    assert barcode.extra == 5


def test_dtyped_index():
    didx = detray.dtyped_index()
    didx.id = 1
    didx.index = 2
    assert didx.id == 1
    assert didx.index == 2
    # i don't know how i produce an invalid value
    # assert not didx.is_invalid
    # assert not didx.is_invalid_id
    # assert not didx.is_invalid_index
    didx.shift(1)
    assert didx.index == 3
    # assert didx.is_invalid
    # assert didx.is_invalid_id
    # assert didx.is_invalid_index
    assert str(didx) == "[id = 1(1) | index = 3]"
