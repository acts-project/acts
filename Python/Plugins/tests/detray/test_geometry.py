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


def test_toy_metadata():
    print(detray.toy.metadata)

    print(detray.toy.metadata.mask_ids.e_rectangle2)
    detray.toy.mask_link(
        detray.toy.metadata.mask_ids.e_rectangle2,
        detray.toy.mask_link.index_type(),
    )
