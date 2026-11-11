import gc
import weakref

import numpy as np
import pytest

import acts

alignment = pytest.importorskip("acts.examples.alignment")


def test_structure_without_frame():
    structure = alignment.AlignableStructure(acts.GeometryIdentifier(volume=1))
    assert not structure.hasTransform
    with pytest.raises(RuntimeError, match="no frame provider"):
        structure.localToGlobalTransform(
            acts.GeometryContext.dangerouslyDefaultConstruct()
        )
    with pytest.raises(ValueError, match="provider is empty"):
        alignment.AlignableStructure(
            acts.GeometryIdentifier(volume=1), transformProvider=None
        )


def test_structure_fixed_frame():
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    frame = acts.Transform3(
        acts.Vector3(12.0, -8.0, 30.0),
        acts.RotationMatrix3(*(acts.Vector3(*rotation[:, i]) for i in range(3))),
    )
    structure = alignment.AlignableStructure(
        acts.GeometryIdentifier(volume=1), transform=frame
    )
    del frame
    gc.collect()
    assert structure.hasTransform
    resolved = structure.localToGlobalTransform(
        acts.GeometryContext.dangerouslyDefaultConstruct()
    )
    np.testing.assert_allclose(resolved.translation, [12.0, -8.0, 30.0])
    np.testing.assert_allclose(
        [[resolved.rotation[i, j] for j in range(3)] for i in range(3)], rotation
    )


def test_structure_owns_frame_provider():
    class Provider:
        def __init__(self):
            self.frame = acts.Transform3(acts.Vector3(1.0, 2.0, 3.0))

        def __call__(self, gctx):
            assert isinstance(gctx, acts.GeometryContext)
            return self.frame

    provider = Provider()
    retained = weakref.ref(provider)
    structure = alignment.AlignableStructure(
        acts.GeometryIdentifier(volume=1), transformProvider=provider
    )
    del provider
    gc.collect()
    assert retained() is not None
    context = acts.GeometryContext.dangerouslyDefaultConstruct()
    first = structure.localToGlobalTransform(context)
    retained().frame = acts.Transform3(acts.Vector3(4.0, 5.0, 6.0))
    second = structure.localToGlobalTransform(context)
    np.testing.assert_allclose(first.translation, [1.0, 2.0, 3.0])
    np.testing.assert_allclose(second.translation, [4.0, 5.0, 6.0])
    del structure
    gc.collect()
    assert retained() is None


@pytest.mark.parametrize("initial_kind", ["missing", "fixed"])
def test_structure_set_transform(initial_kind):
    context = acts.GeometryContext.dangerouslyDefaultConstruct()
    structure = alignment.AlignableStructure(acts.GeometryIdentifier(volume=1))
    fixed_translation = [12.0, -8.0, 30.0]
    provider_translation = [-4.0, 5.0, 6.0]

    with pytest.raises(ValueError, match="provider is empty"):
        structure.setTransform(transformProvider=None)
    assert not structure.hasTransform

    if initial_kind == "fixed":
        frame = acts.Transform3(acts.Vector3(*fixed_translation))
        structure.setTransform(transform=frame)
        del frame
        gc.collect()
        assert structure.hasTransform
        np.testing.assert_allclose(
            structure.localToGlobalTransform(context).translation, fixed_translation
        )

    class Provider:
        def __call__(self, gctx):
            assert isinstance(gctx, acts.GeometryContext)
            return acts.Transform3(acts.Vector3(*provider_translation))

    provider = Provider()
    retained = weakref.ref(provider)
    structure.setTransform(transformProvider=provider)
    del provider
    gc.collect()
    assert retained() is not None
    assert structure.hasTransform
    np.testing.assert_allclose(
        structure.localToGlobalTransform(context).translation, provider_translation
    )

    with pytest.raises(ValueError, match="provider is empty"):
        structure.setTransform(transformProvider=None)
    assert retained() is not None
    np.testing.assert_allclose(
        structure.localToGlobalTransform(context).translation, provider_translation
    )

    structure.setTransform(transform=acts.Transform3(acts.Vector3(*fixed_translation)))
    gc.collect()
    assert retained() is None
    with pytest.raises(ValueError, match="provider is empty"):
        structure.setTransform(transformProvider=None)
    assert structure.hasTransform
    np.testing.assert_allclose(
        structure.localToGlobalTransform(context).translation, fixed_translation
    )
