import pytest

from helpers import alignmentEnabled


@pytest.mark.skipif(not alignmentEnabled, reason="Alignment module is not set up")
def test_alignment_mask_binding():
    from acts.examples.alignment import AlignmentMask

    combined = AlignmentMask.Center0 | AlignmentMask.Rotation1
    assert combined.value == AlignmentMask.Center0.value | AlignmentMask.Rotation1.value


@pytest.mark.skipif(not alignmentEnabled, reason="Alignment module is not set up")
def test_iteration_state_rejected_assignment_preserves_config():
    from acts.examples.alignment import AlignmentAlgorithmConfig, AlignmentMask

    cfg = AlignmentAlgorithmConfig()
    cfg.iterationState = {
        0: AlignmentMask.Center0,
        1: AlignmentMask.Rotation1.value,
    }
    original = dict(cfg.iterationState)

    with pytest.raises(TypeError):
        cfg.iterationState = {0: AlignmentMask.All, 2: "bad"}
    assert dict(cfg.iterationState) == original

    with pytest.raises(TypeError):
        cfg.iterationState = {"nope": AlignmentMask.Center0}
    assert dict(cfg.iterationState) == original

    cfg.iterationState = {3: AlignmentMask.All}
    # iterationState stores the six alignment degrees of freedom. All is 0xFF.
    stored_dofs = (
        AlignmentMask.Center0
        | AlignmentMask.Center1
        | AlignmentMask.Center2
        | AlignmentMask.Rotation0
        | AlignmentMask.Rotation1
        | AlignmentMask.Rotation2
    )
    assert cfg.iterationState == {3: stored_dofs.value}
