import pytest

import acts


def test_navigation_policy_factory():

    policy = (
        acts.NavigationPolicyFactory.make()
        .add(acts.TryAllNavigationPolicy)
        .add(
            acts.SurfaceArrayNavigationPolicy,
            acts.SurfaceArrayNavigationPolicy.Config(
                layerType=acts.SurfaceArrayNavigationPolicy.LayerType.Disc,
                bins=(10, 10),
            ),
        )
    )

    policy._buildTest()

    policy = acts.NavigationPolicyFactory.make().add(acts.TryAllNavigationPolicy)

    policy._buildTest()


def test_navigation_policy_factory_build_empty():
    policy = acts.NavigationPolicyFactory.make()

    with pytest.raises(RuntimeError):
        policy._buildTest()


def test_try_all_arguments():
    acts.NavigationPolicyFactory.make().add(
        acts.TryAllNavigationPolicy, acts.TryAllNavigationPolicy.Config(sensitives=True)
    )


def test_detray_grid_export_policy():
    cfg = acts.DetrayGridExportPolicy.Config(
        vertexZHalfRange=100 * acts.UnitConstants.mm, fixedWindow=True
    )
    assert cfg.vertexZHalfRange == 100 * acts.UnitConstants.mm
    assert cfg.fixedWindow

    surfaceArrayConfig = acts.SurfaceArrayNavigationPolicy.Config(
        layerType=acts.SurfaceArrayNavigationPolicy.LayerType.Disc,
        bins=(10, 10),
    )

    # Finds the earlier surface array sibling
    acts.NavigationPolicyFactory.make().add(acts.TryAllNavigationPolicy).add(
        acts.SurfaceArrayNavigationPolicy, surfaceArrayConfig
    ).add(acts.DetrayGridExportPolicy, cfg)._buildTest()

    # Needs a surface array sibling before it
    with pytest.raises(ValueError):
        acts.NavigationPolicyFactory.make().add(acts.TryAllNavigationPolicy).add(
            acts.DetrayGridExportPolicy, cfg
        )._buildTest()

    with pytest.raises(ValueError):
        acts.NavigationPolicyFactory.make().add(acts.DetrayGridExportPolicy, cfg).add(
            acts.SurfaceArrayNavigationPolicy, surfaceArrayConfig
        )._buildTest()
