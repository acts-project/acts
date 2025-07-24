import pytest
import acts


def test_navigation_policy_factory():

    policy = (
        acts.MyNavigationPolicyFactory.make()
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

    policy = acts.MyNavigationPolicyFactory.make().add(acts.TryAllNavigationPolicy)

    policy._buildTest()


def test_navigation_policy_factory_build_empty():
    policy = acts.MyNavigationPolicyFactory.make()

    with pytest.raises(RuntimeError):
        policy._buildTest()


def test_navigation_policy_factory_add_multiple():
    with pytest.raises(ValueError):
        (
            acts.MyNavigationPolicyFactory.make()
            .add(acts.TryAllNavigationPolicy)
            .add(acts.TryAllNavigationPolicy)
        )


def test_try_all_arguments():
    acts.MyNavigationPolicyFactory.make().add(
        acts.TryAllNavigationPolicy, acts.TryAllNavigationPolicy.Config(sensitives=True)
    )
