import pytest

import acts

import acts.examples


def test_navigation_policy_factory():

    policy = (
        acts.NavigationPolicyFactory()
        .add(acts.TryAllPortalNavigationPolicy)
        .add(acts.TryAllSurfaceNavigationPolicy)
    )

    policy._buildTest()

    policy = (
        acts.NavigationPolicyFactory()
        .add(acts.TryAllPortalNavigationPolicy)
        .add(acts.TryAllSurfaceNavigationPolicy)
    )

    policy._buildTest()


def test_navigation_policy_factory_build_empty():
    policy = acts.NavigationPolicyFactory()

    with pytest.raises(RuntimeError):
        policy._buildTest()


def test_navigation_policy_factory_add_multiple():
    with pytest.raises(ValueError):
        (
            acts.NavigationPolicyFactory()
            .add(acts.TryAllPortalNavigationPolicy)
            .add(acts.TryAllSurfaceNavigationPolicy)
            .add(acts.TryAllPortalNavigationPolicy)
        )
