# Track Finding

:::{todo}
Write CKF documentation
:::

(ckf_core)=
## Combinatorial Kalman Filter

% Functional approach
% Implementation as an actor
% Limitations (no inward filtering)
% Requirement to copy out after smoothing if smoothed track states are desired

## Machine-Learning based Track Finding

There is a lot of research ongoing about machine-learning based approaches to Track Finding. Because these are not yet stable and bullet-prove, there are no such algorithms distributed with the core library. However, there exists a [plugin](exatrkxplugin), that implements the *Exa.TrkX* algorithm in ACTS.
