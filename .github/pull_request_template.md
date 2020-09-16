PLEASE FOLLOW THE CHECKLIST BELOW WHEN CREATING A NEW PULL REQUEST. THE
CHECKLIST IS FOR YOUR INFORMATION AND MUST BE REMOVED BEFORE SUBMITTING THE PULL
REQUEST.

## Checklist

- [ ] Does the PR title follow the `<prefix>: title` scheme?

    The prefix must be one of:

    - `fix`: for a bugfix
    - `feat`: for a new feature
    - `refactor`: for an improvement of an existing feature
    - `perf`, `test`: for performance- or test-related changes
    - `docs`: for documentation-related changes
    - `build`, `ci`, `chore`: as appropriated for infrastructure changes

- [ ] Does this modify the public API as defined in `docs/versioning.rst`?

    - [ ] Does the PR title contain a `!` to indicate a breaking change?
    - [ ] Is there section starting with `BREAKING CHANGE:` in the PR body
        that explains the breaking change?

- [ ] Is the PR ready to be merged?

    - [ ] If not: is it marked as work-in-progress using the `WIP` label?

- [ ] Is the PR assigned to a milestone? This should be `next` unless you
    target a specific release.

- [ ] Does this PR close an existing issue?

    - [ ] Is the issue correctly linked so it will be automatically closed
        upon successful merge (See closing keywords link in the sidebar)?
    - [ ] Does the PR milestone match the issue milestone?
