PLEASE FOLLOW THE CHECKLIST BELOW WHEN CREATING A NEW PULL REQUEST. THE
CHECKLIST IS FOR YOUR INFORMATION AN MUST BE REMOVED BEFORE SUBMITTING THE PULL
REQUEST.

## Checklist

- [ ] Does the PR title follow the *Conventional Commits* guideline?
- [ ] Is the PR assigned to a milestone?
- [ ] Does this PR close an existing issue?

    - [ ] If yes: is the issue correctly linked so it will be automatically
        closed upon successful merge (See closing keywords link in the sidebar)?
    - [ ] If no: why not?

- [ ] Does this modify the public API as defined in `docs/versioning.rst`?

    - [ ] Does the PR title contain a `!` to indicate a breaking change using
          *Conventional Commits*?
    - [ ] Is there an entry that describes the user-facing changes in the
          public changelog in `docs/changelog.rst`?

- [ ] Does this modify a user-facing API that is not part of the public API as
    defined in `docs/versioning.rst`?

    - [ ] Is there an entry that describes the user-facing changes in the
        public changelog in `docs/changelog.rst`?
