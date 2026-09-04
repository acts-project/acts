(howto_format)=

# Source code formatting


## Setting up pre-commit hooks

Acts comes with a [`pre-commit`](https://pre-commit.com/) configuration file to enable pre-commit hooks.
In order to use them, one has to install the `pre-commit` package. At the time of writing, the `pre-commit v4.5.0` works with the Acts hooks. It could be that your system has a `pre-commit` installation with an older version, i.e. some `LCG` views use an old version that is incompatible with Acts. The best way is to have `pre-commit` installed in a venv and use that at time of commit. Here are some instruction on how to set it up and make sure to call it.

- Create a python virtual environment and activate it
```console
python -m venv ~/.venvs/commit-env
source ~/.venvs/precommit-env/bin/activate
```

- Install pre-commit inside the venv and verify the version
```console
pip install --upgrade pip
pip install pre-commit
pre-commit --version
```

- Optional: Setup a wrapper around pre-commit to use a clean `PYTHONPATH`. This avoids
conflicts with external setups, like lcg-views.
```console

# Change script name
mv $VIRTUAL_ENV/bin/pre-commit $VIRTUAL_ENV/bin/pre-commit-real

# Make a wrapper to run with clearn PYTHONPATH
cat > $VIRTUAL_ENV/bin/pre-commit << 'EOF'
#!/bin/bash
# Clean PYTHONPATH only for this process
PYTHONPATH="" exec "$VIRTUAL_ENV/bin/python" -m pre_commit "$@"
EOF

# Make it executable
chmod +x $VIRTUAL_ENV/bin/pre-commit
```

One can then setup an external setup, e.g LCG, for compilation and running. This step can
be avoided if one uses two shells, one for developing and one for committing.

- Go the acts repo and, while inside the venv, install the git hook
```console
pre-commit install
```
This adds a hook `.git/hooks/pre-commit` which automatically calls `pre-commit`

- At this point, `git commit` should pick up the pre-commit and run the formatting tools and apply corrections.
If some files were modified, then just `git add` and `git commit` again to pick the changes.

### Execution
```console
cd acts/
pre-commit run --all-files
```

## C++ formatting: `clang-format`

Code formatting is handled by
[`clang-format`](https://clang.llvm.org/docs/ClangFormat.html). A configuration
file is available in the repository root at `.clang-format` and should be used
to automatically format the code. Many editors / IDEs support `clang-format`
and also format-on-save actions.

The ACTS CI system will automatically check code formatting using the provided
`clang-format` configuration and will notify incompatible formatting.

To simplify this, a script located in `CI/check_format` can be used like:

```console
$ CI/check_format $SOURCE_DIR
```

In some cases, different `clang-format` versions will result in slightly
different outputs. In many cases, this is accepted by the CI. However, it is
recommended to use the same major version of `clang-format` to perform local
formatting. Options to obtain a compatible `clang-format` version
are to use your package manager (e.g. Ubuntu distributions usually offer a set of
versions to install), or to use statically linked binaries from
[here](https://github.com/muttleyxd/clang-tools-static-binaries)[^1] and use them with:

```console
CLANG_FORMAT_BINARY=<path/to/clang-format> CI/check_format $SOURCE_DIR
```

You can also download the required changes by clicking on *Summary* on the top left-hand
portion of the CI job and scrolling down to the bottom of the page (see *Changed*).
However, it is suggested to run the `CI/check_format` locally before committing, to not
clog the shared resources with repeated checks.

## Python formatting

Formatting of the Python source code uses the library
[`black`](https://github.com/psf/black). To run it, you can locally install the
`black` package. You can use `pip` to install it:

```console
$ pip install black
$ black <source>
```

:::{tip}
It is **strongly recommended** to use a [virtual
environment](https://realpython.com/python-virtual-environments-a-primer/) for
this purpose! For example, run

```console
$ python -m venv venv
$ source venv/bin/activate
```

and then install and use black. You can also use a tool like
[`pipx`](https://github.com/pypa/pipx) to simplify this.
:::

[^1]: This repository is external to the ACTS project, so proceed with caution!
