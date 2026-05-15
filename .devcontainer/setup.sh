#!/bin/bash

WORKSPACE=$1

git config --global --add safe.directory $WORKSPACE

pushd $WORKSPACE
pre-commit install
# Signing does not seem to work
git config --add commit.gpgsign false
popd

cat >> $HOME/.local/bin/configure_acts <<EOF
echo "+ cmake -S $WORKSPACE -B $WORKSPACE/build_devcontainer --preset github-ci"
cmake -S $WORKSPACE -B $WORKSPACE/build_devcontainer --preset github-ci
EOF
chmod +x $HOME/.local/bin/configure_acts

cat >> $HOME/.local/bin/build_acts <<EOF
echo "+ cmake --build $WORKSPACE/build_devcontainer \$@"
cmake --build $WORKSPACE/build_devcontainer \$@
EOF
chmod +x $HOME/.local/bin/build_acts

cat > /etc/motd <<EOF
=============== ACTS dev container with dependencies ===============
- If you want to run with ODD, you'll need to run:
    git submodule init && git submodule update
- Configure:
    configure_acts
- Build:
    build_acts
  OR
    cmake --build build_devcontainer
- Run:
    source build_devcontainer/this_acts_withdeps.sh
    acts/Examples/Scripts/Python/full_chain_odd.py -n1
EOF
