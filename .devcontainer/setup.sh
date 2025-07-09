#!/bin/bash

WORKSPACE=$1

git config --global --add safe.directory $WORKSPACE

pushd $WORKSPACE
pre-commit install
popd

echo "export WORKSPACE=$WORKSPACE" >> $HOME/.bashrc

cat >> $HOME/.bashrc <<EOF
function configure_acts() {
    echo "+ cmake -S $WORKSPACE -B $WORKSPACE/build_devcontainer --preset github-ci"
    cmake -S $WORKSPACE -B $WORKSPACE/build_devcontainer --preset github-ci
}
export configure_acts

function build_acts() {
    echo "+ cmake --build $WORKSPACE/build_devcontainer"
    cmake --build $WORKSPACE/build_devcontainer
}
export build_acts
EOF


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
