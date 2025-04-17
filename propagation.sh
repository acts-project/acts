#!/bin/bash
set -e

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}"   )" &> /dev/null && pwd   )

cd $SCRIPT_DIR

source $SCRIPT_DIR/build/this_acts_withdeps.sh
source $SCRIPT_DIR/.venv/bin/activate

exec python3 Examples/Scripts/Python/propagation.py --geo $1
