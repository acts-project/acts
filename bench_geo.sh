#!/bin/bash

N=1000
exec hyperfine -w1 \
  -n gen1 "Examples/Scripts/Python/propagation.py --geo gen1 -n $N" \
  -n gen3 "Examples/Scripts/Python/propagation.py --geo gen3 -n $N" \
