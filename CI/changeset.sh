#!/bin/bash

if [ ! -z "$1" ]; then
  # The base ref is not none, so we are a pull request and we can just find
  # the changeset between this branch and main.
  git diff --name-only origin/$1...HEAD > $4
else
  # Otherwise we are a push run, so we will compare the before and after.
  git diff --name-only $2 $3 > $4
fi
