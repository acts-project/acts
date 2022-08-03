#!/bin/bash

find . -type f -name "*.hpp" -print0 | parallel -0 clang-tidy -p build {} 2>/dev/null
