#!/bin/bash

find . -type f \( -name "*.hpp" -or -name "*.ipp" \) -and -not -path "*build*" -and -not -path "*install*" -print0 | parallel -0 clang-tidy -p build {} 2>/dev/null
