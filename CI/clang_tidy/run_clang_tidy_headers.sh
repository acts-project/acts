#!/bin/bash

function execute() {
    clang-tidy --config .clang-tidy-headers "$@" || true
}

export -f execute

find . -type f \( -name "*.hpp" -or -name "*.ipp" \) -and -not -path "*build*" -and -not -path "*install*" -print0 | parallel -0 execute -p build {} 2>/dev/null
