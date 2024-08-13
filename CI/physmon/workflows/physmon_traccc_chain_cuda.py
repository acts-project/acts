#!/usr/bin/env python3

from physmon_traccc_common import runTracccChain

# With Traccc seeding.
runTracccChain("cuda", False)

# With Acts seeding.
runTracccChain("cuda", True)