#!/usr/bin/env python3

from physmon_traccc_common import runTracccChain

# With Traccc seeding.
runTracccChain("host", False)

# With Acts seeding.
runTracccChain("host", True)