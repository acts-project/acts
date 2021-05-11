#!/usr/bin/env python3
import sys

from headwind.spec import CollectorResult, Run, Metric

import pandas as pd

assert len(sys.argv) > 1
input_file = sys.argv[1]


metrics = []

df = pd.read_csv(input_file)
df.max_rss = df.max_rss / 1024 / 1024 / 1024 # GB

mdf = df.sort_values(by="max_rss", ascending=False)
for row in mdf.head(10).itertuples():
    metrics.append(Metric(name="max_rss_"+row.file, value=row.max_rss, unit="GB", group="compile_max_rss"))

tdf = df.sort_values(by="time", ascending=False)
for row in tdf.head(10).itertuples():
    metrics.append(Metric(name="time_"+row.file, value=row.time, unit="seconds", group="compile_time"))

sys.stdout.write(CollectorResult(metrics=metrics).json(indent=2))