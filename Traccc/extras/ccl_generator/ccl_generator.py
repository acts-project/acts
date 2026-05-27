#!/usr/bin/env python3


import argparse
import numpy
import math
import scipy.stats
import random
import collections
import csv


def neighbourhood(p, m):
    return [
        (x, y)
        for (x, y) in [
            (p[0] + 1, p[1]),
            (p[0] - 1, p[1]),
            (p[0] + 1, p[1] - 1),
            (p[0] - 1, p[1] - 1),
            (p[0] + 1, p[1] + 1),
            (p[0] - 1, p[1] + 1),
            (p[0], p[1] - 1),
            (p[0], p[1] + 1),
        ]
        if x >= 0 and y >= 0 and x < m and y < m
    ]


def rand(d):
    return int(round(d.rvs()))


def generate_file(name, Md, Hd, size, N):
    print("Generating file %s..." % name)

    with open(name, "w") as f:
        h = 0
        w = csv.DictWriter(
            f,
            fieldnames=[
                "geometry_id",
                "measurement_id",
                "channel0",
                "channel1",
                "timestamp",
                "value",
            ],
            lineterminator="\n",
        )
        w.writeheader()

        for m in range(N):
            acts = 0
            hits = rand(Md)

            points = set()

            for q in range(hits):
                cells = rand(Hd)
                center = (random.randint(0, size - 1), random.randint(0, size - 1))

                cands = [center]
                seen = set()

                for c in range(cells):
                    if not cands:
                        break

                    p = random.choice(cands)

                    if p not in points:
                        w.writerow(
                            {
                                "geometry_id": m,
                                "measurement_id": h,
                                "channel0": p[0],
                                "channel1": p[1],
                                "timestamp": 0,
                                "value": random.uniform(0.1, 1.0),
                            }
                        )

                    acts += 1

                    seen.add(p)
                    cands += neighbourhood(p, size)
                    cands = [x for x in cands if x not in seen and x not in points]

                h += 1
                points |= seen


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate some CCL examples")

    parser.add_argument(
        "-C",
        type=int,
        help="number of files to generate",
    )
    parser.add_argument(
        "-N",
        type=int,
        default=2500,
        help="number of modules in the detector",
    )
    parser.add_argument(
        "-S",
        type=int,
        default=655,
        help="size of each module (SxS)",
    )

    parser.add_argument(
        "--Mm",
        type=float,
        default=1.65,
        help="mean hits per module",
    )
    parser.add_argument(
        "--Ms",
        type=float,
        default=0.95,
        help="scale of hits per module",
    )

    parser.add_argument(
        "--Hm",
        type=float,
        default=1.80,
        help="mean cells per hit",
    )
    parser.add_argument(
        "--Hs",
        type=float,
        default=0.89,
        help="scale of cells per hit",
    )

    parser.add_argument(
        "-o",
        type=str,
        default="ccl",
        help="output name",
    )

    args = parser.parse_args()

    hits_dist = scipy.stats.lognorm(args.Ms, scale=math.exp(args.Mm))
    cell_dist = scipy.stats.lognorm(args.Hs, scale=math.exp(args.Hm))

    print("Number of modules:               %d" % args.N)
    print("Module size:                     %dx%d" % (args.S, args.S))
    print("Total number of pixels:          %d" % (args.N * args.S * args.S))
    print()
    print(
        "Dist parameters hits per module: %s(μ=%.2f, σ=%.2f)"
        % (hits_dist.dist.name, args.Mm, args.Ms)
    )
    print(
        "Dist parameters cells per hit:   %s(μ=%.2f, σ=%.2f)"
        % (cell_dist.dist.name, args.Hm, args.Hs)
    )
    print()
    print("Expected hits per module:        %.2f" % hits_dist.mean())
    print("Expected total hits:             %.0f" % (args.N * hits_dist.mean()))
    print("Expected cells per hit:          %.2f" % cell_dist.mean())
    print(
        "Expected cells per module:       %.0f" % (cell_dist.mean() * hits_dist.mean())
    )
    print(
        "Expected total cells:            %.0f"
        % (args.N * cell_dist.mean() * hits_dist.mean())
    )
    print(
        "Expected density:                %.5f%%"
        % ((100 * hits_dist.mean() * cell_dist.mean()) / (args.S * args.S))
    )
    print()

    if args.C is None:
        generate_file(("%s.csv" % args.o), hits_dist, cell_dist, args.S, args.N)
    else:
        for i in range(args.C):
            generate_file(
                "%s_%010d.csv" % (args.o, i), hits_dist, cell_dist, args.S, args.N
            )
