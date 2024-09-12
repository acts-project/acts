import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
import uproot
import seaborn as sns

p = argparse.ArgumentParser()

p.add_argument(
    "-i",
    "--input",
    nargs="+",
    type=str,
    default="",
    help="Input file(s) with propagation summary",
)

p.add_argument(
    "-l",
    "--label",
    nargs="+",
    type=str,
    default="",
    help="Input label(s) for the propagation summary",
)

p.add_argument(
    "-c",
    "--color",
    nargs="+",
    type=str,
    default="",
    help="Input label(s) for the propagation summary",
)

p.add_argument(
    "-m",
    "--marker",
    nargs="+",
    type=str,
    default=["o"],
    choices=[
        "o",
        "s",
        "D",
        "v",
        "^",
        "<",
        ">",
        "p",
        "P",
        "*",
        "h",
        "H",
        "+",
        "x",
        "X",
        "d",
    ],
    help="Input label(s) for the propagation summary",
)

p.add_argument(
    "-o",
    "--output",
    type=str,
    default="",
    help="Output file base name",
)

args = p.parse_args()

try:

    assert len(args.input) == len(args.label) == len(args.marker)
    
    fig, ax = plt.subplots(1, 1, figsize=(11, 10))

    eta_bins = np.linspace(-4, 4, 100)
    eta_centers = 0.5 * (eta_bins[:-1] + eta_bins[1:])
    eta_width = 8.0 / 100

    for irfile, label, marker in zip(args.input, args.label, args.marker):
        # load the tree
        tree = uproot.open(args.input[ir] + ":" + "propagation_summary")
        # get the numpy arrays
        eta = tree["eta"].array(library="np")
        phi = tree["phi"].array(library="np")
        sens = tree["nSensitives"].array(library="np")
        portals = tree["nPortals"].array(library="np")
        materials = tree["nMaterials"].array(library="np")

        df = pd.DataFrame(
            {
                "eta": eta,
                "phi": phi,
                "sens": sens,
                "portals": portals,
                "materials": materials,
            }
        )

        df["eta_bin"] = np.digitize(eta, bins=eta_bins)

        # grouby bin, so we can calculate stuff
        eta_binned = df.groupby("eta_bin")
        sens_result = eta_binned["sens"].agg(["mean", "sem"])
        sens_result["eta"] = eta_centers
        sens_result["xerr"] = eta_width / 2

        sens_result.plot(
            x="eta",
            y="mean",
            xerr="xerr",
            yerr="sem",
            linestyle="none",
            capsize=1,
            marker=args.marker[ir],
            label=args.label[ir],
            ax=ax,
        )
        ax.set_xlabel(r"$\eta$")
        ax.set_ylabel("Avg. Number of sensitive modules / track")

    fig.show()

except:
    print("The number of input files and labels must match")
    exit()
