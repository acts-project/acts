import pandas
import matplotlib.pyplot
import pathlib
import numpy
import scipy

# A list of plots, described as tuples of (category, name, type), where
# category is any of "seeding", "fitting", and "finding"; name is a tree name
# in the output ROOT file, and type is any of "eff", "prod", or "hist".
PLOT_NAMES = [
    ("seeding", "seeding_trackeff_vs_pT", "eff"),
    ("seeding", "seeding_trackeff_vs_eta", "eff"),
    ("seeding", "seeding_trackeff_vs_phi", "eff"),
    ("finding", "finding_trackeff_vs_pT", "eff"),
    ("finding", "finding_trackeff_vs_eta", "eff"),
    ("finding", "finding_trackeff_vs_phi", "eff"),
    ("finding", "finding_nDuplicated_vs_eta", "prof"),
    ("finding", "finding_nFakeTracks_vs_eta", "prof"),
    ("finding", "ndf", "hist"),
    ("finding", "pval", "hist"),
    ("finding", "purity", "hist"),
    ("finding", "completeness", "hist"),
    ("fitting", "res_d0", "hist"),
    ("fitting", "res_z0", "hist"),
    ("fitting", "res_phi", "hist"),
    ("fitting", "res_qop", "hist"),
    ("fitting", "res_qopT", "hist"),
    ("fitting", "res_qopz", "hist"),
    ("fitting", "res_theta", "hist"),
    ("fitting", "pull_d0", "hist"),
    ("fitting", "pull_z0", "hist"),
    ("fitting", "pull_phi", "hist"),
    ("fitting", "pull_qop", "hist"),
    ("fitting", "pull_theta", "hist"),
    ("fitting", "ndf", "hist"),
    ("fitting", "pval", "hist"),
]

# Ratio plots encoded as tuples of two plots, where each individual plot is of
# the form (category, name) as described above.
RATIO_PLOTS = [
    (
        ("seeding", "seeding_trackeff_vs_pT"),
        ("finding", "finding_trackeff_vs_pT"),
    ),
    (
        ("seeding", "seeding_trackeff_vs_eta"),
        ("finding", "finding_trackeff_vs_eta"),
    ),
    (
        ("seeding", "seeding_trackeff_vs_phi"),
        ("finding", "finding_trackeff_vs_phi"),
    ),
]


def make_plots(plot_candidate_names, output_dir, fig_kwargs=None, file_base=None):
    color_spacing = max(len(plot_candidate_names.keys()), 2)

    files_to_copy = []

    if fig_kwargs is None:
        fig_kwargs = {}

    if file_base is None:
        file_base = ""
    else:
        file_base = file_base + "_"

    for data_cat, plot, plot_type in PLOT_NAMES:
        dfs = {
            k: pandas.read_csv(d / data_cat / (plot + ".csv"))
            for (k, (d, _)) in plot_candidate_names.items()
        }

        fig = matplotlib.pyplot.figure(**fig_kwargs)
        ax = fig.subplots()

        for k, v in dfs.items():
            v["center"] = (v["bin_left"] + v["bin_right"]) / 2
            v["width"] = v["center"] - v["bin_left"]
            if "trackeff" in plot:
                dfs[k] = v[v["ntotal"] > 0]

        if plot_type == "eff":
            ax.set_ylabel("Efficiency")
        elif plot == "finding_nDuplicated_vs_eta":
            ax.set_ylabel("Duplicate rate")
        elif plot == "finding_nFakeTracks_vs_eta":
            ax.set_ylabel("Fake rate")
        else:
            ax.set_ylabel("Normalized entries")

        plot_normal = False

        if "vs_eta" in plot:
            ax.set_xlabel("$\\eta$")
        elif "vs_phi" in plot:
            ax.set_xlabel("$\\phi$")
        elif "vs_pT" in plot:
            ax.set_xlabel("$p_T$ (GeV)")
        elif "pval" in plot:
            ax.set_xlabel("$p$")
        elif "ndf" in plot:
            ax.set_xlabel("NDF")
        elif "completeness" in plot:
            ax.set_xlabel("Completeness")
        elif "purity" in plot:
            ax.set_xlabel("Purity")
        elif "res_d0" == plot:
            ax.set_xlabel("Residual $d_0$")
        elif "res_z0" == plot:
            ax.set_xlabel("Residual $z_0$")
        elif "res_phi" == plot:
            ax.set_xlabel("Residual $\\phi$")
        elif "res_qop" == plot:
            ax.set_xlabel("Residual $q/p$")
        elif "res_qopT" == plot:
            ax.set_xlabel("Residual $q/p_T$")
        elif "res_theta" == plot:
            ax.set_xlabel("Residual $\\theta$")
        elif "res_qopz" == plot:
            ax.set_xlabel("Residual $q/p_z$")
        elif "pull_d0" == plot:
            plot_normal = True
            ax.set_xlabel("Pull $d_0$")
        elif "pull_z0" == plot:
            plot_normal = True
            ax.set_xlabel("Pull $z_0$")
        elif "pull_phi" == plot:
            plot_normal = True
            ax.set_xlabel("Pull $\\phi$")
        elif "pull_qop" == plot:
            plot_normal = True
            ax.set_xlabel("Pull $q/p$")
        elif "pull_theta" == plot:
            plot_normal = True
            ax.set_xlabel("Pull $\\theta$")

        if data_cat == "seeding":
            base_color = 0 * color_spacing
        elif data_cat == "finding":
            base_color = 1 * color_spacing
        else:
            base_color = 2 * color_spacing

        kwargs = {}
        plot_kwargs = {k: {} for k in dfs.keys()}
        scale_factors = {k: 1 for k in dfs.keys()}

        if plot_type == "hist" or plot_type == "prof":
            kwargs["drawstyle"] = "steps-mid"
        else:
            kwargs["fmt"] = "."
            for k in dfs.keys():
                plot_kwargs[k]["xerr"] = dfs[k]["width"]

        if plot_type == "eff":
            xkey = "efficiency"
        elif plot_type == "prof":
            xkey = "value"
        elif plot_type == "hist":
            xkey = "ntotal"
            for k in dfs.keys():
                scale_factors[k] = dfs[k][xkey].sum()

        labels = {k: plot_candidate_names[k][1] for k in dfs.keys()}

        if "pull_" in plot or "res_" in plot:
            for k in dfs.keys():
                mean = numpy.average(dfs[k]["center"], weights=dfs[k][xkey])
                std = numpy.sqrt(
                    numpy.average((dfs[k]["center"] - mean) ** 2, weights=dfs[k][xkey])
                )
                labels[k] = labels[k] + "; $\\mu = %.3f$, $\\sigma = %.3f$" % (
                    mean,
                    std,
                )

        if plot_normal:
            x = numpy.linspace(-5, 5, 200)
            ax.plot(x, scipy.stats.norm.pdf(x, 0, 1), label="Ideal", color="black")

        if plot_type == "hist":
            for k in dfs.keys():
                scale_factors[k] *= (dfs[k]["bin_right"] - dfs[k]["bin_left"])[0]

        for i, k in enumerate(dfs.keys()):
            ax.errorbar(
                dfs[k]["center"],
                dfs[k][xkey] / scale_factors[k],
                yerr=(
                    dfs[k]["err_low"] / scale_factors[k],
                    dfs[k]["err_high"] / scale_factors[k],
                ),
                capsize=3,
                label=labels[k],
                color="C%d" % (base_color + i),
                **kwargs,
                **plot_kwargs[k],
            )

        first_key = list(dfs.keys())[0]
        ax.set_xlim(
            xmin=dfs[first_key]["bin_left"].min(),
            xmax=dfs[first_key]["bin_right"].max(),
        )

        ax.legend()
        fig.tight_layout()
        n_data_cat = data_cat
        fig.savefig(
            pathlib.Path(output_dir) / ("%s%s_%s.png" % (file_base, n_data_cat, plot))
        )

        matplotlib.pyplot.close()

        files_to_copy.append(
            "%s/%s%s_%s.png" % (output_dir, file_base, n_data_cat, plot)
        )

    for (from_data_cat, from_plot), (to_data_cat, to_plot) in RATIO_PLOTS:
        dfs = {
            k: (
                pandas.read_csv(d / from_data_cat / (from_plot + ".csv")),
                pandas.read_csv(d / to_data_cat / (to_plot + ".csv")),
            )
            for (k, (d, _)) in plot_candidate_names.items()
        }

        fig = matplotlib.pyplot.figure(**fig_kwargs)
        ax = fig.subplots()

        masks = {}
        ratios = {}
        yerr_lows = {}
        yerr_highs = {}

        for k in dfs.keys():
            dfs[k][0]["center"] = (dfs[k][0]["bin_left"] + dfs[k][0]["bin_right"]) / 2
            dfs[k][0]["width"] = dfs[k][0]["center"] - dfs[k][0]["bin_left"]

            masks[k] = dfs[k][0]["ntotal"] > 0
            ratios[k] = (
                dfs[k][1][masks[k]]["efficiency"] / dfs[k][0][masks[k]]["efficiency"]
            )

            yerr_lows[k] = ratios[k] * numpy.sqrt(
                (dfs[k][0][masks[k]]["err_low"] / dfs[k][0][masks[k]]["efficiency"])
                ** 2
                + (dfs[k][1][masks[k]]["err_low"] / dfs[k][1][masks[k]]["efficiency"])
                ** 2
            )
            yerr_highs[k] = ratios[k] * numpy.sqrt(
                (dfs[k][0][masks[k]]["err_high"] / dfs[k][0][masks[k]]["efficiency"])
                ** 2
                + (dfs[k][1][masks[k]]["err_high"] / dfs[k][1][masks[k]]["efficiency"])
                ** 2
            )

        ax.set_ylabel("Relative efficiency")

        if "vs_eta" in from_plot:
            ax.set_xlabel("$\\eta$")
        elif "vs_phi" in from_plot:
            ax.set_xlabel("$\\phi$")
        elif "vs_pT" in from_plot:
            ax.set_xlabel("$p_T$ (GeV)")

        if from_data_cat == "seeding":
            base_color = 3 * color_spacing
        elif from_data_cat == "finding":
            base_color = 4 * color_spacing
        else:
            base_color = 5 * color_spacing

        for i, k in enumerate(dfs.keys()):
            ax.errorbar(
                dfs[k][0][masks[k]]["center"],
                ratios[k],
                yerr=(yerr_lows[k], yerr_highs[k]),
                capsize=3,
                label=plot_candidate_names[k][1],
                color="C%d" % (base_color + i),
                fmt=".",
                xerr=dfs[k][0][masks[k]]["width"],
            )

        first_key = list(dfs.keys())[0]
        ax.set_xlim(
            xmin=dfs[first_key][0][masks[first_key]]["bin_left"].min(),
            xmax=dfs[first_key][0][masks[first_key]]["bin_right"].max(),
        )

        ax.legend()
        fig.tight_layout()

        vs_name = "%svs%s" % (from_data_cat, to_data_cat)

        fig.savefig(
            pathlib.Path(output_dir) / ("%s%s_%s.png" % (file_base, vs_name, from_plot))
        )
        matplotlib.pyplot.close()

        files_to_copy.append(
            "%s/%s%s_%s.png" % (output_dir, file_base, n_data_cat, plot)
        )

    return files_to_copy
