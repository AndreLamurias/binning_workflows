import argparse
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from tabulate import tabulate

# import pyupset as pyu
from upsetplot import from_contents, plot, UpSet
import venn
import itertools
from sklearn.manifold import TSNE

# from graph_functions import plot_embs
import pickle
import os

SEED = 0

colors = [
    "tab:blue",
    "tab:orange",
    "tab:green",
    "tab:red",
    "tab:purple",
    "tab:brown",
    "tab:pink",
    "tab:gray",
    "tab:olive",
    "tab:cyan",
    "black",
    "silver",
    "maroon",
    "fuchsia",
    "lime",
    "olive",
    "yellow",
    "navy",
    "teal",
    "steelblue",
    "darkred",
    "darkgreen",
    "darkblue",
    "darkorange",
    "lightpink",
    "lightgreen",
    "lightblue",
    "crimson",
    "darkviolet",
    "tomato",
    "tan",
]
markers = [4, 5, 6, 7, 0, 2]


def scatter_hist(x, y, ax, ax_histx, ax_histy, min_comp=50):
    min_comp = 10
    max_cont = 50
    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y)

    # now determine nice limits by hand:
    breakpoint()
    binwidth = 0.25
    xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    lim = (int(xymax / binwidth) + 1) * binwidth

    bins = np.arange(min_comp, lim + binwidth, binwidth)
    ax_histx.hist(x, bins=bins)
    ax_histy.hist(y, bins=bins, orientation="horizontal")


def normalize_method_name(binname):
    if "graph.bin" in binname:
        methodname = "graphbin"
    else:
        methodname = binname.split(".")[0]
    if methodname == "graph" or methodname.startswith("graphemb"):
        methodname = "graphemb"
        binname = "graphemb." + ".".join(binname.split(".")[1:])
    return methodname, binname


def read_drep_cdb(drep_outfile, bins_per_method, methods="all"):
    clusters = {}
    method_to_bins = {}  # which bin sets were identified by each method
    method_to_bins_hq = {}
    missing_bins = []
    with open(drep_outfile, "r") as f:
        next(f)
        for line in f:
            values = line.strip().split(",")
            binname = values[0]
            methodname, binname = normalize_method_name(binname)
            if methods != "all" and methodname not in methods.split(","):
                continue
            cluster = values[-1]
            if cluster not in clusters:
                clusters[cluster] = {}
            clusters[cluster][methodname] = binname
            if methodname not in method_to_bins:
                method_to_bins[methodname] = set()
                method_to_bins_hq[methodname] = set()
            method_to_bins[methodname].add(cluster)
            if binname not in bins_per_method[methodname]:
                print("missing", binname, "checkm")
                missing_bins.append(binname)
                continue
            if bins_per_method[methodname][binname]["comp"] > 90 and bins_per_method[methodname][binname]["cont"] < 5:
                method_to_bins_hq[methodname].add(cluster)
    return clusters, method_to_bins, method_to_bins_hq, missing_bins


def read_checkm(checkm_file, min_comp=50, max_cont=10):

    bins_per_method = {}
    comps = {}
    conts = {}
    with open(checkm_file, "r") as f:
        next(f)
        for line in f:
            binname, comp, cont, strain = line.split(",")
            if "graph.bin" in binname:
                methodname = "graphbin"
            else:
                methodname = binname.split(".")[0]
            if methodname == "graph" or methodname.startswith("graphemb"):
                methodname = "graphemb"
            if methodname not in bins_per_method:
                bins_per_method[methodname] = {}
            if float(comp) >= min_comp and float(cont) <= max_cont:
                bins_per_method[methodname][binname] = {
                    "comp": float(comp),
                    "cont": float(cont),
                    "strain": float(strain),
                }
            comps[binname] = float(comp)
            conts[binname] = float(cont)
    return bins_per_method, comps, conts


def print_stats_table(bins_per_method, out=sys.stdout, max_cont=5):
    rows = [100, 90, 80, 70, 60, 50]
    bar_plot_data = []
    print()
    print("cont<{}".format(max_cont), file=out)
    for j, row in enumerate(rows):
        if j == 0:
            continue
        bar_plot_data.append(["comp>{}".format(rows[j])])
        for i, m in enumerate(bins_per_method):
            # clog_name = log_names[i]
            cont = [bins_per_method[m][b]["cont"] for b in bins_per_method[m]]
            comp = [bins_per_method[m][b]["comp"] for b in bins_per_method[m]]

            filtered_bins = [y for x, y in zip(comp, cont) if x > rows[j] and y < max_cont]
            if len(filtered_bins) > 0:
                avg_cont = round(sum(filtered_bins) / len(filtered_bins), 2)
            else:
                avg_cont = 0
            bar_plot_data[j - 1].append("{} ({})".format(len(filtered_bins), avg_cont))
            # bar_plot_data.append([len([x for x, y in zip(comp, cont) if x >= cutoff and y < 5]) for cutoff in rows])

    print(tabulate(bar_plot_data, headers=list(bins_per_method.keys())), file=out)
    print("# bins (avg contamination)", file=out)


def read_contig2bin(contig2bin_file):
    bins = {}
    with open(contig2bin_file, "r") as f:
        print("opening", contig2bin_file)
        for line in f:
            contig, bin = line.strip().split("\t")
            contig = contig.split()[0]
            base_name = contig2bin_file.split("/")[-1].split(".")[0]
            binname = base_name + "." + str(bin) + ".fa"
            methodname, binname = normalize_method_name(binname)
            if binname not in bins:
                bins[binname] = set()
            bins[binname].add(contig)
    return bins


def compare_clusters(clusters, bins_per_method, missing_bins, out=sys.stdout, target="graphemb"):
    unique_clusters = []
    more_complete_clusters = []
    for c in clusters:
        if target in clusters[c] and clusters[c][target] not in missing_bins:
            if len(clusters[c]) == 1:
                unique_clusters.append(c)
            elif clusters[c] not in missing_bins:
                graphemb_comp = bins_per_method[target][clusters[c][target]]["comp"]
                graphemb_cont = bins_per_method[target][clusters[c][target]]["cont"]
                others_comp = [
                    bins_per_method[m][clusters[c][m]]["comp"]
                    for m in clusters[c]
                    if m != target and clusters[c][m] not in missing_bins
                ]
                others_cont = [
                    bins_per_method[m][clusters[c][m]]["cont"]
                    for m in clusters[c]
                    if m != target and clusters[c][m] not in missing_bins
                ]
                if all(i < graphemb_comp for i in others_comp):  # check only for comp
                    more_complete_clusters.append(c)
                # for method in clusters[c]:
                #    if clusters[c][method] in bins_per_method[method]:
                #        print("    " + method + ": " + str(bins_per_method[method][clusters[c][method]]))
                #    else:
                #        print(clusters[c][method])

    print("clusters unique to {}".format(target), file=out)
    print(
        unique_clusters,
        [(clusters[c], bins_per_method[target][clusters[c][target]]) for c in unique_clusters],
        file=out,
    )
    print("clusters more complete graphemb", file=out)
    print(more_complete_clusters, [clusters[c] for c in more_complete_clusters], file=out)


def get_overlaps(clusters, bins_per_method, method_to_bins, noplots, basedir, comp, cont):

    unique_bins = set()
    # TODO pick labels for venn
    # methods_to_plot = ["graphemb", "maxbin2", "vamb", "metabat2"]
    # methods_to_plot = list(method_to_bins.keys())[:4]
    # breakpoint()
    # labels = venn.get_labels([method_to_bins[m] for m in methods_to_plot], fill=["number", "percent"])
    # fig, ax = venn.venn4(labels, names=methods_to_plot)
    # plt.savefig(basedir + "venn.png")
    # if not noplots:
    #    plt.show()
    # else:
    #    plt.close()
    # labels = venn.get_labels([method_to_bins_hq[m] for m in methods_to_plot], fill=["number", "percent"])
    # fig, ax = venn.venn4(labels, names=methods_to_plot)
    # if not args.noplots:
    #    plt.show()
    # else:
    #    plt.close()
    # write UpSet format
    out = open(basedir + "bins.tsv", "w")
    contents = {m: [] for m in method_to_bins}
    for c in clusters:
        methods_with_hq_bin = set()
        for m in method_to_bins.keys():
            # check if any are HQ bins in this cluster
            if (
                m in clusters[c]
                and clusters[c][m] in bins_per_method[m]  # TODO why would this be false
                and bins_per_method[m][clusters[c][m]]["comp"] > comp
                and bins_per_method[m][clusters[c][m]]["cont"] < cont
            ):
                methods_with_hq_bin.add(clusters[c][m])
                # write to file if method m is in the cluster
                if m in clusters[c]:
                    contents[m].append(c)
        # print and save to file cluster where only one method has HQ bin
        if len(methods_with_hq_bin) == 1:
            best_method = methods_with_hq_bin.pop()
            print(
                "unique method with HQ bin",
                best_method,
                "scores:",
                str([(m, bins_per_method[m][clusters[c][m]]) for m in method_to_bins.keys() if m in clusters[c]]),
                "cluster",
                c,
            )
            print(
                best_method,
                "\t",
                "\t".join(
                    [
                        m + ":" + str(bins_per_method[m][clusters[c][m]])
                        for m in method_to_bins.keys()
                        if m in clusters[c]
                    ]
                ),
                file=out,
            )
            unique_bins.add(best_method)
    upset = from_contents(contents)
    UpSet(upset, subset_size="count", show_counts=True, show_percentages=True, sort_by="cardinality").plot()
    plt.show()
    UpSet(upset, subset_size="count", show_counts=True, show_percentages=True, sort_by="cardinality").plot()
    plt.savefig(basedir + "upset.png")
    plt.close()
    print(",".join([m for m in unique_bins if m.startswith("graphemb")]))
    return unique_bins


def plot_scatter(bins_per_method, methods, noplots):
    hist_bin_size = 10
    # scatter plot each method
    # based on https://stackoverflow.com/a/49469428
    fig, axs = plt.subplots(2, (len(bins_per_method) // 2) + 1)
    for i, m in enumerate(bins_per_method):
        if methods != "all" and m not in methods.split(","):
            continue
        cont = [bins_per_method[m][b]["cont"] for b in bins_per_method[m]]
        comp = [bins_per_method[m][b]["comp"] for b in bins_per_method[m]]
        print("plotting", m, i, colors[i])
        # axs[i % 2, i // 2].scatter(comp, cont, c=colors[i], label=m, alpha=0.5)
        # axs[i % 2, i // 2].set_title(m)
        # definitions for the axes
        gs = gridspec.GridSpec(3, 3)
        ax_main = plt.subplot(gs[1:3, :2])
        # ax_main.set(xlabel="Completeness", ylabel="Contamination")
        ax_main.set_title(m)
        ax_xDist = plt.subplot(gs[0, :2], sharex=ax_main)
        ax_yDist = plt.subplot(gs[1:3, 2], sharey=ax_main)
        ax_main.scatter(comp, cont)  # , marker=".")
        ax_main.set(xlabel="Completeness", ylabel="y Contamination")

        ax_xDist.hist(comp, bins=hist_bin_size, align="mid")
        ax_xDist.set(ylabel="count")
        ax_xDist.label_outer()
        ax_xCumDist = ax_xDist.twinx()
        ax_xCumDist.hist(
            comp, bins=hist_bin_size, cumulative=-1, histtype="step", density=False, color="r", align="mid"
        )
        ax_xCumDist.tick_params("y", colors="r")
        ax_xCumDist.set_ylabel("cumulative", color="r")

        ax_yDist.hist(cont, bins=hist_bin_size, orientation="horizontal", align="mid")
        ax_yDist.set(xlabel="count")
        ax_yDist.label_outer()
        ax_yCumDist = ax_yDist.twiny()
        ax_yCumDist.hist(
            cont,
            bins=hist_bin_size,
            cumulative=True,
            histtype="step",
            density=False,
            color="r",
            align="mid",
            orientation="horizontal",
        )
        ax_yCumDist.tick_params("x", colors="r")
        ax_yCumDist.set_xlabel("cumulative", color="r")
        if not noplots:
            plt.show()
        else:
            plt.close()
        # breakpoint()
        # plt.scatter(comp, cont, c=colors[i], label=m, alpha=0.5)
    plt.figure()
    for i, m in enumerate(bins_per_method):
        if methods != "all" and m not in methods.split(","):
            continue
        cont = [bins_per_method[m][b]["cont"] for b in bins_per_method[m]]
        comp = [bins_per_method[m][b]["comp"] for b in bins_per_method[m]]
        print("plotting", m, i, colors[i])
        plt.scatter(comp, cont, c=colors[i], label=m, alpha=0.5, marker=markers[i], s=100)
        plt.legend()
    # for ax in axs.flat:
    #    ax.set(xlabel="Completeness", ylabel="Contamination")
    # Hide x labels and tick labels for top plots and y ticks for right plots.
    # for ax in axs.flat:
    #    ax.label_outer()
    # plt.legend()
    if not noplots:
        plt.show()
    else:
        plt.close()


# first use checkm_drep_fa.tsv to get bins of each approach
def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--basedir", help="basedir")
    parser.add_argument("--checkm", help="CheckM output")
    parser.add_argument("--drep", help="dRep output (Cdb.csv)")
    parser.add_argument("--noplots", action="store_true", help="Do not show plots")
    parser.add_argument("--methods", help="Methods to include, or all, comma separated", default="all")
    parser.add_argument("--bins", nargs="+")

    args = parser.parse_args()
    checkm_file = os.path.join(args.basedir, args.checkm)
    drep_outfile = os.path.join(args.basedir, args.drep)  # Cdb.csv
    bins_per_method, comps, conts = read_checkm(checkm_file)
    # load embs and plot clusters
    bins = {}
    for binfile in args.bins:
        bins.update(read_contig2bin(binfile))
    for m in bins_per_method:
        for b in bins_per_method[m]:
            if b in bins:
                bins_per_method[m][b]["size"] = len(bins[b])

    # plot_scatter(bins_per_method, args.methods, args.noplots)
    fout = open(args.basedir + "/hq_stats.tsv", "w")
    # print_stats_table(bins_per_method, max_cont=1)
    # print_stats_table(bins_per_method, max_cont=2)
    print_stats_table(bins_per_method, max_cont=5)
    print_stats_table(bins_per_method, max_cont=5, out=fout)
    print_stats_table(bins_per_method, max_cont=10)
    fout.close()
    clusters, method_to_bins, method_to_bins_hq, missing_bins = read_drep_cdb(
        drep_outfile, bins_per_method, args.methods
    )

    fout = open(args.basedir + "/relevant_bins.txt", "w")
    compare_clusters(clusters, bins_per_method, missing_bins, target="graphemb", out=fout)
    fout.close()

    # todo: also get bin composition
    print("all overlaps")
    get_overlaps(clusters, bins_per_method, method_to_bins, args.noplots, args.basedir + "all_", comp=50, cont=10)
    print("HQ overlaps")
    unique_bins = get_overlaps(
        clusters, bins_per_method, method_to_bins_hq, args.noplots, args.basedir + "hq_", comp=90, cont=5
    )
    # for m in method_to_bins:
    #    print(m)
    #    print(",".join(method_to_bins_hq[m]))

    # breakpoint()
    cluster_to_contig = {}
    cluster_to_bins = {c: [bins[clusters[c][m]] for m in clusters[c] if clusters[c][m] in bins] for c in clusters}
    for c in cluster_to_bins:
        cluster_to_contig[c] = []
        for m in cluster_to_bins[c]:
            cluster_to_contig[c] += m


if __name__ == "__main__":
    main()
