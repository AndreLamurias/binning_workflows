import sys
import os
import glob
from pathlib import Path
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
# run mash on every pair of genomes of a dir


def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    #import pdb; pdb.set_trace()
    for s in ax.spines.values():
        s.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


def main():
    genomes_dir = sys.argv[1] # dir where genome files are stored
    genomes_filelist = sys.argv[2] # read genomes list from this file, first column
    action = sys.argv[3] # write dists or draw
    dist_type = sys.argv[4] # mash, ani, or quast
    if dist_type not in ("mash", "ani", "quast_fraction", "quast_indels", "quast_mismatches"):
        print("dist type not found")
        sys.exit()
    if genomes_dir[-1] != "/":
        genomes_dir += "/"
    genomes_filelist = sys.argv[2] # read genomes list from this file, first column
    with open(genomes_filelist, 'r') as f:
        next(f) # skip header
        all_genomes = [genomes_dir + l.split()[0] + ".fasta" for l in f]

    #all_genomes = []
    # get genomes
    #for f in Path(genomes_dir).glob('*.fasta'):
    #    all_genomes.append(str(f))

    print(all_genomes)
    genome_names = {g:g.split("/")[-1].split(".")[0] for g in all_genomes}
    dists_dict = {}

    #    print(output)
    if action == "write":
    #get mash dist for each genome
        # generate mash sketch outside of this script
        if dist_type == "mash":
            for g in all_genomes:
                output = os.popen(f'mash sketch {g}').read()
        for g1 in all_genomes:
            dists_dict[genome_names[g1]] = {}
            for g2 in all_genomes:
                #TODO:
                #stream = os.popen(f'mash dist {g1}.msh {g2}.msh')
                if dist_type == "mash":
                    stream = os.popen(f'mash dist {g1}.msh {g2}.msh')
                    output = stream.readlines()
                    for l in output:
                        if l.startswith(genomes_dir):
                            values = l.split("\t")
                            #genome_pair = genome_names[values[1]]
                            dist = float(values[2])
                            dists_dict[genome_names[g1]][genome_names[g2]] = dist
                elif dist_type == "ani":
                    stream = os.popen(f'rm -f dists.tsv; fastANI -r {g1} -q {g2} -o dists.tsv').read()
                    with open("dists.tsv") as f: # check if file is empty
                        filetext = f.read().strip()
                        if len(filetext) == 0:
                            dists_dict[genome_names[g1]][genome_names[g2]] = 0
                        else:
                            score = float(filetext.split("\t")[-3])
                            dists_dict[genome_names[g1]][genome_names[g2]] = score
                elif "quast" in dist_type:
                    stream = os.popen(f"rm -rf quast_temp/; quast.py {g1} -R {g2} -o quast_temp").read()
                    with open("quast_temp/report.tsv") as f:
                        dists_dict[genome_names[g1]][genome_names[g2]] = -1
                        for l in f:
                            if "Genome fraction" in l and dist_type == "quast_fraction":
                                dists_dict[genome_names[g1]][genome_names[g2]] = float(l.strip().split("\t")[-1])
                            elif "# mismatches per 100 kbp" in l and dist_type == "quast_mismatches":
                                dists_dict[genome_names[g1]][genome_names[g2]] = float(l.strip().split("\t")[-1])
                            elif "# indels per 100 kbp" in l and dist_type == "quast_indels":
                                dists_dict[genome_names[g1]][genome_names[g2]] = float(l.strip().split("\t")[-1])

            if len(dists_dict[genome_names[g1]]) != len(all_genomes):
                import pdb; pdb.set_trace()
        #print(dists_dict)
        genome_names = list(genome_names.values())
        dists_table = np.empty((len(genome_names), len(genome_names)))
        for ig1, g1 in enumerate(genome_names):
            for ig2, g2 in enumerate(genome_names):
                if g2 not in dists_dict[g1]:
                    import pdb; pdb.set_trace()
                dists_table[ig1][ig2] = dists_dict[g1][g2]

        # write to file
        filename = "{}_{}_dists.np".format(genomes_filelist, dist_type)
        with open(filename, 'wb') as f:
            np.save(f, dists_table)
    elif action == "draw":
        # load from file
        filename = "{}_{}_dists.np".format(genomes_filelist, dist_type)
        with open(filename, 'rb') as f:
            dists_table = np.load(f)
        plt.rcParams.update({'font.size': 4})
        fig, ax = plt.subplots(figsize=(12.8,9.6))

        im, cbar = heatmap(dists_table, genome_names, genome_names, ax=ax,
                        cmap="YlGn", cbarlabel=dist_type)
        texts = annotate_heatmap(im, valfmt="{x:.2f}")

        fig.tight_layout()
        plt.savefig("{}_{}_heatmap.png".format(genomes_filelist, dist_type), dpi=200)

if __name__ == "__main__":
    main()
