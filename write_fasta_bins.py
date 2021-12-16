import sys
import pdb
import argparse

from pathlib import Path

parser = argparse.ArgumentParser(description="Write bins from TSV file to fasta")
parser.add_argument("--contigs2bin", help="File mapping contig names to bins")
parser.add_argument("--assembly", help="Assembly file with contig sequences")
parser.add_argument("--outputdir", help="Output dir to write bins")
parser.add_argument("--prefix", help="Prefix to add to file names", default="")
parser.add_argument("--filter", help="Exclude bins smaller than this size", default=200000)
args = parser.parse_args()

if args.prefix != "" and not args.prefix.endswith("."):
    args.prefix = args.prefix + "."

print("read assembly file")
contig_seqs = {}
with open(args.assembly, "r") as f:
    for line in f:
        if line.startswith(">"):
            contig_name = line[1:].strip().split(" ")[0]
            contig_seqs[contig_name] = ""
        else:
            contig_seqs[contig_name] += line.strip()
print("read {} sequences".format(len(contig_seqs)))
print("read contig2bin file")
bin2contig = {}
with open(args.contigs2bin, "r") as f:
    for line in f:
        contig_name, bin_name = line.strip().split("\t")
        contig_name = contig_name.split(" ")[0]
        if bin_name not in bin2contig:
            bin2contig[bin_name] = []
        bin2contig[bin_name].append(contig_name)
print("read {} bins".format(len(bin2contig)))

print("write bins fasta")
bin_dir = Path(args.outputdir)
bin_dir.mkdir(parents=True, exist_ok=True)

skipped = 0
[f.unlink() for f in bin_dir.glob("*.fa") if f.is_file()]
for bin in bin2contig:

    if sum([len(contig_seqs.get(c, "")) for c in bin2contig[bin]]) < float(args.filter):
        skipped += 1
        continue
    with open(str(bin_dir / (args.prefix + bin + ".fa")), "w") as f:
        for contig in bin2contig[bin]:
            if len(contig_seqs.get(contig, "")) < 1:
                skipped += 1
                continue
            f.write(">" + contig + "\n")
            f.write(contig_seqs[contig] + "\n")
print("skipped", skipped)
