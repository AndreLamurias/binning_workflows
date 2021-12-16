# Collection of scripts to compare binning approaches

This repository contains shell and python scripts to run, evaluate, and compare different metagenomic binners.
The scripts assume that the correct versions of the binners are installed.


## Summary of scripts
Most of the scripts are based on the `module` framework to work with specific version of bioinformatic tools.
If you have a different setup, you have to adapt the scripts to your needs.

### Preprocessing

These scripts are intended to assemble and prepare ONT reads, as well as polished the assembly with short reads.
A workflow to simulate reads is also provided.

### Working with real ONT reads

Use the `ont_pipeline.sh` script.

```bash
./ont_pipeline.sh samplename path_to_np_reads_file path_to_ilm_reads_dir njobs
```

This script will run metafly assembly, filter, polish, run checkM on all contigs and run MetaBAT with abundances calculated on the
short reads samples.
The polishing is done on the flye assembly graph edges, so that the assembly graph does not have to be converted into a 
contig assembly graph.
That can be acomplished with `generate_contig_graph.py`.

The output is saved to `flye_assembly_<samplename>`.

### Simulating data

Use the `simulate.sh` script to simulate ONT reads.
It takes as input a file with the genome file names and how many reads to simulate from that genome, output directory, and path to genome files, and number of jobs, for example:

```bash
./simulate.sh composition.tsv simulation_output/ genomes_dir/ 10
```

This script generates a fastq for each genome using badread.
It also uses parallel to simulate multiple genomes at the same time.
Then you can use `process_reads.sh`, since the polishing steps of `ont_pipeline.sh` are not necessary.

```bash
./process_reads.sh simulation_output/ genomes.fasta 10
```

`genomes.fasta` is used to map the contigs back to the genomes.
It is the concatenation of all genome files in `genomes_dir/`.
This script also renames the reads so that they contain the name of the original genome.

### Comparing binners

We provide a script to run several binners on the same assembly, and then run dRep and DASTool. 
Once again this script is based on the `module` framework.
This script is not meant to be run linearly, as you may not want to install every binner or you may have different installation configurations or results from previous experiments.
It is meant to show examples of how to run every binner and obtain comparable outputs.
We provide examples for MetaBAT2, MaxBin2, VAMB, GraphBin, SemiBin, dRep and DASTool.




