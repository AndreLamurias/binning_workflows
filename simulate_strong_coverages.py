import sys
from pathlib import Path
import numpy as np

SEED = 1

mu_p = 1
sigma_p = 0.125
k_p = 1
theta_p = 1
alpha = 1

target_bases = 7.5 * (10 ** 9)

# args:
# 1 - genomes TSV file with species to strains mapping
# 2 - dir with genome files in fasta format
# 3 - output abundances files

species_to_strains = {}
# open mapping file
total_strains = 0
with open(sys.argv[1], "r") as f:
    next(f)  # skip header
    for line in f:
        values = line.split("\t")
        if values[0] not in species_to_strains:
            species_to_strains[values[0]] = []
        species_to_strains[values[0]].append(values[3])
        total_strains += 1

# open dir
genome_files = []
genomes_dir = Path(sys.argv[2]).glob("*.fasta")
for g in genomes_dir:
    genome_files.append(g.stem)
for species in species_to_strains:
    for strain in species_to_strains[species]:
        if strain not in genome_files:
            print(species, strain)
print(total_strains, "genomes from ", len(species_to_strains), "species")

# simulate species
rng = np.random.default_rng(SEED)
species_mean = rng.normal(mu_p, sigma_p, len(species_to_strains))
species_std = rng.gamma(k_p, theta_p, len(species_to_strains))
species_to_y = {}

species_to_y = rng.lognormal(species_mean, species_std)
species_to_y = [y / sum(species_to_y) for y in species_to_y]
# simulate strains
output_file = open(sys.argv[3], "w")
output_file.write("filename\tcov\n")
strain_covs = {}
for i, species in enumerate(species_to_strains):
    strain_p = rng.dirichlet([alpha] * len(species_to_strains[species]))
    strain_cov = strain_p * species_to_y[i]
    for istrain, strain in enumerate(species_to_strains[species]):
        strain_covs[strain] = int(strain_cov[istrain] * target_bases)
        output_file.write("\t".join([strain, str(strain_covs[strain])]) + "\n")
print(sum(strain_covs.values()))
