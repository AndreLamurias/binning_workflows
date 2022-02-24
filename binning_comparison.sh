# metabat, vamb, and GraphMB outputs should already exist
set -x
# assembly dir
jobs=$1

graphbin_path=""
semibin_path=""
gtdb_path="../SemiBin_GTDB/GTDB"
TMPDIR="../temp/"
#assembly_name="edges.fasta"
#depth_name="edges_depth.txt"
#bam_suffix="*edges.bam"

assembly_name="polished.fasta"
depth_name="polished_depth.txt"
bam_suffix="*polished.bam"
assembly_graph="assembly_graph.gfa"

# binning
# metabat
# ls metabat_bins/

#jgi_summarize_bam_contig_depths --outputDepth $depth_name $bam_suffix
#metabat2 -i $assembly_name -a edges_depth.txt -o metabat_bins/bin
cd metabat_bins
grep ">" *.fa | sed -e "s/bin\.//" -e "s/\.fa:>/ /" | awk '{ print $2 "\t" $1}' > ../metabat.tsv; 
rename 's/^bin\./metabat.bin./' *.fa
cd ..


#vamb
cd vamb_output/bins/
rename "s/fna$/fa/" *.fna
rename 's/^(\d+)\.fa$/vamb.$1.fa/' *.fa
grep ">" *.fa | sed -e "s/bin\.//" -e "s/\.fa:>/ /" | awk '{ print $2 "\t" $1}' > ../../vamb.tsv; 
cd ../..

# graphbin
set +x
module load Miniconda3/4.9.2-foss-2020b
set -x
conda init bash
eval "$(conda shell.bash hook)"
conda activate $graphbin_path
#cd $1
tr '\t' ',' < metabat.tsv > metabat.csv
../GraphBin/graphbin --assembler flye --graph $assembly_graph --contigs $assembly_name --binned metabat.csv --output graphbin_output/
cd graphbin_output
tr ',' '\t' < graphbin_output.csv > graphbin_output.tsv
cd ..
python ../write_fasta_bins.py --contigs2bin graphbin_output/graphbin_output.tsv --assembly $assembly_name --outputdir graphbin_bins/ --prefix graphbin --filter 500000
cp graphbin_output/graphbin_output.tsv graphbin.tsv

# run maxbin
cut -f1,3 $depth_name | tail -n+2 > maxbin_cov.tsv
mkdir maxbin_results
set +x
module load MaxBin/2.2.7-foss-2018a-Perl-5.26.1
set -x
pwd
run_MaxBin.pl -contig $assembly_name -abund maxbin_cov.tsv -out maxbin_results/maxbin -thread $jobs
#
cd maxbin_results/; mkdir bins; mv *.fasta bins/; cd bins/; rename "s/fasta$/fa/" *.fasta; cd ../..
cd maxbin_results/bins/; grep ">" *.fa | sed -e "s/bin\.//" -e "s/\.fa:>/ /" | awk '{ print $2 "\t" $1}' > ../../maxbin.tsv; cd ../..
# evaluate maxbin


# run semibin
set +x
module load Miniconda3/4.9.2-foss-2020b
set -x
conda init bash
conda activate $semibin_path
SemiBin single_easy_bin -i polished.fasta -b *.bam  -o semibin_output -t $jobs --recluster -r $gtdb_path --environment ocean
cd semibin_output/output_recluster_bins/; grep ">" *.fa | sed -e "s/bin\.//" -e "s/\.fa:>/ /" | awk '{ print $2 "\t" $1}' > ../../semibin.tsv; cd ..
python ../write_fasta_bins.py --contigs2bin semibin.tsv --assembly polished.fasta --outputdir semibin_bins/ --prefix semibin --filter 500000


# copy all bins to same dir
mkdir all_bins/
cp metabat_bins/* all_bins/
cp vamb_output/bins/* all_bins/
cp maxbin_results/bins/* all_bins/
cp graphbin_bins/* all_bins/

set +x
module load CheckM/1.1.2-foss-2018a-Python-3.6.4
checkm lineage_wf -x fa -t $jobs --tmpdir $TMPDIR --pplacer_threads 10 --reduced_tree --tab_table all_bins checkm -f checkm/checkm.tsv

module load DAS_Tool/1.1.3-foss-2020b-R-4.0.3-Python-3.8.6
DAS_Tool -i metabat.tsv,vamb.tsv,graphbin.tsv,maxbin.tsv -l metabat,vamb,graphbin,maxbin -c $assembly_name -o dastool_out --search_engine diamond -t 20 --write_bins 1 --db_directory $TMPDIR
