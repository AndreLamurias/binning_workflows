set -x
checkmdir=$1
newbins=$2
outdir=$3
jobs=$4
tempdir=/srv/MA/users/andrel/temp/
#run from inside assembly dir
rm all_bins/graphemb.*
rm graphemb_bins/*
sed '/^graphemb/d' $checkmdir/checkm.tsv > $checkmdir/checkm_base.tsv
python3 /srv/MA/users/andrel/write_fasta_bins.py --contigs2bin $newbins --assembly polished.fasta --outputdir graphemb_bins/ --prefix graphemb --filter 500000
rm -r checkm_graphemb/
#module load checkm
set +x
module load CheckM/1.1.2-foss-2018a-Python-3.6.4
set -x
checkm lineage_wf -x fa -t $jobs --tmpdir $tempdir --pplacer_threads 10 --reduced_tree --tab_table graphemb_bins/ checkm_graphemb -f checkm_graphemb/checkm.tsv
cp graphemb_bins/* all_bins/

#head -n1 maxbin_results/checkm/checkm.tsv > maxbin_results/checkm/checkm_all.tsv
#tail -n+2 checkm_all_bins/checkm.tsv >> maxbin_results/checkm/checkm_all.tsv
mkdir $outdir
cp $checkmdir/checkm_base.tsv $outdir/checkm.tsv
tail -n+2 checkm_graphemb/checkm.tsv >>  $outdir/checkm.tsv
set +x
module purge
set -x
cd $outdir
cut -f1,12,13,14 checkm.tsv | sed 's/\t/,/g' | sed 's/Bin Id/genome/ ; s/Completeness/completeness/ ; s/Contamination/contamination/ ; s/Strain heterogeneity/strain_heterogeneity/' > checkm_drep.tsv
sed 's/,/.fa,/' checkm_drep.tsv > checkm_drep_fa.tsv
sed -i 's/genome.fa/genome/' checkm_drep_fa.tsv
cd ..
set +x
module load dRep/2.3.2-foss-2018a-Python-3.6.4
set -x  
dRep dereplicate $outdir -g all_bins/*.fa -p $jobs -comp 50 -con 10 -sa 0.95 --genomeInfo $outdir/checkm_drep_fa.tsv
#module purge



#DAS_Tool -i metabat.tsv,vamb.tsv,graphbin.tsv,$newbins -l metabat,vamb,graphbin,graphemb -c polished.fasta -o $outdir/dastool --search_engine diamond -t $jobs --write_bins 1
