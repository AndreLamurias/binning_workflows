# run flye on WWTP, polish, calculate abundances, bin for 1 plant

samplename=$1
np_reads_file=$2
ilm_reads_path=$3
jobs=$4

ASSEMBLY=1
FILTERING=1
POLISHING=1
EVAL1=1
BINMETABAT=1

FLYEV="Flye/2.9-GCC-10.2.0"
SEQTKV="seqtk/1.3-foss-2018a"
MINIMAPV="Minimap2/2.17-foss-2020b"
RACONV="Racon/1.3.3-pikachu-foss-2018a"
CHECKMV="CheckM/1.1.2-foss-2018a-Python-3.6.4"
METABATV="MetaBAT/2.12.1-foss-2018a"
SAMTOOLSV="SAMtools/1.11-foss-2020b"
TEMPDIR="./temp/"


if [ $ASSEMBLY -eq 1 ] 
then
    # run flye
    echo "---------------------------------------------------"
    echo "-------------ASSEMBLY------------------------------"
    echo "---------------------------------------------------"


    #FLYE="../Flye/bin/flye" 
    module load $FLYEV
    FLYE="flye"

    $FLYE --nano-raw $np_reads_file -o flye_assembly_$samplename --threads $jobs --meta
    echo flye_assembly_$samplename
fi
cd flye_assembly_$samplename

if [ $FILTERING -eq 1 ]
then
    echo "---------------------------------------------------"
    echo "-------------FILTERING-----------------------------"
    echo "---------------------------------------------------"
    module load $SEQTKV
    awk '/^S/{print ">"$2"\n"$3}' assembly_graph.gfa | fold > edges.fasta
    #filter <1k
    seqtk seq -L 1000 edges.fasta > edges_filt.fasta
    #cat edges_filt.fasta | awk '{ if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")} print $0 > filename }'
    #mkdir edges
    #mv edge_* edges/
fi

module load $MINIMAPV

if [ $POLISHING -eq 1 ]
then
    #polish
    echo "---------------------------------------------------"
    echo "-------------POLISHING-----------------------------"
    echo "---------------------------------------------------"
    
    module load $RACONV
    # map NP reads to contigs
    echo "RACON 1"
    minimap2 -I 64GB -t $jobs -d edges.mmi edges_filt.fasta # make index
    minimap2 -I 64GB -t $jobs -ax map-ont edges.mmi $np_reads_file   > $(basename $np_reads_file .fq)_edges.sam
    racon -u --include-unpolished -t $jobs $np_reads_file $(basename $np_reads_file .fq)_edges.sam edges_filt.fasta > racon_1.fasta

    echo "RACON 2"
    minimap2 -I 64GB -t $jobs -d racon_1.mmi racon_1.fasta # make index
    minimap2 -I 64GB -t $jobs -ax map-ont racon_1.mmi $np_reads_file   > $(basename $np_reads_file .fq)_racon1.sam
    racon -u --include-unpolished -t $jobs $np_reads_file $(basename $np_reads_file .fq)_racon1.sam racon_1.fasta > racon_2.fasta

    echo "RACON 3"
    minimap2 -I 64GB -t $jobs -d racon_2.mmi racon_2.fasta # make index
    minimap2 -I 64GB -t $jobs -ax map-ont racon_2.mmi $np_reads_file   > $(basename $np_reads_file .fq)_racon2.sam
    racon -u --include-unpolished -t $jobs $np_reads_file $(basename $np_reads_file .fq)_racon2.sam racon_2.fasta > racon_3.fasta

    minimap2 -I 64GB -t 30 -d racon_3.mmi racon_3.fasta

    # polish with short reads
    echo "RACON 4 (short reads)"
    cat ${ilm_reads_path}/${samplename}_*ilmtrim.fq > ${samplename}_ilm.fq
    python ../join_reads.py ${samplename}_ilm.fq > ${samplename}_ilm_single.fq
    minimap2 -I 64GB -t 50 -ax sr racon_3.mmi ${samplename}_ilm_single.fq  > ilm_all_racon3.sam
    racon -u --include-unpolished -t 30 ${samplename}_ilm_single.fq ilm_all_racon3.sam racon_3.fasta > racon_4.fasta

    #medaka_consensus -t $jobs -i $np_reads_file -d racon_3.fasta -o medaka -m r941_flip213 -o medaka_consensus1
    #medaka_consensus -t $jobs -i $np_reads_file -d medaka_consensus1/consensus.fasta -o medaka -m r941_flip213 -o medaka_consensus2
    #mv medaka_consensus2/consensus.fasta polished.fasta
fi

if [ $EVAL1 -eq 1 ]
then
    echo "---------------------------------------------------"
    echo "-------------EVAL ASSEMBLY-------------------------"
    echo "---------------------------------------------------"
    module load $CHECKMV
    mv racon_4.fasta polished.fasta
    mkdir edges
    cd edges; cat ../polished.fasta | awk '{ if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")} print $0 > filename }'; cd ..
    find edges/ -name "* *" -type f | rename 's/ /_/g'
    # evaluate edges
    checkm taxonomy_wf --tmpdir $TEMPDIR -t $jobs -x fa domain Bacteria edges/ checkm_edges/
    checkm qa --tmpdir $TEMPDIR -t $jobs checkm_edges/Bacteria.ms checkm_edges/ -f checkm_edges_polished_results.txt --tab_table -o 2
fi

if [[ $BINMETABAT = 1 ]]
then
    echo "---------------------------------------------------"
    echo "-------------METABAT-------------------------------"
    echo "---------------------------------------------------"
    # map NP reads to contigs
    echo "-------------DEPTH----------------------"
    module load $MINIMAPV
    module load $SAMTOOLSV
    minimap2 -I 64GB -t $jobs -d polished.mmi polished.fasta # make index
    minimap2 -I 64GB -t $jobs -ax map-ont polished.mmi $np_reads_file > $(basename $np_reads_file .fq)_polished.sam
    samtools sort $(basename $np_reads_file .fq)_polished.sam > $(basename $np_reads_file .fq)_polished.bam 

    for f in ${ilm_reads_path}/${samplename}*.fq;  do 
        minimap2 -I 64GB -t $jobs -ax sr polished.mmi $f > $(basename $f .fq)_polished.sam
        samtools sort $(basename $f .fq)_polished.sam > $(basename $f .fq)_polished.bam
    done;

    module load $METABATV
    # calculate depth
    jgi_summarize_bam_contig_depths --outputDepth polished_depth.txt *polished.bam
    #jgi_summarize_bam_contig_depths --percentIdentity 97 --outputDepth ilm.tsv il.bam
    #jgi_summarize_bam_contig_depths --percentIdentity 90 --outputDepth np.tsv np.bam
    echo "-------------METABAT---------------------------"
    metabat2 -i polished.fasta -a polished_depth.txt -o metabat_bins_polished/bin

    echo "-------------METABAT EVAL----------------------"
    checkm lineage_wf --tmpdir $TEMPDIR-t $jobs -x fa metabat_bins_polished/ checkm_metabat_polished
    checkm qa --tmpdir $TEMPDIR -t $jobs checkm_metabat_polished/lineage.ms checkm_metabat_polished/ -f checkm_metabat_polished_results.txt --tab_table -o 2
    checkm taxonomy_wf --tmpdir $TEMPDIR -t $jobs -x fa domain Bacteria metabat_bins_polished/ checkm_metabat_polished_bacteria/
    checkm qa --tmpdir $TEMPDIR -t $jobs checkm_metabat_polished_bacteria/Bacteria.ms checkm_metabat_polished_bacteria/ -f checkm_metabat_bacteria_polished_results.txt --tab_table -o 2
    #vamb --outdir path/to/outdir --fasta /path/to/catalogue.fna.gz --jgi edges_depth.txt --minfasta 20000
fi
