
#module restore processreads
#COMPFILE=$1
OUTPUT=$1
GENOMEFILE=$2
jobs=$3

FLYEV="Flye/2.9-GCC-10.2.0"
SEQTKV="seqtk/1.3-foss-2018a"
MINIMAPV="Minimap2/2.17-foss-2020b"
CHECKMV="CheckM/1.1.2-foss-2018a-Python-3.6.4"

TEMPDIR="./temp/"

#TODO add NanoPlot
#rename reads
module load $SEQTKV
for readsfile in $OUTPUT/*reads.fastq; do
    # get filename
    genomename=$(basename $readsfile .fastq)
    echo $genomename
    seqtk rename $readsfile  $genomename > $OUTPUT/${genomename}_renamed.fastq
done


# step 2: filter and assemble everything together with flye
module load $FLYEV
flye --nano-raw $OUTPUT/*_reads_renamed.fastq -o flye_assembly_$OUTPUT --threads $jobs --meta

# get edge sequences
cd flye_assembly_${OUTPUT}
awk '/^S/{print ">"$2"\n"$3}' assembly_graph.gfa | fold > edges.fasta
cat edges.fasta | awk '{ if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")} print $0 > filename }'
mkdir edges
mv edge_* edges/

# evaluate edges
module load $CHECKMV
checkm taxonomy_wf --tmpdir $TEMPDIR -t $jobs -x fa domain Bacteria edges/ checkm_edges/

module load $MINIMAPV
# run minimap of edge sequences to genomes to get correct genome
minimap2 -t $jobs -ax asm20 ../$GENOMEFILE edges.fasta > ${OUTPUT%/}.sam
#-ax asm20

#run minimap to map reads to contigs
minimap2 -I 64GB -t $jobs -d edges.mmi edges.fasta; # make index
for readsfile in ../$OUTPUT/*renamed.fastq; do
    # get filename
    genomename=$(basename $readsfile .fastq)
    minimap2 -I 64GB -t $jobs -ax map-ont edges.mmi $readsfile > ${genomename}_${OUTPUT%/}.sam
    cut -f3 ${genomename}_${OUTPUT%/}.sam | grep "edge" | sort | uniq -c > ${genomename}_contigs.txt
done

python ../get_contig_labels_from_reads.py ./
mkdir reads_mapping
mv *_contigs.txt reads_mapping/
mv *reads*.sam reads_mapping/
