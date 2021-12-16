# data simulation
# simulate reads according to metagenome simulation file

# requires:
# parallel
# badread

module load parallel/20190122-foss-2018a
module load badread/0.2.0-foss-2020b-Python-3.8.6
#GENOMESDIR=genomes/zymobiomics/
#GENOMESDIR=genomes/strong1/
GENOMESDIR=$3
COMPFILE=$1
OUTPUT=$2
mkdir $2
jobs=$4


#analyse genomes
# python compare_genomes.py genomes/zymobiomics/ zymobiomics_simulation2.tsv mash


# for each line in compfile, simulate reads from first column fasta file
# skip header
all_genomes=()
all_lengths=()
all_seeds=()
all_outputs=()
i=0
while read fastadir comp 
do
    echo $GENOMESDIR/${fastadir}.fasta
    #read_length=$(bc -l <<<"${comp}*${LENGTH}")
    read_length=$comp
    echo ${read_length%%.*}
    all_genomes+=($fastadir)
    if [[ read_length = 0 ]]; then
        read_length="$((10 + $RANDOM % 80))x"
    fi
    all_lengths+=(${read_length%%.*})
    all_seeds+=($i)
    i=$((i+1))
    all_outputs+=($OUTPUT/${fastadir}_reads.fastq)
    #all_genomes+="--reference $GENOMESDIR/${fastadir}.fasta --quantity $read_length"
done < <(tail $COMPFILE -n +2) 
#badread --version
#echo "${all_genomes[@]}"
#echo "${all_lengths[@]}"
#echo "${all_seeds[@]}"
#echo "${all_outputs[@]}"
parallel --link -j$jobs badread simulate --length 10000,7000  --identity 98,99.9,5 --reference $GENOMESDIR/{1}.fasta --quantity {2}  '>' {3} ::: ${all_genomes[@]} ::: ${all_lengths[@]} ::: ${all_outputs[@]} 



