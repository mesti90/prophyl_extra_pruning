# CONFIGURE JOB
## UNZIP REFERENCE FILE - MANUAL

jobname="test"
threads=10
gubbins_iterations=100

# uncomment to run workflow locally, set full path to your docker directory
#dockerdir="/home/tamas/Programs/docker"

# uncomment to run workflow on hpc, set full path to your own hpc docker directory (optional)
#dockerdir="/node8_R10/stamas/docker"

# CREATE GENOME LIST

singularity exec "${dockerdir}/r-packages_latest.sif" \
  Rscript ./scripts/create_read_pairs.R $jobname

## SNIPPY

# pull docker image
(cd $dockerdir && singularity pull docker://staphb/snippy)

# create output directory
[ ! -d "./output/${jobname}/snippy" ] && mkdir -p "./output/${jobname}/snippy"

#
reffile=$(cd "./input/${jobname}/reference_genome" && dir)
while read col1 col2
do
  singularity exec "${dockerdir}/snippy_latest.sif" snippy \
    --cpus 10 \
    --outdir "./output/${jobname}/snippy/$col1" \
    --ref "./input/${jobname}/reference_genome/${reffile}" \
    --R1 "./input/${jobname}/raw_reads/$col1" \
    --R2 "./input/${jobname}/raw_reads/$col2" \
    --force
done < "./output/${jobname}/read_pairs.tsv"

## MERGE GENOMES

singularity exec "${dockerdir}/r-packages_latest.sif" \
  Rscript ./scripts/merge_genomes.R $jobname

# GUBBINS

# pull docker image
(cd $dockerdir && singularity pull docker://nanozoo/gubbins)

# create output directory
[ ! -d "./output/${jobname}/gubbins" ] && mkdir -p "./output/${jobname}/gubbins"

singularity exec "${dockerdir}/gubbins_latest.sif" run_gubbins.py \
  --tree_builder fasttree \
  --threads $threads \
  --iterations $gubbins_iterations \
  -d \
  "./output/${jobname}/snippy/consensus.subs.fasta"
