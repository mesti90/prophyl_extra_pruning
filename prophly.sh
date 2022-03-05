# CONFIGURE JOB
## UNZIP REFERENCE FILE - MANUAL

jobname="test"
threads=10
gubbins_iterations=100

# uncomment to run workflow locally, set full path to your docker directory
dockerdir="/home/tamas/Programs/docker"

# uncomment to run workflow on hpc, set full path to your own hpc docker directory (optional)
#dockerdir="/node8_R10/stamas/docker"

# pull docker images
[ ! -d "$dockerdir/staphb" ] && mkdir -p "$dockerdir/staphb"
(cd $dockerdir/staphb && singularity pull docker://staphb/snippy)

[ ! -d "$dockerdir/nanozoo" ] && mkdir -p "$dockerdir/nanozoo"
(cd $dockerdir/nanozoo && singularity pull docker://nanozoo/gubbins)

[ ! -d "$dockerdir/evolbioinfo" ] && mkdir -p "$dockerdir/evolbioinfo"
(cd $dockerdir/evolbioinfo && singularity pull docker://evolbioinfo/treetime:v0.7.4)

[ ! -d "$dockerdir/stitam" ] && mkdir -p "$dockerdir/stitam"
(cd $dockerdir/stitam && singularity pull docker://stitam/r-packages)

# CREATE GENOME LIST

singularity exec "${dockerdir}/stitam/r-packages_latest.sif" \
  Rscript ./scripts/create_read_pairs.R $jobname

## SNIPPY

# create output directory
[ ! -d "./jobs/${jobname}/output/snippy" ] && mkdir -p "./jobs/${jobname}/output/snippy"

# run snippy
reffile=$(cd "./jobs/${jobname}/reference_genome" && dir)
while read col1 col2
do
  singularity exec "${dockerdir}/staphb/snippy_latest.sif" snippy \
    --cpus 10 \
    --outdir "./jobs/${jobname}/output/snippy/$col1" \
    --ref "./jobs/${jobname}/reference_genome/${reffile}" \
    --R1 "./jobs/${jobname}/raw_reads/$col1" \
    --R2 "./jobs/${jobname}/raw_reads/$col2" \
    --force
done < "./jobs/${jobname}/output/read_pairs.tsv"

# run snippy for single end reads

reffile=$(cd "./jobs/${jobname}/reference_genome" && dir)
files=$(cd "./jobs/${jobname}/raw_reads"; dir)
for i in $files
do
  singularity exec "${dockerdir}/staphb/snippy_latest.sif" snippy \
    --cpus 10 \
    --outdir "./jobs/${jobname}/output/snippy/$col1" \
    --ref "./jobs/${jobname}/reference_genome/${reffile}" \
    --se "./jobs/${jobname}/raw_reads/$i" \
    --force
done


## MERGE GENOMES

singularity exec "${dockerdir}/stitam/r-packages_latest.sif" \
  Rscript ./scripts/merge_genomes.R $jobname

# GUBBINS

# create output directory
[ ! -d "./jobs/${jobname}/output/gubbins" ] && mkdir -p "./jobs/${jobname}/output/gubbins"

singularity exec "${dockerdir}/nanozoo/gubbins_latest.sif" run_gubbins.py \
  --tree_builder fasttree \
  --threads $threads \
  --iterations $gubbins_iterations \
  -d \
  "./jobs/${jobname}/output/snippy/consensus.subs.fasta"
  
singularity exec "${dockerdir}/stitam/r-packages_latest.sif" \
  Rscript ./scripts/move_gubbins_results.R $jobname
  
# TREEDATER

singularity exec "${dockerdir}/stitam/r-packages_latest.sif" \
  Rscript ./scripts/run_treedater.R $jobname 0.2
  
# TREETIME

singularity exec "${dockerdir}/evolbioinfo/treetime_v0.7.4.sif" treetime \
  --tree "./jobs/${jobname}/output/gubbins/consensus.subs.node_labelled.final_tree.tre" \
  --aln "./jobs/${jobname}/output/gubbins/consensus.subs.filtered_polymorphic_sites.fasta" \
  --dates "./jobs/${jobname}/dates.tsv" \
  --name-column name \
  --outdir "./jobs/${jobname}/output/treetime"
