# CONFIGURE JOB
## UNZIP REFERENCE FILE - MANUAL

jobname="test"
threads=30
gubbins_iterations=100
bootstrap_replicates=1000

# uncomment to run workflow locally, set full paths
#dockerdir="/home/tamas/Programs/docker"
#sourcedir="/home/tamas/BRC/aci"

# uncomment to run workflow on hpc, set full paths
dockerdir="/node8_R10/stamas/docker"
sourcedir="/node8_R10/kintses_lab/aci"

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

singularity exec --bind $sourcedir "${dockerdir}/stitam/r-packages_latest.sif" \
  Rscript ./scripts/prep_snippy_input.R $jobname $sourcedir

# SNIPPY

# create output directories
[ ! -d "./jobs/${jobname}/output/snippy" ] && mkdir -p "./jobs/${jobname}/output/snippy"
[ ! -d "./jobs/${jobname}/output/snippy/paired" ] && mkdir -p "./jobs/${jobname}/output/snippy/paired"
[ ! -d "./jobs/${jobname}/output/snippy/single" ] && mkdir -p "./jobs/${jobname}/output/snippy/single"
[ ! -d "./jobs/${jobname}/output/snippy/contig" ] && mkdir -p "./jobs/${jobname}/output/snippy/contig"

reffile=$(cd "./jobs/${jobname}/reference_genome" && dir)

# run snippy for paired reads
while read col1 col2 col3
do
  singularity exec --bind $sourcedir "${dockerdir}/staphb/snippy_latest.sif" snippy \
    --cpus $threads \
    --outdir "./jobs/${jobname}/output/snippy/paired/$col1" \
    --ref "./jobs/${jobname}/reference_genome/${reffile}" \
    --R1 $col2 \
    --R2 $col3 \
    --force
done < "./jobs/${jobname}/output/snippy/paired_reads.tsv"

# run snippy for single end reads
while read col1 col2
do
  singularity exec --bind $sourcedir "${dockerdir}/staphb/snippy_latest.sif" snippy \
    --cpus $threads \
    --outdir "./jobs/${jobname}/output/snippy/single/$col1" \
    --ref "./jobs/${jobname}/reference_genome/${reffile}" \
    --se $col2 \
    --force
done < "./jobs/${jobname}/output/snippy/single_reads.tsv"

# run snippy for contigs
while read col1 col2
do
  singularity exec --bind $sourcedir "${dockerdir}/staphb/snippy_latest.sif" snippy \
    --cpus $threads \
    --outdir "./jobs/${jobname}/output/snippy/contig/$col1" \
    --ref "./jobs/${jobname}/reference_genome/${reffile}" \
    --ctgs $col2 \
    --force
done < "./jobs/${jobname}/output/snippy/contigs.tsv"

singularity exec "${dockerdir}/stitam/r-packages_latest.sif" \
  Rscript ./scripts/merge_snippy_output.R $jobname

# GUBBINS

singularity exec "${dockerdir}/nanozoo/gubbins_latest.sif" run_gubbins.py \
  --tree_builder fasttree \
  --threads $threads \
  --iterations $gubbins_iterations \
  -d \
  "./jobs/${jobname}/output/snippy/consensus.subs.fasta"

singularity exec "${dockerdir}/stitam/r-packages_latest.sif" \
  Rscript ./scripts/move_gubbins_output.R $jobname
  
# IQTREE BOOTSTRAP

singularity exec "${dockerdir}/stitam/r-packages_latest.sif" \
  Rscript ./scripts/prep_iqtree_input.R $jobname

singularity exec "${dockerdir}/staphb/iqtree_latest.sif" iqtree \
  -t "./jobs/${jobname}/output/iqtree/consensus.subs.node_labelled.final_tree.tre" \
  -s "./jobs/${jobname}/output/iqtree/consensus.subs.filtered_polymorphic_sites.fasta" \
  -nt $threads \
  -bb $bootstrap_replicates \
  -wbtl
 
# TREEDATER

singularity exec "${dockerdir}/stitam/r-packages_latest.sif" \
  Rscript ./scripts/run_treedater.R $jobname 0.2
  
## TREETIME

#singularity exec "${dockerdir}/evolbioinfo/treetime_v0.7.4.sif" treetime \
  #--tree "./jobs/${jobname}/output/gubbins/consensus.subs.node_labelled.final_tree.tre" \
  #--aln "./jobs/${jobname}/output/gubbins/consensus.subs.filtered_polymorphic_sites.fasta" \
  #--dates "./jobs/${jobname}/dates.tsv" \
  #--name-column name \
  #--outdir "./jobs/${jobname}/output/treetime"
  
# PASTML

singularity exec "${dockerdir}/evolbioinfo/pastml_latest.sif pastml" \
  -t "./jobs/${jobname}/output/treedater/treedater_tree_with_time.nwk" \
  -d "./jobs/${jobname}/treemeta.tsv" \
  -c "country" \
  -o "./jobs/${jobname}/output/pastml" \
  --work_dir "./jobs/${jobname}/output/pastml" \
  --offline  \
