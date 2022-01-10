# CONFIGURE JOB

jobname="test"

# uncomment to run workflow locally, set full path to your docker directory
#dockerdir="/home/tamas/Programs/docker"

# uncomment to run workflow on hpc, set full path to your own hpc docker directory (optional)
#dockerdir="/node8_R10/stamas/docker"

# CREATE GENOME LIST

singularity exec "${dockerdir}/r-packages_latest.sif" \
  Rscript ./scripts/create_read_pairs.R $jobname

