#!/bin/bash

# Check if at least three arguments were provided
if [ "$#" -lt 3 ]; then
  echo "Usage: $0 <query_dir> <target_dir> <output_file> [min_seq_id] [cluster_threshold]"
  exit 1
fi

# Command line arguments
QUERY_DIR="$1"
TARGET_DIR="$2"
OUTPUT_FILE="$3"

# Optional parameters with default values
MIN_SEQ_ID="${4:-0.8}"
CLUSTER_THRESHOLD="${5:-0.8}"

# Other variables
TMP_DIR="tmp"
FOLDSEEK_PATH="./foldseek/bin/foldseek"

# Submit the Slurm job using a heredoc
sbatch << EOF
#!/bin/bash
#SBATCH --time=10:0:0
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=foldseek_easysearch

# Set the number of threads for Foldseek
export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK

# Run Foldseek multimer search with user-specified directories and parameters
$FOLDSEEK_PATH easy-search \
    "$QUERY_DIR" "$TARGET_DIR" "$OUTPUT_FILE" "$TMP_DIR" \
    --format-output "query,target,pident,fident,theader,tset,taxid,taxname,bits,rmsd,prob" \
    --alignment-type 2 \
    --mask-bfactor-threshold 0.4 \
    --min-seq-id "$MIN_SEQ_ID" \
    --cov-mode 2 \
    -s 10.0 \
    --lddt-threshold 0.2 \
    --tmscore-threshold 0.2 \
    -c "$CLUSTER_THRESHOLD" \
    --threads \$SLURM_CPUS_PER_TASK 

# Sort the results file by the first column and update it in place
sort -k1,1 "$OUTPUT_FILE" -o "$OUTPUT_FILE"

EOF

