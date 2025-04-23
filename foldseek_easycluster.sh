#!/bin/bash

# Default values
INPUT_DIR=""
OUTPUT_DIR=""
IDENTITY="0.85"
COVERAGE="0.85"
TMP_DIR="tmp"
FOLDSEEK_PATH="./foldseek/bin/foldseek"

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input_dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        --output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -i|--identity)
            IDENTITY="$2"
            shift 2
            ;;
        -c|--coverage)
            COVERAGE="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check if required arguments are provided
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: $0 --input_dir <input_directory> --output <output_directory> [-i <identity>] [-c <coverage>]"
    echo "Default values: -i $IDENTITY -c $COVERAGE"
    exit 1
fi

# Submit the Slurm job
sbatch << EOF
#!/bin/bash
#SBATCH --time=5:0:0
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=foldseek_cluster

# Set the number of threads for Foldseek
export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK

# Run Foldseek easy-cluster
$FOLDSEEK_PATH easy-cluster \
    $INPUT_DIR $OUTPUT_DIR $TMP_DIR \
    --min-seq-id $IDENTITY \
    -c $COVERAGE \
    -s 10.0 \
    --cov-mode 0 \
    --alignment-type 2 \
    --cluster-mode 0 \
    --cluster-steps 6 \
    --cluster-reassign 1 \
    --lddt-threshold 0.4 \
    --tmscore-threshold 0.4 \
    --mask-bfactor-threshold 0.5 \
    --threads \$SLURM_CPUS_PER_TASK

EOF
