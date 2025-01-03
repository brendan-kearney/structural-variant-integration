#!/bin/bash
#SBATCH --job-name=run_sv_tools    # Job name
#SBATCH --mem-per-cpu=32G         # Memory per CPU
#SBATCH --array=0-11              # Array range for 12 samples
#SBATCH --output=logs/%A_%a.out   # Standard output log
#SBATCH --error=logs/%A_%a.err    # Standard error log

# =======================================================================
# Script: run_sv_tools.sh
# Description: Runs four structural variation (SV) tools: BreakDancer,
#              LUMPY, Delly, and Pindel, on multiple samples using SLURM.
#
# Prerequisites:
# - SLURM workload manager
# - Conda environments for LUMPY and base dependencies
# - Sample list file (one sample name per line)
#
# Usage:
# 1. Update paths and variables in the "Paths and Variables" section.
# 2. Submit the job: `sbatch run_sv_tools.sh`
# =======================================================================

# ================================
# Paths and Variables
# ================================
WORKDIR="/data/reddylab/cegs_ccgr/analyses/shah_acgc/pipeline_gatk/SV_analysis"  # Main working directory
INPUTDIR="/data/reddylab/cegs_ccgr/analyses/shah_acgc/pipeline_gatk/finished_bams"
GENOME_FILE="/data/reddylab/btk20/pipeline_gatk/hg38_genome/hg38.fa"
CONDA_ENV_LUMPY="/data/reddylab/btk20/software/conda_env/py-popgen" # Conda environment for LUMPY
BASE_CONDA_ENV="/data/reddylab/btk20/software/envs/my_env"         # Base Conda environment
SAMPLE_LIST="/data/reddylab/cegs_ccgr/analyses/shah_acgc/pipeline_gatk/sample_names.txt"

# ================================
# Derived Variables
# ================================
SAMPLE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLE_LIST")   # Select sample based on array ID
SAMPLEDIR="$WORKDIR/$SAMPLE"
STATS_FILE="$INPUTDIR/bam_stats/${SAMPLE}.bam.stats"
BAM_FILE="$INPUTDIR/${SAMPLE}.hg38.bam"
OUTPUT_FILE="$SAMPLEDIR/${SAMPLE}.breakdancer.out"
CONFIG_FILE="$SAMPLEDIR/${SAMPLE}.config"
PINDEL_CONFIG="$SAMPLEDIR/pindel_config.txt"

# ================================
# Ensure Directories Exist
# ================================
mkdir -p "$SAMPLEDIR"
mkdir -p logs

# ================================
# Extract Statistics from Stats File
# ================================
avg_ln=$(grep "average length:" "$STATS_FILE" | awk -F: '{print $2}' | awk '{print $1}')
size_mean=$(grep "insert size average:" "$STATS_FILE" | awk -F: '{print $2}' | awk '{print $1}')
stdev=$(grep "insert size standard deviation:" "$STATS_FILE" | awk -F: '{print $2}' | awk '{print $1}')

# ================================
# Define Tool Functions
# ================================
run_breakdancer() {
    echo "Starting BreakDancer for $SAMPLE..."
    "$WORKDIR/bam2cfg.pl" -q 35 "$BAM_FILE" > "$CONFIG_FILE"
    breakdancer-max "$CONFIG_FILE" > "$OUTPUT_FILE"
    awk '!/^#/ && $2 <= $5' "$OUTPUT_FILE" > "$SAMPLEDIR/${SAMPLE}.breakdancer.filtered.out"
    echo "BreakDancer completed for $SAMPLE."
}

run_lumpy() {
    echo "Starting LUMPY for $SAMPLE..."
    source activate "$CONDA_ENV_LUMPY" || echo "Failed to activate Conda environment for LUMPY"
    "$WORKDIR/run_histo.sh" "$BAM_FILE" "$SAMPLEDIR/${SAMPLE}.insert_size_metrics.txt" 5 "$avg_ln"
    lumpy -mw 4 \
        -pe id:$SAMPLE,bam_file:$BAM_FILE,histo_file:$SAMPLEDIR/${SAMPLE}.insert_size_metrics.txt,mean:$size_mean,stdev:$stdev,read_length:$avg_ln,min_non_overlap:101,discordant_z:5,back_distance:10,min_mapping_threshold:20,weight:1 \
        > "$SAMPLEDIR/${SAMPLE}.lumpy.vcf"
    source deactivate
    echo "LUMPY completed for $SAMPLE."
}

run_delly() {
    echo "Starting Delly for $SAMPLE..."
    delly call -g "$GENOME_FILE" "$BAM_FILE" > "$SAMPLEDIR/${SAMPLE}.delly.vcf"
    echo "Delly completed for $SAMPLE."
}

run_pindel() {
    echo "Starting Pindel for $SAMPLE..."
    echo -e "${BAM_FILE}\t${size_mean}\t${SAMPLE}" > "$PINDEL_CONFIG"
    pindel -f "$GENOME_FILE" -i "$PINDEL_CONFIG" -c ALL -o "$SAMPLEDIR/sample_$SAMPLE"
    echo "Pindel completed for $SAMPLE."
}

# ================================
# Run Tools in Parallel
# ================================
source activate "$BASE_CONDA_ENV" || echo "Failed to activate base Conda environment"

run_breakdancer &
PID_BREAKDANCER=$!

run_lumpy &
PID_LUMPY=$!

run_delly &
PID_DELLY=$!

run_pindel &
PID_PINDEL=$!

wait $PID_BREAKDANCER
wait $PID_LUMPY
wait $PID_DELLY
wait $PID_PINDEL

echo "All tools completed for sample: $SAMPLE"
