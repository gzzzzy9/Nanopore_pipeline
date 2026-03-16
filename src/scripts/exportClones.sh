#!/usr/bin/sh
#SBATCH -n 10
#SBATCH -p CPU2
#SBATCH -o ../../../slurms_out/mixcr_exportClones_20250930_combined.txt

# activate conda env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate Nanopore

# User-specific params
while getopts "b:i:o:" opt; do
  case $opt in
    b) BATCH="$OPTARG" ;;
    i) DATA_DIR="$OPTARG" ;;
    o) RESULT_DIR="$OPTARG" ;;
    *) echo "Usage: $0 -b batch_id -i input_dir -o output_dir" && exit 1 ;;
  esac
done

# Metadata file
TABLE_FILE="$DATA_DIR/sample_table.csv"
if [[ ! -f $TABLE_FILE ]]; then
    echo "Error: sample_table.csv do not exist!"
    exit 1
fi

tail -n +2 $TABLE_FILE | while IFS=, read -r library barcode isotype species; do
    SAMPLE_DIR=$RESULT_DIR/$library/$barcode
    # export tsv file
    mixcr exportClones \
        -cloneId \
        --drop-default-fields -targetSequences -defaultAnchorPoints \
        -uniqueTagCount Molecule -uniqueTagFraction Molecule -readCount \
        -nFeature CDR3 -allAAFeatures \
        -vHit -dHit -jHit -cHit \
        -vHitScore -dHitScore -jHitScore -cHitScore \
        -vBestIdentityPercent -dBestIdentityPercent -jBestIdentityPercent -cBestIdentityPercent \
        -allNMutations FR1Begin FR3End substitutions \
        -allNMutationsCount FR1Begin FR3End substitutions \
        -allAAMutations FR1Begin FR3End substitutions \
        -allAAMutationsCount FR1Begin FR3End substitutions \
        -isProductive FR1 -isProductive FR2 -isProductive FR3 \
        -isProductive CDR1 -isProductive CDR2 -isProductive CDR3 \
        --force-overwrite \
        $SAMPLE_DIR/mixcr_assemble.clna \
        $SAMPLE_DIR/mixcr_assemble.tsv

    # Build VDJ clone
    python src/utils/build_vdj_clones.py \
    --result-dir "$RESULT_DIR" \
    --library "$library" \
    --barcode "$barcode" \
    --isotype "$isotype" \
    --output "$SAMPLE_DIR/vdj_clone_info.csv"
done