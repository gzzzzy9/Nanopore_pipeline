#!/usr/bin/sh
#SBATCH -n 10
#SBATCH -p CPU2
#SBATCH -o pipeline_out.txt

# activate conda env
source ~/miniconda3/etc/profile.d/conda.sh
conda activate Nanopore

# Specify path for local or slurm execution
if [ -n "$SLURM_SUBMIT_DIR" ]; then
    SCRIPT_DIR="$SLURM_SUBMIT_DIR"
else
    SCRIPT_DIR=$(cd "$(dirname "$0")"; pwd)
fi
SRC_DIR="$SCRIPT_DIR/src"

# set global params
assemble_by="VDJRegion" # align param，please set to 'VDJRegion' if mutation analyses are needs，otherwise set to 'CDR3'.
# RESULT_DIR="/home/huazl/data1/GZY/Nanopore_results" # Directory where the output is stored
# DATA_DIR="/home/huazl/data1/GZY/Nanopore_data" # Directory where the raw data is stored
trimmingQualityThreshold=5 # align param
error_rate=0.2 # cutadapt param

# User-specific params
while getopts "b:i:o:" opt; do
  case $opt in
    b) BATCH="$OPTARG" ;;
    i) DATA_DIR="$OPTARG" ;;
    o) RESULT_DIR="$OPTARG" ;;
    *) echo "Usage: $0 -b batch_name -i input_dir -o output_dir" && exit 1 ;;
  esac
done

# Resolve soft link
DATA_DIR=$(readlink -f "$DATA_DIR")
RESULT_DIR=$(readlink -f "$RESULT_DIR")

# Define 3' primer dict for primer filter
declare -A mouse_3_primers human_3_primers
mouse_3_primers[IgM]="ACCAGAACAACACTGAAGTCAT"
mouse_3_primers[IgG]="GTCAARGGYTAYTTCCCTGAG"
mouse_3_primers[TRB]="AAAGGCTACCCTCGTGTGCTT"
mouse_3_primers[IgK]="GTGCTTCTTGAACAACTTCTACC"
human_3_primers[IgA]="GCTTCTTCCCCCAGGAGCCA"
human_3_primers[IgM]="GATACGAGCAGCGTGGCCGT"
human_3_primers[IgG]="CCCTGGGCTGCCTGGTCAAG"
human_3_primers[IgK]="GCTGAATAACTTCTATCCCA"
human_3_primers[IgL]="AGTGACTTCTACCCGGGAGC"
human_3_primers[TRB]="CACCCAAAAGGCCACACTGG"

# Metadata file
TABLE_FILE="$DATA_DIR/sample_table.csv"
if [[ ! -f $TABLE_FILE ]]; then
    echo "Error: sample_table.csv do not exist!"
    exit 1
fi

tail -n +2 $TABLE_FILE | while IFS=, read -r library barcode isotype species; do
    # # remove \r and \n 
    isotype=${isotype//$'\r'/}
    species=${species//$'\r'/}
    sample_name="$library - $barcode - $isotype - $species"
    echo "----------------------------------------"
    echo "($(date)) Processing $sample_name"
    # # set local params
    # if [[ $species == "Mouse" ]]; then
    #     align_library="Mus_B6"
    #     align_s="mmu_b6"
    # fi
    # if [[ $species == "Human" ]]; then
    #     align_library="Human"
    #     align_s="homo_sapiens"
    # fi
    # # merge *.gz into a single pass.fastq
    # SAMPLE_DATA_DIR=$DATA_DIR/$library/$barcode
    # if [ ! -f "$SAMPLE_DATA_DIR/pass.fastq" ]; then
    #     cat $SAMPLE_DATA_DIR/pass/*.gz | gunzip -c > $SAMPLE_DATA_DIR/pass.fastq
    # fi
    SAMPLE_DIR=$RESULT_DIR/$library/$barcode
    # mkdir -p $SAMPLE_DIR
    # # Calculate raw reads
    # echo "$sample_name Raw reads: $(($(wc -l < $SAMPLE_DATA_DIR/pass.fastq) / 4))"
    # # 1. Filter reads based on length
    # if [ ! -f "$SAMPLE_DIR/filtered_length.fastq" ]; then
    #     filtlong --min_length 500 --max_length 1000 $SAMPLE_DATA_DIR/pass.fastq > $SAMPLE_DIR/filtered_length.fastq 2>/dev/null
    # fi 
    # # Calculate filtered reads
    # echo "$sample_name Length filtered reads: $(($(wc -l < $SAMPLE_DIR/filtered_length.fastq) / 4))"

    # # 2. Align reads with MiXCR and filter reads lacking V gene or CDR3
    # echo "Alignment using params: library: ${align_library}, species: ${align_s}, assembled by ${assemble_by}, quality threshold: $trimmingQualityThreshold"
    # mixcr align --preset generic-amplicon \
    #     --library $align_library \
    #     -s $align_s \
    #     --rna \
    #     -MtrimmingQualityThreshold=$trimmingQualityThreshold \
    #     --floating-left-alignment-boundary \
    #     --floating-right-alignment-boundary C \
    #     --force-overwrite \
    #     --assemble-clonotypes-by $assemble_by \
    #     --not-aligned-R1 $SAMPLE_DIR/na1.fastq \
    #     --drop-non-CDR3-alignments \
    #     -OvParameters.geneFeatureToAlign=VRegionWithP \
    #     --report $SAMPLE_DIR/align_report_1.txt \
    #     $SAMPLE_DIR/filtered_length.fastq \
    #     $SAMPLE_DIR/mixcr_alignment_1.vdjca 
    
    # mixcr exportAlignments --force-overwrite --not-covered-as-empty $SAMPLE_DIR/mixcr_alignment_1.vdjca $SAMPLE_DIR/mixcr_alignment_1.tsv

    # # Execute MiXCR filter 
    # python $SRC_DIR/filter_by_MiXCR.py \
    #     $SAMPLE_DIR/mixcr_alignment_1.tsv \
    #     $SAMPLE_DIR/filtered_length.fastq \
    #     $SAMPLE_DIR/filtered_length_mixcr.fastq
    
    # echo "$sample_name MiXCR filtered reads: $(($(wc -l < $SAMPLE_DIR/filtered_length_mixcr.fastq) / 4))"

    # # 3. Filter reads lacking 5' or 3' primers
    # if [[ $species == "Mouse" ]]; then
    #     echo "$sample_name applys $species 3' primer (${mouse_3_primers[$isotype]}) filter in error rate $error_rate"
    #     cutadapt --report=minimal -e $error_rate \
    #         -g "AGTCGGTAGGATAGCGAGTA...${mouse_3_primers[$isotype]}" \
    #         --rc --discard-untrimmed \
    #         -o $SAMPLE_DIR/filtered_length_mixcr_primer.fastq $SAMPLE_DIR/filtered_length_mixcr.fastq
    # fi
    # if [[ $species == "Human" ]]; then
    #     echo "$sample_name applys $species 3' primer (${human_3_primers[$isotype]}) filter in error rate $error_rate"
    #     cutadapt --report=minimal -e $error_rate \
    #         -g "AGTCGGTAGGATAGCGAGTA...${human_3_primers[$isotype]}" \
    #         --rc --discard-untrimmed \
    #         -o $SAMPLE_DIR/filtered_length_mixcr_primer.fastq $SAMPLE_DIR/filtered_length_mixcr.fastq
    # fi
    # echo "$sample_name Primer filtered reads: $(($(wc -l < $SAMPLE_DIR/filtered_length_mixcr_primer.fastq) / 4))"

    # # 4. Align reads with MiXCR again
    # echo "Alignment with UMI using params: library: ${align_library}, species: ${align_s}, assembled by ${assemble_by}, quality threshold: $trimmingQualityThreshold"
    # mixcr align --preset generic-amplicon-with-umi \
    #     --library $align_library \
    #     -s $align_s \
    #     --rna \
    #     -MtrimmingQualityThreshold=$trimmingQualityThreshold \
    #     --floating-left-alignment-boundary \
    #     --floating-right-alignment-boundary C \
    #     --force-overwrite \
    #     --assemble-clonotypes-by $assemble_by \
    #     --tag-pattern "(UMI:TNNNNTNNNNTNNNNTCTT)(R1:*)" \
    #     --tag-max-budget 10 \
    #     --tag-parse-unstranded \
    #     --not-aligned-R1 $SAMPLE_DIR/na2.fastq \
    #     --drop-non-CDR3-alignments \
    #     -OvParameters.geneFeatureToAlign=VRegionWithP \
    #     -OsaveOriginalReads=true \
    #     --report $SAMPLE_DIR/align_report_2.txt \
    #     $SAMPLE_DIR/filtered_length_mixcr_primer.fastq \
    #     $SAMPLE_DIR/mixcr_alignment_2.vdjca

    # mixcr exportAlignments --force-overwrite --not-covered-as-empty $SAMPLE_DIR/mixcr_alignment_2.vdjca $SAMPLE_DIR/mixcr_alignment_2.tsv    

    # mixcr refineTagsAndSort \
    #     --force-overwrite \
    #     $SAMPLE_DIR/mixcr_alignment_2.vdjca \
    #     $SAMPLE_DIR/mixcr_alignment.corrected.vdjca

    # mixcr exportAlignments --force-overwrite --not-covered-as-empty $SAMPLE_DIR/mixcr_alignment.corrected.vdjca $SAMPLE_DIR/mixcr_alignment.corrected.tsv

    # mixcr assemble -a \
    #     --report $SAMPLE_DIR/assemble_report.txt \
    #     --force-overwrite \
    #     $SAMPLE_DIR/mixcr_alignment.corrected.vdjca \
    #     $SAMPLE_DIR/mixcr_assemble.clna 

    # mixcr exportAlignments --force-overwrite --not-covered-as-empty $SAMPLE_DIR/mixcr_assemble.clna $SAMPLE_DIR/mixcr_assemble_alignment.tsv

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
    # export reads that failed in clone assembly
    mixcr exportReadsForClones -f \
        --id -1 \
        $SAMPLE_DIR/mixcr_assemble.clna \
        $SAMPLE_DIR/na3.fastq

    # Build VDJ clone
    python src/utils/build_vdj_clones.py \
    --result-dir "$RESULT_DIR" \
    --library "$library" \
    --barcode "$barcode" \
    --isotype "$isotype" \
    --output "$SAMPLE_DIR/vdj_clone_info.csv"
done

# echo "Pipeline finished. Generating summary..."
# LOG_FILE=$(scontrol show job $SLURM_JOB_ID | grep -oP 'StdOut=\K\S+')
# python "$SRC_DIR/utils/parsers.py" --batch "$BATCH" --input "$LOG_FILE" --output "$RESULT_DIR/summary.csv"