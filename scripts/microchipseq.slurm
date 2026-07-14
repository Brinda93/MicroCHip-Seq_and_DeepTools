#!/bin/bash -l
#SBATCH --job-name=microChIP_v2
#SBATCH --output=/gpfs/Labs/Dovat/bxp5423/analysis_microchip_v2/logs/ChIP.%A_%a.out
#SBATCH --error=/gpfs/Labs/Dovat/bxp5423/analysis_microchip_v2/logs/ChIP.%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=compute
#SBATCH --account=lab_dovat
#SBATCH --array=1-40

set -euo pipefail

echo "Job started: $(date)"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"

# =========================
# NEW CLEAN PROJECT DIR
# =========================
BASE="/gpfs/Labs/Dovat/Yen/YT20260413"
OUT="/gpfs/Labs/Dovat/bxp5423/analysis_microchip_v2"

RAW="${OUT}/01_raw"
ALIGN="${OUT}/02_align"
DEDUP="${OUT}/03_dedup"
FILTER="${OUT}/04_filter"
BW="${OUT}/05_bigwig"
PEAKS="${OUT}/06_peaks"
QC="${OUT}/07_qc"
LOGS="${OUT}/logs"

mkdir -p $RAW $ALIGN $DEDUP $FILTER $BW $PEAKS $QC $LOGS

THREADS=8

# =========================
# GENOME
# =========================
BOWTIE2_INDEX="/gpfs/Labs/Dovat/bxp5423/genome/Genome/Hg38/hg38.bowtie_index"
CHROMSIZES="/gpfs/Labs/Dovat/bxp5423/genome/Genome/Hg38/hg38.chrom.sizes"
BLACKLIST="/gpfs/Labs/Dovat/bxp5423/genome/Genome/Hg38/hg38-blacklist.v2.bed"

# =========================
# MODULES
# =========================
module load bowtie2/2.4.2
module load samtools/1.19.2
module load macs2/2.2.9.1
module load bedtools
module load deeptools

export PATH=/ri/shared/modules7/UCSCTools/406/bin/:$PATH

# =========================
# SAMPLE DETECTION (ROBUST)
# =========================
cd "$BASE"

SAMPLES=($(find . -maxdepth 1 -name "*_S*_R1_001.fastq.gz" \
    | sed 's|.*/||' \
    | cut -d'_' -f1 \
    | sort -u))

echo "Total samples detected: ${#SAMPLES[@]}"

if (( SLURM_ARRAY_TASK_ID > ${#SAMPLES[@]} )); then
    echo "No sample for this array ID → exit"
    exit 0
fi

SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}

echo "Processing sample: $SAMPLE"

# =========================
# SAFE FASTQ HANDLING
# =========================
R1=$(ls ${SAMPLE}_S*_R1_001.fastq.gz 2>/dev/null | head -n 1)
R2=$(ls ${SAMPLE}_S*_R2_001.fastq.gz 2>/dev/null | head -n 1)

if [[ -z "$R1" || -z "$R2" ]]; then
    echo "ERROR: Missing FASTQ pair for $SAMPLE → skipping"
    exit 0
fi

echo "R1: $R1"
echo "R2: $R2"

# =========================
# TRIMMING (optional)
# =========================
TRIM1="${QC}/${SAMPLE}_R1.fastq.gz"
TRIM2="${QC}/${SAMPLE}_R2.fastq.gz"

if command -v fastp >/dev/null 2>&1; then
    fastp -i $R1 -I $R2 \
        -o $TRIM1 -O $TRIM2 \
        --qualified_quality_phred 15 \
        --length_required 25 \
        --thread $THREADS \
        2> ${LOGS}/${SAMPLE}_fastp.log
else
    TRIM1=$R1
    TRIM2=$R2
fi

# =========================
# ALIGNMENT
# =========================
SAM="${ALIGN}/${SAMPLE}.sam"
BAM="${ALIGN}/${SAMPLE}.bam"

bowtie2 -x $BOWTIE2_INDEX \
    --very-sensitive \
    -p $THREADS \
    -1 $TRIM1 -2 $TRIM2 \
    -S $SAM \
    2> ${LOGS}/${SAMPLE}_bowtie2.log

samtools view -@ $THREADS -bS $SAM > $BAM
rm -f $SAM

# =========================
# SORT + FIXMATE
# =========================
NS="${DEDUP}/${SAMPLE}.nsort.bam"
FS="${DEDUP}/${SAMPLE}.fix.bam"
SS="${DEDUP}/${SAMPLE}.sorted.bam"

samtools sort -n -@ $THREADS -o $NS $BAM
samtools fixmate -m -@ $THREADS $NS $FS
samtools sort -@ $THREADS -o $SS $FS

rm -f $NS $FS $BAM

# =========================
# DEDUP (micro-ChIP SAFE)
# =========================
DEDUP_BAM="${DEDUP}/${SAMPLE}.dedup.bam"

samtools markdup -s -@ $THREADS $SS $DEDUP_BAM
samtools index $DEDUP_BAM
rm -f $SS

# =========================
# FILTER
# =========================
FILTER_BAM="${FILTER}/${SAMPLE}.filtered.bam"

samtools view -@ $THREADS \
    -b -f 2 -F 1804 -q 20 \
    $DEDUP_BAM > $FILTER_BAM

samtools index $FILTER_BAM
FINAL=$FILTER_BAM

# =========================
# BLACKLIST REMOVE
# =========================
if [[ -f $BLACKLIST ]]; then
    FINAL2="${FILTER}/${SAMPLE}.final.bam"

    bedtools intersect -v \
        -abam $FINAL \
        -b $BLACKLIST > $FINAL2

    samtools index $FINAL2
    FINAL=$FINAL2
fi

# =========================
# QC
# =========================
samtools flagstat $FINAL > ${QC}/${SAMPLE}_flagstat.txt

# =========================
# BIGWIG
# =========================
bamCoverage -b $FINAL \
    -o ${BW}/${SAMPLE}.bw \
    --normalizeUsing CPM \
    --binSize 10 \
    --extendReads \
    -p $THREADS

# =========================
# PEAK CALLING
# =========================
echo "Peak calling: $SAMPLE"

if [[ "$SAMPLE" == *input* || "$SAMPLE" == *Input* ]]; then
    TREAT=$FINAL
    CONTROL=""
else
    BASE_ID=$(echo $SAMPLE | sed 's/_H3.*//;s/_Lamin.*//')
    CONTROL=$(ls ${FILTER}/*${BASE_ID}*input*.bam 2>/dev/null | head -n 1 || true)
    TREAT=$FINAL
fi

# safer MACS2 thresholds for micro-ChIP
if [[ "$SAMPLE" == *H3K27ac* || "$SAMPLE" == *H3K4me3* || "$SAMPLE" == *H3K122ac* ]]; then

    macs2 callpeak \
        -t $TREAT \
        -c ${CONTROL:-} \
        -f BAMPE \
        -g hs \
        -n ${SAMPLE}_narrow \
        --outdir $PEAKS \
        --keep-dup all \
        -q 0.05 \
        2> ${LOGS}/${SAMPLE}_macs2.log

else
    macs2 callpeak \
        -t $TREAT \
        -c ${CONTROL:-} \
        -f BAMPE \
        -g hs \
        -n ${SAMPLE}_broad \
        --outdir $PEAKS \
        --broad \
        --broad-cutoff 0.1 \
        --keep-dup all \
        2> ${LOGS}/${SAMPLE}_macs2.log
fi

echo "DONE $SAMPLE : $(date)"
