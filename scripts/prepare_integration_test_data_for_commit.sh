#!/usr/bin/env bash
set -euo pipefail
#
# prepare_integration_test_data_for_commit.sh
#
# Prepares minimal integration test data that can be committed to GitHub.
# This creates a small subset of test data (3 regions instead of 10) to keep
# file sizes manageable for version control.
#
# The prepared data includes:
#   - Sub-region BAMs (~1-2MB each) for HG002, HG003, HG004
#   - Region-specific reference FASTA (~10KB)
#   - BED file with test regions
#
# Usage:
#   bash scripts/prepare_integration_test_data_for_commit.sh
#
# After running this script, you can commit the integration_test/denovo_data/
# directory to version control.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUTDIR="${REPO_DIR}/integration_test/denovo_data"

# Use only 3 regions to minimize file size for committing
NUM_REGIONS=3
PADDING=500  # Smaller padding to reduce file size
THREADS=4

# GRCh38 BAM URLs
HG002_BAM_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam"
HG003_BAM_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.GRCh38.300x.bam"
HG004_BAM_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.GRCh38.300x.bam"

# GRCh38 reference for extracting regions
GRCH38_REF_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz"

# Minimal curated regions - just 3 diverse chromosomes
MINIMAL_REGIONS=(
    "chr1:16570000-16571000"    # chr1 - small 1kb region
    "chr7:100940000-100941000"  # chr7 - small 1kb region
    "chr20:30410000-30411000"   # chr20 - small 1kb region
)

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

# Check dependencies
check_deps() {
    local missing=()
    command -v samtools &>/dev/null || missing+=("samtools")
    if [[ ${#missing[@]} -gt 0 ]]; then
        die "Missing required tools: ${missing[*]}"
    fi
    log "Found samtools version: $(samtools --version | head -1 | awk '{print $2}')"
}

log "============================================================"
log "Preparing Minimal Integration Test Data for Git Commit"
log "============================================================"
log "Output:      ${OUTDIR}"
log "Regions:     ${NUM_REGIONS} (minimal for small file size)"
log "Padding:     ${PADDING}bp"
log "============================================================"

check_deps

# Clean up any existing data
log "Cleaning up existing data..."
rm -rf "${OUTDIR}"
mkdir -p "${OUTDIR}/bams" "${OUTDIR}/reference"

# Create BED file with minimal regions
BED_FILE="${OUTDIR}/denovo_regions.bed"
log "Creating minimal BED file..."
for region in "${MINIMAL_REGIONS[@]}"; do
    chr=$(echo "${region}" | cut -d: -f1)
    start=$(echo "${region}" | cut -d: -f2 | cut -d- -f1)
    end=$(echo "${region}" | cut -d: -f2 | cut -d- -f2)
    echo -e "${chr}\t${start}\t${end}\t${region}"
done > "${BED_FILE}"

log "  Created ${BED_FILE} with $(wc -l < "${BED_FILE}") regions"

# Build region string
REGIONS=$(awk '{print $1":"$2"-"$3}' "${BED_FILE}" | tr '\n' ' ')
log "  Regions: ${REGIONS}"

# Download sub-region BAMs
log ""
log "=== Downloading minimal sub-region BAMs ==="

for sample in HG002 HG003 HG004; do
    bam_var="${sample}_BAM_URL"
    bam_url="${!bam_var}"
    out_bam="${OUTDIR}/bams/${sample}.GRCh38.denovo_regions.bam"

    log "  Extracting ${sample}..."
    samtools view -b -h -@ "${THREADS}" "${bam_url}" ${REGIONS} -o "${out_bam}"
    samtools index "${out_bam}"

    size=$(du -h "${out_bam}" | cut -f1)
    count=$(samtools view -c "${out_bam}")
    log "    Size: ${size}, Reads: ${count}"
done

# Download reference for regions
log ""
log "=== Downloading reference for selected regions ==="

REF_DIR="${OUTDIR}/reference"
FULL_REF="${REF_DIR}/GRCh38_full.fa"
FULL_REF_GZ="${REF_DIR}/GRCh38_full.fa.gz"
REGION_REF="${REF_DIR}/GRCh38_denovo_regions.fa"

# Check if we can use samtools faidx on remote reference (faster if supported)
# Otherwise download and extract
log "  Downloading GRCh38 reference..."

if command -v wget &>/dev/null; then
    wget -q --show-progress -O "${FULL_REF_GZ}" "${GRCH38_REF_URL}"
else
    curl -# -L -o "${FULL_REF_GZ}" "${GRCH38_REF_URL}"
fi

log "  Decompressing..."
gunzip -k "${FULL_REF_GZ}"

log "  Indexing..."
samtools faidx "${FULL_REF}"

log "  Extracting regions..."
samtools faidx "${FULL_REF}" ${REGIONS} > "${REGION_REF}"
samtools faidx "${REGION_REF}"

# Create BWA index if available
if command -v bwa &>/dev/null; then
    log "  Creating BWA index..."
    bwa index "${REGION_REF}" 2>/dev/null || log "  (BWA indexing skipped)"
fi

# Clean up large files
log "  Cleaning up large reference files..."
rm -f "${FULL_REF}" "${FULL_REF_GZ}" "${FULL_REF}.fai"

log ""
log "  Region reference: ${REGION_REF} ($(du -h "${REGION_REF}" | cut -f1))"

# Generate README
log ""
log "=== Generating README ==="

cat > "${OUTDIR}/README.md" << 'EOF'
# RUFUS De Novo Mutation Integration Test Data

## Overview

This directory contains minimal sub-region BAM files from the GIAB Ashkenazi Trio
(HG002/HG003/HG004) for integration testing. The data is small enough to be
committed to version control.

## Samples

| Sample | Role    | Description |
|--------|---------|-------------|
| HG002  | Proband | Son (NA24385) - contains de novo mutations |
| HG003  | Father  | Father (NA24149) |
| HG004  | Mother  | Mother (NA24143) |

## Data Source

- **Alignment**: NHGRI Illumina 300x WGS, novoalign, GRCh38
- **Reference**: GRCh38 (GCA_000001405.15, no alt analysis set)
- **Regions**: 3 representative 1kb regions for minimal file size

## File Structure

```
denovo_data/
├── bams/
│   ├── HG002.GRCh38.denovo_regions.bam
│   ├── HG002.GRCh38.denovo_regions.bam.bai
│   ├── HG003.GRCh38.denovo_regions.bam
│   ├── HG003.GRCh38.denovo_regions.bam.bai
│   ├── HG004.GRCh38.denovo_regions.bam
│   └── HG004.GRCh38.denovo_regions.bam.bai
├── reference/
│   ├── GRCh38_denovo_regions.fa
│   └── GRCh38_denovo_regions.fa.fai
├── denovo_regions.bed
└── README.md
```

## Running the Integration Test

```bash
cd RUFUS/integration_test
bash run_denovo_integration_test.sh
```

## Regenerating Test Data

To regenerate this minimal test data:

```bash
bash scripts/prepare_integration_test_data_for_commit.sh
```

For more comprehensive testing with more regions:

```bash
bash scripts/download_denovo_integration_test_data.sh -q -n 10
```

## References

- GIAB Ashkenazi Trio: https://www.nist.gov/programs-projects/genome-bottle
- Zook JM et al. Nature Biotechnology, 2019
EOF

log "  README written to: ${OUTDIR}/README.md"

# Summary
log ""
log "============================================================"
log "Preparation complete!"
log "============================================================"
log ""
log "Files ready for commit:"
find "${OUTDIR}" -type f -exec ls -lh {} \; 2>/dev/null | awk '{print "  " $NF ": " $5}'
log ""
log "Total size:"
du -sh "${OUTDIR}"
log ""
log "To commit this data:"
log "  git add integration_test/denovo_data/"
log "  git commit -m 'Add minimal integration test data for de novo mutation detection'"
log "============================================================"

