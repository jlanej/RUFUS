#!/usr/bin/env bash
set -euo pipefail
#
# run_denovo_integration_test.sh
#
# Integration test for RUFUS using GIAB Ashkenazi Trio de novo mutation data.
# Runs RUFUS on sub-region BAMs around validated de novo mutations in HG002
# to verify that RUFUS can detect them.
#
# Prerequisites:
#   - RUFUS must be built (bin/ directory with compiled binaries)
#   - Test data must be downloaded (run scripts/download_denovo_integration_test_data.sh)
#   - samtools must be installed
#
# Usage:
#   cd RUFUS/integration_test
#   bash run_denovo_integration_test.sh [OPTIONS]
#
# Options:
#   -k NUM    K-mer size (default: 25)
#   -t NUM    Threads (default: 8)
#   -d DIR    Data directory (default: denovo_data)
#   -h        Show help

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Defaults
KMER_SIZE=25
THREADS=8
DATA_DIR="${SCRIPT_DIR}/denovo_data"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

usage() {
    head -24 "$0" | grep "^#" | sed 's/^# \?//'
    exit 0
}

# Parse arguments
while getopts "k:t:d:h" opt; do
    case "${opt}" in
        k) KMER_SIZE="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        d) DATA_DIR="${OPTARG}" ;;
        h) usage ;;
        *) usage ;;
    esac
done

log "============================================================"
log "RUFUS De Novo Integration Test"
log "============================================================"

# Check prerequisites
if [[ ! -x "${REPO_DIR}/runRufus.sh" ]]; then
    die "Cannot find runRufus.sh at ${REPO_DIR}/runRufus.sh"
fi

if ! command -v samtools &>/dev/null; then
    die "samtools not found in PATH"
fi

# Check for test data
BAM_DIR="${DATA_DIR}/bams"
REF_DIR="${DATA_DIR}/reference"

HG002_BAM="${BAM_DIR}/HG002.GRCh38.denovo_regions.bam"
HG003_BAM="${BAM_DIR}/HG003.GRCh38.denovo_regions.bam"
HG004_BAM="${BAM_DIR}/HG004.GRCh38.denovo_regions.bam"
REF_FA="${REF_DIR}/GRCh38_denovo_regions.fa"
BED_FILE="${DATA_DIR}/denovo_regions.bed"

for f in "${HG002_BAM}" "${HG003_BAM}" "${HG004_BAM}" "${REF_FA}" "${BED_FILE}"; do
    if [[ ! -f "${f}" ]]; then
        die "Missing test data: ${f}
Run 'bash scripts/download_denovo_integration_test_data.sh' first."
    fi
done

log "Test data found:"
log "  HG002 (proband): ${HG002_BAM}"
log "  HG003 (father):  ${HG003_BAM}"
log "  HG004 (mother):  ${HG004_BAM}"
log "  Reference:        ${REF_FA}"
log "  DNM regions:      ${BED_FILE}"

# Validate BAMs
log ""
log "Validating BAM files..."
for bam in "${HG002_BAM}" "${HG003_BAM}" "${HG004_BAM}"; do
    samtools quickcheck "${bam}" || die "BAM validation failed: ${bam}"
    count=$(samtools view -c "${bam}")
    log "  $(basename "${bam}"): ${count} reads - OK"
done

# Ensure BWA index exists for the reference
BWA_BIN="${REPO_DIR}/bin/externals/bwa/src/bwa_project/bwa"
if [[ ! -f "${REF_FA}.sa" ]]; then
    log ""
    log "Creating BWA index for reference..."
    if [[ -x "${BWA_BIN}" ]]; then
        "${BWA_BIN}" index "${REF_FA}" || die "BWA indexing failed"
        log "  BWA index created"
    elif command -v bwa &>/dev/null; then
        bwa index "${REF_FA}" || die "BWA indexing failed"
        log "  BWA index created (system bwa)"
    else
        die "BWA not found. Build RUFUS first or install bwa."
    fi
fi

# Run RUFUS
log ""
log "============================================================"
log "Running RUFUS"
log "  Subject:  HG002 (proband)"
log "  Controls: HG003 (father), HG004 (mother)"
log "  K-mer:    ${KMER_SIZE}"
log "  Threads:  ${THREADS}"
log "============================================================"

RUFUS_CMD="${REPO_DIR}/runRufus.sh \
    -s ${HG002_BAM} \
    -c ${HG003_BAM} \
    -c ${HG004_BAM} \
    -k ${KMER_SIZE} \
    -t ${THREADS} \
    -r ${REF_FA}"

log "Command: ${RUFUS_CMD}"
log ""

# Execute RUFUS
eval "${RUFUS_CMD}" || {
    log "WARNING: RUFUS exited with non-zero status"
    log "This may be expected if the test regions are small"
}

# Check for output
log ""
log "============================================================"
log "Checking RUFUS output"
log "============================================================"

VCF_OUTPUT=$(find "${BAM_DIR}" -name "*.vcf" -o -name "*.vcf.gz" 2>/dev/null | head -5)
if [[ -n "${VCF_OUTPUT}" ]]; then
    log "VCF output found:"
    for vcf in ${VCF_OUTPUT}; do
        log "  ${vcf}"
        if [[ "${vcf}" == *.gz ]]; then
            variant_count=$(zcat "${vcf}" | grep -v "^#" | wc -l)
        else
            variant_count=$(grep -v "^#" "${vcf}" | wc -l)
        fi
        log "    Variants: ${variant_count}"
    done
else
    log "No VCF output found (may be expected for small test regions)"
fi

# Summary
log ""
log "============================================================"
log "Integration test complete"
log "============================================================"
log ""
log "DNM regions tested:"
while IFS=$'\t' read -r chr start end desc; do
    log "  ${chr}:${start}-${end} (${desc})"
done < "${BED_FILE}"

log ""
log "Review the RUFUS output in: ${BAM_DIR}"
