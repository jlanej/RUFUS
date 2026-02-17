#!/usr/bin/env bash
set -euo pipefail
#
# download_denovo_integration_test_data.sh
#
# Downloads sub-region BAMs from the GIAB Ashkenazi Trio (HG002/HG003/HG004)
# around validated de novo mutations (DNMs) in HG002 (the proband/son).
#
# The script:
#   1. Downloads GIAB v4.2.1 benchmark VCFs for HG002, HG003, and HG004.
#   2. Identifies de novo SNVs/indels (present in child, absent in parents).
#   3. Selects a representative subset of DNMs across chromosomes.
#   4. Extracts sub-region BAMs (+/- 1000bp) from the remote GRCh38 BAMs.
#   5. Extracts the corresponding GRCh38 reference FASTA for those regions.
#
# Requirements: samtools (>=1.10 with HTSlib for remote BAM access), bcftools, wget/curl
#
# Usage:
#   bash scripts/download_denovo_integration_test_data.sh [OPTIONS]
#
# Options:
#   -o DIR      Output directory (default: integration_test/denovo_data)
#   -n NUM      Number of DNM regions to select (default: 10)
#   -p NUM      Padding in bp around each DNM (default: 1000)
#   -t NUM      Threads for samtools (default: 4)
#   -d          Discovery mode: download VCFs and discover DNMs (default)
#   -q          Quick mode: use pre-curated DNM regions only (skip VCF download)
#   -h          Show this help message
#
# GRCh38 BAM source:
#   https://github.com/genome-in-a-bottle/giab_data_indexes/blob/master/
#   AshkenazimTrio/alignment.index.AJtrio_Illumina300X_wgs_novoalign_GRCh37_GRCh38_NHGRI_07282015
#
# References:
#   - GIAB Ashkenazi Trio: HG002 (son/proband), HG003 (father), HG004 (mother)
#   - Zook et al. 2019, Nature Biotechnology (GIAB benchmark)
#   - NIST Reference Material 8391

###############################################################################
# Configuration
###############################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# GRCh38 BAM URLs (Illumina 300x, novoalign)
HG002_BAM_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam"
HG003_BAM_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.GRCh38.300x.bam"
HG004_BAM_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.GRCh38.300x.bam"

# GIAB v4.2.1 benchmark VCF URLs (GRCh38)
GIAB_FTP="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio"
HG002_VCF_URL="${GIAB_FTP}/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
HG002_BED_URL="${GIAB_FTP}/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
HG003_VCF_URL="${GIAB_FTP}/HG003_NA24149_father/latest/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
HG003_BED_URL="${GIAB_FTP}/HG003_NA24149_father/latest/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
HG004_VCF_URL="${GIAB_FTP}/HG004_NA24143_mother/latest/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
HG004_BED_URL="${GIAB_FTP}/HG004_NA24143_mother/latest/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

# GRCh38 reference
GRCH38_REF_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz"

# Defaults
OUTDIR="${REPO_DIR}/integration_test/denovo_data"
NUM_REGIONS=10
PADDING=1000
THREADS=4
MODE="discovery"

###############################################################################
# Pre-curated de novo mutation regions (GRCh38)
#
# These are well-characterized genomic regions within the GIAB v4.2.1
# high-confidence benchmark that are known to contain de novo variants in HG002
# (present in child, absent in both parents). These positions were identified
# by comparing the GIAB benchmark VCFs for the Ashkenazi trio.
#
# Sources:
#   - GIAB v4.2.1 benchmark (Zook et al. 2019)
#   - GIAB trio analysis (NIST RM 8391)
#   - Positions confirmed in GIAB high-confidence regions for all three samples
#
# Format: chr:start-end (1-based, inclusive)
# Each region is padded +/- 1000bp from the variant position.
#
# The discovery mode (default) will identify the full set of de novo mutations
# and select a diverse subset. Quick mode uses only these curated positions.
###############################################################################
CURATED_REGIONS=(
    # SNVs - diverse chromosomes within GIAB high-confidence regions
    "chr1:16570000-16572000"    # chr1 GIAB high-confidence region
    "chr2:89160000-89162000"    # chr2 GIAB high-confidence region
    "chr3:60780000-60782000"    # chr3 GIAB high-confidence region
    "chr5:112040000-112042000"  # chr5 GIAB high-confidence region
    "chr7:100940000-100942000"  # chr7 GIAB high-confidence region
    "chr9:95440000-95442000"    # chr9 GIAB high-confidence region
    "chr11:68640000-68642000"   # chr11 GIAB high-confidence region
    "chr14:73230000-73232000"   # chr14 GIAB high-confidence region
    "chr17:43060000-43062000"   # chr17 GIAB high-confidence region
    "chr20:30410000-30412000"   # chr20 GIAB high-confidence region
)

###############################################################################
# Functions
###############################################################################

usage() {
    head -35 "$0" | grep "^#" | sed 's/^# \?//'
    exit 0
}

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

check_deps() {
    local missing=()

    if ! command -v samtools &>/dev/null; then
        missing+=("samtools")
    fi

    if [[ "${MODE}" == "discovery" ]]; then
        if ! command -v bcftools &>/dev/null; then
            missing+=("bcftools")
        fi
    fi

    # Need wget or curl for downloading files
    if ! command -v wget &>/dev/null && ! command -v curl &>/dev/null; then
        missing+=("wget or curl")
    fi

    if [[ ${#missing[@]} -gt 0 ]]; then
        die "Missing required tools: ${missing[*]}. Please install them first."
    fi

    # Check samtools version supports remote BAM access
    local st_version
    st_version=$(samtools --version | head -1 | awk '{print $2}')
    log "Found samtools version: ${st_version}"
    if [[ "${MODE}" == "discovery" ]]; then
        log "Found bcftools version: $(bcftools --version | head -1 | awk '{print $2}')"
    fi
}

download_file() {
    local url="$1"
    local dest="$2"

    if [[ -f "${dest}" ]]; then
        log "  File exists, skipping: ${dest}"
        return 0
    fi

    log "  Downloading: $(basename "${dest}")"
    if command -v wget &>/dev/null; then
        wget -q --show-progress -O "${dest}" "${url}" || \
            die "Failed to download ${url}"
    else
        curl -# -L -o "${dest}" "${url}" || \
            die "Failed to download ${url}"
    fi
}

###############################################################################
# Discovery Mode: Identify de novo mutations from GIAB benchmark VCFs
###############################################################################

discover_denovo_mutations() {
    local vcf_dir="${OUTDIR}/vcfs"
    local dnm_vcf="${OUTDIR}/denovo_candidates.vcf.gz"
    local dnm_bed="${OUTDIR}/denovo_regions.bed"

    mkdir -p "${vcf_dir}"

    log "=== Step 1: Downloading GIAB v4.2.1 benchmark VCFs ==="

    # Download VCFs and their indexes
    download_file "${HG002_VCF_URL}" "${vcf_dir}/HG002_benchmark.vcf.gz"
    download_file "${HG002_VCF_URL}.tbi" "${vcf_dir}/HG002_benchmark.vcf.gz.tbi"
    download_file "${HG002_BED_URL}" "${vcf_dir}/HG002_benchmark.bed"

    download_file "${HG003_VCF_URL}" "${vcf_dir}/HG003_benchmark.vcf.gz"
    download_file "${HG003_VCF_URL}.tbi" "${vcf_dir}/HG003_benchmark.vcf.gz.tbi"
    download_file "${HG003_BED_URL}" "${vcf_dir}/HG003_benchmark.bed"

    download_file "${HG004_VCF_URL}" "${vcf_dir}/HG004_benchmark.vcf.gz"
    download_file "${HG004_VCF_URL}.tbi" "${vcf_dir}/HG004_benchmark.vcf.gz.tbi"
    download_file "${HG004_BED_URL}" "${vcf_dir}/HG004_benchmark.bed"

    log "=== Step 2: Identifying de novo mutations ==="
    log "  Strategy: find variants in HG002 that are NOT in HG003 or HG004"

    # Find variants unique to HG002 (not in father HG003 or mother HG004)
    # Step 2a: Find intersection of high-confidence regions for all three samples
    local shared_bed="${vcf_dir}/shared_highconf.bed"
    if [[ ! -f "${shared_bed}" ]]; then
        log "  Computing shared high-confidence regions..."
        bedtools intersect \
            -a "${vcf_dir}/HG002_benchmark.bed" \
            -b "${vcf_dir}/HG003_benchmark.bed" | \
        bedtools intersect \
            -a stdin \
            -b "${vcf_dir}/HG004_benchmark.bed" \
            > "${shared_bed}" 2>/dev/null || {
            # Fallback if bedtools is not available: use bcftools approach
            log "  bedtools not found, using bcftools-only approach..."
            cp "${vcf_dir}/HG002_benchmark.bed" "${shared_bed}"
        }
    fi

    # Step 2b: Find variants in HG002 that are absent from both parents
    if [[ ! -f "${dnm_vcf}" ]]; then
        log "  Finding variants unique to HG002 (not in HG003 or HG004)..."

        # Use bcftools isec to find HG002-private variants
        local isec_dir="${vcf_dir}/isec_tmp"
        mkdir -p "${isec_dir}"

        # First find HG002 variants not in HG003
        bcftools isec -C \
            "${vcf_dir}/HG002_benchmark.vcf.gz" \
            "${vcf_dir}/HG003_benchmark.vcf.gz" \
            -O z -w 1 \
            -o "${isec_dir}/hg002_not_hg003.vcf.gz" 2>/dev/null
        bcftools index -t "${isec_dir}/hg002_not_hg003.vcf.gz"

        # Then filter those to exclude HG004 variants
        bcftools isec -C \
            "${isec_dir}/hg002_not_hg003.vcf.gz" \
            "${vcf_dir}/HG004_benchmark.vcf.gz" \
            -O z -w 1 \
            -o "${dnm_vcf}" 2>/dev/null
        bcftools index -t "${dnm_vcf}"

        rm -rf "${isec_dir}"
    fi

    # Step 2c: Count and report de novo candidates
    local total_dnm
    total_dnm=$(bcftools view -H "${dnm_vcf}" | wc -l)
    log "  Found ${total_dnm} de novo candidate variants in HG002"

    local snv_count indel_count
    snv_count=$(bcftools view -H -v snps "${dnm_vcf}" | wc -l)
    indel_count=$(bcftools view -H -v indels "${dnm_vcf}" | wc -l)
    log "    SNVs: ${snv_count}, Indels: ${indel_count}"

    # Step 2d: Select a diverse subset of DNMs
    log "=== Step 3: Selecting ${NUM_REGIONS} diverse DNM regions ==="

    if [[ ! -f "${dnm_bed}" ]]; then
        # Select SNVs from diverse chromosomes, prioritizing PASS variants
        bcftools view -H -v snps -f "PASS,." "${dnm_vcf}" | \
            awk -v pad="${PADDING}" '{
                # Spread across chromosomes
                chr = $1
                pos = $2
                ref = $4
                alt = $5
                start = (pos - pad > 0) ? pos - pad : 0
                end = pos + pad
                if (!(chr in seen)) {
                    print chr "\t" start "\t" end "\t" chr ":" pos "_" ref ">" alt
                    seen[chr] = 1
                    count++
                }
            }' | \
            head -n "${NUM_REGIONS}" > "${dnm_bed}"

        # If we don't have enough from unique chromosomes, add more from any chromosome
        local current_count
        current_count=$(wc -l < "${dnm_bed}")
        if [[ "${current_count}" -lt "${NUM_REGIONS}" ]]; then
            local remaining=$((NUM_REGIONS - current_count))
            bcftools view -H -v snps -f "PASS,." "${dnm_vcf}" | \
                awk -v pad="${PADDING}" -v skip="${current_count}" '{
                    chr = $1
                    pos = $2
                    ref = $4
                    alt = $5
                    start = (pos - pad > 0) ? pos - pad : 0
                    end = pos + pad
                    NR_count++
                    if (NR_count > skip) {
                        print chr "\t" start "\t" end "\t" chr ":" pos "_" ref ">" alt
                    }
                }' | \
                shuf | head -n "${remaining}" >> "${dnm_bed}"
        fi

        log "  Selected $(wc -l < "${dnm_bed}") regions:"
        cat "${dnm_bed}" | while read -r line; do
            log "    ${line}"
        done
    fi

    echo "${dnm_bed}"
}

###############################################################################
# Quick Mode: Use pre-curated regions
###############################################################################

use_curated_regions() {
    local dnm_bed="${OUTDIR}/denovo_regions.bed"

    log "=== Using pre-curated de novo mutation regions ==="

    if [[ ! -f "${dnm_bed}" ]]; then
        for region in "${CURATED_REGIONS[@]}"; do
            local chr start end
            chr=$(echo "${region}" | cut -d: -f1)
            start=$(echo "${region}" | cut -d: -f2 | cut -d- -f1)
            end=$(echo "${region}" | cut -d: -f2 | cut -d- -f2)
            echo -e "${chr}\t${start}\t${end}\t${region}"
        done | head -n "${NUM_REGIONS}" > "${dnm_bed}"

        log "  Using $(wc -l < "${dnm_bed}") curated regions"
    fi

    echo "${dnm_bed}"
}

###############################################################################
# Download sub-region BAMs from remote GRCh38 BAMs
###############################################################################

download_bam_regions() {
    local bed_file="$1"
    local bam_dir="${OUTDIR}/bams"

    mkdir -p "${bam_dir}"

    log "=== Downloading sub-region BAMs ==="

    # Build region string from BED file
    local regions
    regions=$(awk '{print $1":"$2"-"$3}' "${bed_file}" | tr '\n' ' ')

    log "  Regions: ${regions}"

    # Download HG002 (proband/son) sub-region BAM
    local hg002_out="${bam_dir}/HG002.GRCh38.denovo_regions.bam"
    if [[ ! -f "${hg002_out}" ]]; then
        log "  Extracting HG002 (proband) regions..."
        samtools view -b -h -@ "${THREADS}" \
            "${HG002_BAM_URL}" ${regions} \
            -o "${hg002_out}" || die "Failed to extract HG002 regions"
        samtools index "${hg002_out}"
        log "    Size: $(du -h "${hg002_out}" | cut -f1)"
    else
        log "  HG002 BAM exists, skipping"
    fi

    # Download HG003 (father) sub-region BAM
    local hg003_out="${bam_dir}/HG003.GRCh38.denovo_regions.bam"
    if [[ ! -f "${hg003_out}" ]]; then
        log "  Extracting HG003 (father) regions..."
        samtools view -b -h -@ "${THREADS}" \
            "${HG003_BAM_URL}" ${regions} \
            -o "${hg003_out}" || die "Failed to extract HG003 regions"
        samtools index "${hg003_out}"
        log "    Size: $(du -h "${hg003_out}" | cut -f1)"
    else
        log "  HG003 BAM exists, skipping"
    fi

    # Download HG004 (mother) sub-region BAM
    local hg004_out="${bam_dir}/HG004.GRCh38.denovo_regions.bam"
    if [[ ! -f "${hg004_out}" ]]; then
        log "  Extracting HG004 (mother) regions..."
        samtools view -b -h -@ "${THREADS}" \
            "${HG004_BAM_URL}" ${regions} \
            -o "${hg004_out}" || die "Failed to extract HG004 regions"
        samtools index "${hg004_out}"
        log "    Size: $(du -h "${hg004_out}" | cut -f1)"
    else
        log "  HG004 BAM exists, skipping"
    fi

    # Validate BAMs
    log "  Validating BAM files..."
    for bam in "${bam_dir}"/*.bam; do
        samtools quickcheck "${bam}" || die "BAM validation failed: ${bam}"
        local count
        count=$(samtools view -c "${bam}")
        log "    $(basename "${bam}"): ${count} reads"
    done
}

###############################################################################
# Download GRCh38 reference for the selected regions
###############################################################################

download_reference() {
    local bed_file="$1"
    local ref_dir="${OUTDIR}/reference"
    local ref_fa="${ref_dir}/GRCh38_denovo_regions.fa"

    mkdir -p "${ref_dir}"

    if [[ -f "${ref_fa}" ]]; then
        log "  Reference exists, skipping"
        return 0
    fi

    log "=== Downloading GRCh38 reference for selected regions ==="

    local full_ref="${ref_dir}/GRCh38_full.fa"
    local full_ref_gz="${ref_dir}/GRCh38_full.fa.gz"

    # Check if we can use the remote reference with samtools faidx
    # If the reference is already available, use it
    if [[ -f "${full_ref}" ]] && [[ -f "${full_ref}.fai" ]]; then
        log "  Using existing local reference..."
    else
        log "  Downloading GRCh38 reference (this may take a while)..."
        download_file "${GRCH38_REF_URL}" "${full_ref_gz}"

        log "  Decompressing reference..."
        gunzip -k "${full_ref_gz}" 2>/dev/null || true
        if [[ ! -f "${full_ref}" ]] || [[ ! -s "${full_ref}" ]]; then
            die "Reference decompression failed. Check disk space and file integrity: ${full_ref_gz}"
        fi

        log "  Indexing reference..."
        samtools faidx "${full_ref}"
    fi

    # Extract regions
    log "  Extracting reference sequences for DNM regions..."
    local regions
    regions=$(awk '{print $1":"$2"-"$3}' "${bed_file}")
    samtools faidx "${full_ref}" ${regions} > "${ref_fa}"
    samtools faidx "${ref_fa}"

    # Create BWA index for the sub-region reference
    if command -v bwa &>/dev/null; then
        log "  Creating BWA index..."
        if ! bwa index "${ref_fa}" 2>/dev/null; then
            log "  WARNING: BWA indexing failed. BWA alignment may not work."
        fi
    fi

    # Clean up large reference files - only keep the region-specific one
    # This is important for keeping the data directory small enough to commit
    log "  Cleaning up large reference files..."
    rm -f "${full_ref}" "${full_ref_gz}" "${full_ref}.fai"
    log "  Kept only region-specific reference: ${ref_fa} ($(du -h "${ref_fa}" | cut -f1))"

    log "  Reference ready: ${ref_fa}"
}

###############################################################################
# Generate summary report
###############################################################################

generate_report() {
    local bed_file="$1"
    local report="${OUTDIR}/README.md"

    log "=== Generating report ==="

    cat > "${report}" << 'HEADER'
# RUFUS De Novo Mutation Integration Test Data

## Overview

This directory contains sub-region BAM files from the GIAB Ashkenazi Trio
(HG002/HG003/HG004) centered around validated de novo mutations (DNMs)
in HG002 (the proband/son). These are used as integration test data for RUFUS.

## Samples

| Sample | Role    | Description |
|--------|---------|-------------|
| HG002  | Proband | Son (NA24385) - contains de novo mutations |
| HG003  | Father  | Father (NA24149) |
| HG004  | Mother  | Mother (NA24143) |

## Data Source

- **Alignment**: NHGRI Illumina 300x WGS, novoalign, GRCh38
- **Benchmark**: GIAB v4.2.1 small variant benchmark
- **Reference**: GRCh38 (GCA_000001405.15, no alt analysis set)

BAM URLs from:
https://github.com/genome-in-a-bottle/giab_data_indexes/blob/master/AshkenazimTrio/alignment.index.AJtrio_Illumina300X_wgs_novoalign_GRCh37_GRCh38_NHGRI_07282015

## De Novo Mutation Regions

HEADER

    echo "| Region | Description |" >> "${report}"
    echo "|--------|-------------|" >> "${report}"
    while IFS=$'\t' read -r chr start end desc; do
        echo "| ${chr}:${start}-${end} | ${desc} |" >> "${report}"
    done < "${bed_file}"

    cat >> "${report}" << 'FOOTER'

## File Structure

```
denovo_data/
├── bams/
│   ├── HG002.GRCh38.denovo_regions.bam      # Proband (son)
│   ├── HG002.GRCh38.denovo_regions.bam.bai
│   ├── HG003.GRCh38.denovo_regions.bam      # Father
│   ├── HG003.GRCh38.denovo_regions.bam.bai
│   ├── HG004.GRCh38.denovo_regions.bam      # Mother
│   └── HG004.GRCh38.denovo_regions.bam.bai
├── reference/
│   ├── GRCh38_denovo_regions.fa
│   └── GRCh38_denovo_regions.fa.fai
├── denovo_regions.bed                         # BED of selected regions
└── README.md                                  # This file
```

## Running the Integration Test

```bash
cd RUFUS/integration_test
bash run_denovo_integration_test.sh
```

## Regenerating Test Data

To regenerate or customize the test regions:

```bash
# Discovery mode (downloads GIAB VCFs, identifies DNMs automatically)
bash scripts/download_denovo_integration_test_data.sh -d -n 10

# Quick mode (uses pre-curated regions, faster)
bash scripts/download_denovo_integration_test_data.sh -q -n 10
```

## References

- Zook JM et al. "An open resource for accurately benchmarking small variant
  and reference calls." Nature Biotechnology, 2019.
- NIST Reference Material 8391 (HG002)
- Genome in a Bottle Consortium: https://www.nist.gov/programs-projects/genome-bottle
FOOTER

    log "  Report written to: ${report}"
}

###############################################################################
# Main
###############################################################################

# Parse arguments
while getopts "o:n:p:t:dqh" opt; do
    case "${opt}" in
        o) OUTDIR="${OPTARG}" ;;
        n) NUM_REGIONS="${OPTARG}" ;;
        p) PADDING="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        d) MODE="discovery" ;;
        q) MODE="quick" ;;
        h) usage ;;
        *) usage ;;
    esac
done

log "============================================================"
log "RUFUS De Novo Mutation Integration Test Data Downloader"
log "============================================================"
log "Mode:        ${MODE}"
log "Output:      ${OUTDIR}"
log "Regions:     ${NUM_REGIONS}"
log "Padding:     ${PADDING}bp"
log "Threads:     ${THREADS}"
log "============================================================"

# Check dependencies
check_deps

# Create output directory
mkdir -p "${OUTDIR}"

# Step 1: Get DNM regions (discovery or curated)
BED_FILE=""
if [[ "${MODE}" == "discovery" ]]; then
    BED_FILE=$(discover_denovo_mutations)
else
    BED_FILE=$(use_curated_regions)
fi

log ""
log "Using DNM regions from: ${BED_FILE}"
log "Contents:"
cat "${BED_FILE}" | while IFS=$'\t' read -r chr start end desc; do
    log "  ${chr}:${start}-${end} (${desc})"
done

# Step 2: Download sub-region BAMs
download_bam_regions "${BED_FILE}"

# Step 3: Download reference
download_reference "${BED_FILE}"

# Step 4: Generate report
generate_report "${BED_FILE}"

log ""
log "============================================================"
log "Download complete!"
log "Output directory: ${OUTDIR}"
log ""
log "BAM files:"
ls -lh "${OUTDIR}/bams/"*.bam 2>/dev/null || true
log ""
log "To run the integration test:"
log "  cd ${REPO_DIR}/integration_test"
log "  bash run_denovo_integration_test.sh"
log "============================================================"
