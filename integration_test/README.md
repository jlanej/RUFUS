# RUFUS De Novo Mutation Integration Test

## Overview

This directory contains scripts and configuration for integration testing RUFUS
using real de novo mutation (DNM) data from the GIAB Ashkenazi Trio:

| Sample | Role    | GIAB ID  | Description |
|--------|---------|----------|-------------|
| HG002  | Proband | NA24385  | Son — contains de novo mutations |
| HG003  | Father  | NA24149  | Father |
| HG004  | Mother  | NA24143  | Mother |

The test uses sub-regions of GRCh38-aligned BAM files (Illumina 300x, novoalign)
centered around validated de novo mutations in HG002.

## Quick Start

### Step 1: Download Test Data

```bash
# From the RUFUS repository root:

# Discovery mode (recommended): automatically identifies de novo mutations
# from GIAB v4.2.1 benchmark VCFs and downloads BAM sub-regions
bash scripts/download_denovo_integration_test_data.sh -d -n 10

# Quick mode (faster): uses pre-curated regions within GIAB high-confidence areas
bash scripts/download_denovo_integration_test_data.sh -q -n 10
```

### Step 2: Run Integration Test

```bash
cd integration_test
bash run_denovo_integration_test.sh
```

## How It Works

### De Novo Mutation Discovery

The download script identifies de novo mutations by comparing the GIAB v4.2.1
benchmark VCFs for all three trio members:

1. **Download** benchmark VCFs for HG002, HG003, and HG004 from GIAB
2. **Intersect** to find variants present in HG002 but absent in both parents
3. **Select** a diverse subset across chromosomes
4. **Extract** BAM sub-regions (±1000bp) from remote GRCh38 BAMs using samtools

### What is a De Novo Mutation?

A de novo mutation is a genetic variant that:
- Is present in the child (HG002) — typically heterozygous
- Is absent in both parents (HG003, HG004) — homozygous reference
- Arose spontaneously during gametogenesis or early embryonic development

The human genome has approximately 50-80 de novo SNVs per generation.

## Data Source

### BAM Files (GRCh38, Illumina 300x)

From the GIAB data index:
```
https://github.com/genome-in-a-bottle/giab_data_indexes/blob/master/
AshkenazimTrio/alignment.index.AJtrio_Illumina300X_wgs_novoalign_GRCh37_GRCh38_NHGRI_07282015
```

| Sample | URL |
|--------|-----|
| HG002  | `ftp://ftp-trace.ncbi.nlm.nih.gov/.../HG002.GRCh38.300x.bam` |
| HG003  | `ftp://ftp-trace.ncbi.nlm.nih.gov/.../HG003.GRCh38.300x.bam` |
| HG004  | `ftp://ftp-trace.ncbi.nlm.nih.gov/.../HG004.GRCh38.300x.bam` |

### Benchmark VCFs (GIAB v4.2.1)

```
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/
├── HG002_NA24385_son/latest/GRCh38/
│   ├── HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
│   └── HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
├── HG003_NA24149_father/latest/GRCh38/
│   └── ...
└── HG004_NA24143_mother/latest/GRCh38/
    └── ...
```

### Reference Genome

```
GRCh38 (GCA_000001405.15, no alt analysis set)
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/
```

## Download Script Options

```
bash scripts/download_denovo_integration_test_data.sh [OPTIONS]

Options:
  -o DIR    Output directory (default: integration_test/denovo_data)
  -n NUM    Number of DNM regions to select (default: 10)
  -p NUM    Padding in bp around each DNM (default: 1000)
  -t NUM    Threads for samtools (default: 4)
  -d        Discovery mode (download VCFs, identify DNMs automatically)
  -q        Quick mode (use pre-curated regions only, skip VCF download)
  -h        Show help
```

## Requirements

- **samtools** ≥ 1.10 (with HTSlib for remote BAM access via HTTPS)
- **bcftools** (for discovery mode VCF processing)
- **wget** or **curl** (for downloading files)
- **bedtools** (optional, improves region intersection in discovery mode)
- Internet access to NCBI FTP servers

## Directory Structure (after download)

```
integration_test/
├── README.md                              # This file
├── run_denovo_integration_test.sh         # Test runner script
└── denovo_data/
    ├── bams/
    │   ├── HG002.GRCh38.denovo_regions.bam      # Proband sub-region BAM
    │   ├── HG002.GRCh38.denovo_regions.bam.bai
    │   ├── HG003.GRCh38.denovo_regions.bam      # Father sub-region BAM
    │   ├── HG003.GRCh38.denovo_regions.bam.bai
    │   ├── HG004.GRCh38.denovo_regions.bam      # Mother sub-region BAM
    │   └── HG004.GRCh38.denovo_regions.bam.bai
    ├── reference/
    │   ├── GRCh38_denovo_regions.fa
    │   └── GRCh38_denovo_regions.fa.fai
    ├── vcfs/                              # (discovery mode only)
    │   ├── HG002_benchmark.vcf.gz
    │   ├── HG003_benchmark.vcf.gz
    │   └── HG004_benchmark.vcf.gz
    ├── denovo_candidates.vcf.gz           # (discovery mode only)
    ├── denovo_regions.bed                 # Selected DNM regions
    └── README.md                          # Auto-generated report
```

## References

- Zook JM et al. "An open resource for accurately benchmarking small variant
  and reference calls." Nature Biotechnology, 2019.
- Zook JM et al. "A robust benchmark for detection of germline large deletions
  and insertions." Nature Biotechnology, 2020.
- NIST Reference Material 8391 — Human DNA for Whole-Genome Variant Assessment
- Genome in a Bottle: https://www.nist.gov/programs-projects/genome-bottle
