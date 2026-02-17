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
