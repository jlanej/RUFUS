# testRun Expected Output

This document describes the expected output from the `testRun` integration test. Running
`bash runTest.sh` from the `testRun/` directory exercises the full RUFUS pipeline on a
small synthetic trio (Child, Mother, Father) and produces the files listed below.

## Test Command

```bash
cd testRun
bash runTest.sh
```

This invokes:

```
runRufus.sh -s Child.bam -c Mother.bam -c Father.bam -k 25 -t 40 \
  -r ../resources/references/small_test_human_reference_v37_decoys.fa
```

---

## Final Output Files

These are the primary deliverables of a successful run.

| File | Description |
|------|-------------|
| `Child.bam.generator.V2.overlap.hashcount.fastq.bam.FINAL.vcf.gz` | Bgzipped VCF containing filtered, sorted, and deduplicated variant calls |
| `Child.bam.generator.V2.overlap.hashcount.fastq.bam.FINAL.vcf.gz.tbi` | Tabix index for the final VCF |

## Expected Variant Call

The final VCF should contain a **single** de novo variant call. The data line
(excluding header lines that begin with `#`) should look exactly like:

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	./Child.bam	Mother.bam	Father.bam
5:177630000	12896	X-DeNovo	T	G	25	PASS	RN=NODE_Child.bam.generator.V2_0_L273_D22:10:12::MH0;MQ=60;cigar=273M;SB=0.454545;CVT=X;HD=-1_-1_-1_-1_-1_19_-1_19_19_-1_-1_-1_-1_-1_20_20_19_-1_-1_18_-1_18_-1_-1_18_-1_-1_;AO=19;VT=X	GT:DP:RO:AO	0/1:39:20:19	0/0:23:23:0	0/0:23:23:0
```

Key fields:

| Field | Value | Meaning |
|-------|-------|---------|
| CHROM | `5:177630000` | Chromosome 5, region starting at position 177630000 |
| POS | `12896` | Position within the region |
| ID | `X-DeNovo` | Variant type — a de novo mutation |
| REF / ALT | `T` / `G` | Reference allele T, alternate allele G |
| QUAL | `25` | Quality score |
| FILTER | `PASS` | Passed all filters |
| Child genotype | `0/1:39:20:19` | Heterozygous — 39 total depth, 20 ref, 19 alt |
| Mother genotype | `0/0:23:23:0` | Homozygous reference — 23 depth, 0 alt |
| Father genotype | `0/0:23:23:0` | Homozygous reference — 23 depth, 0 alt |

---

## Intermediate Files

The pipeline also produces several intermediate files. These are used during
processing and can be removed with `bash clean.sh`. They are listed here for
troubleshooting purposes.

### Jellyfish K-mer Hash Files

Created during the k-mer counting stage for each sample (Child, Mother, Father):

| File pattern | Description |
|-------------|-------------|
| `<sample>.bam.generator.Jhash` | Jellyfish k-mer hash database |
| `<sample>.bam.generator.Jhash.histo` | K-mer frequency histogram |
| `<sample>.bam.generator.Jhash.histo.7.7.out` | Model fitting output |
| `<sample>.bam.generator.Jhash.histo.7.7.model` | Coverage model (determines mutation threshold) |
| `<sample>.bam.generator.Jhash.histo.7.7.dist` | Distribution file from model |
| `<sample>.bam.generator.Jhash.histo.7.7.prob` | Probability file from model |

For the test data, `<sample>` is `Child`, `Mother`, or `Father`.

### Hash List

| File | Description |
|------|-------------|
| `Child.bam.generator.k25_c<N>.HashList` | K-mers unique to the child (mutant hashes). `<N>` is the minimum coverage cutoff determined by the model. |

### Filtered Reads (RUFUS.Filter output)

| File | Description |
|------|-------------|
| `Child.bam.generator.Mutations.Mate1.fastq` | Forward reads containing mutant k-mers |
| `Child.bam.generator.Mutations.Mate2.fastq` | Reverse reads containing mutant k-mers |

### Overlap Assembly (in `TempOverlap/`)

The overlap assembly stage builds contigs from filtered reads:

| File | Description |
|------|-------------|
| `Child.bam.generator.V2.overlap.fastqd` | Assembled contigs in FASTQD format |
| `Child.bam.generator.V2.overlap.fastq` | Assembled contigs in FASTQ format |
| `Child.bam.generator.V2.overlap.hashcount.fastq` | Contigs annotated with hash counts |
| `Child.bam.generator.V2.overlap.hashcount.fastq.bam` | Contigs aligned to the reference |
| `Child.bam.generator.V2.overlap.hashcount.fastq.bam.bai` | BAM index |

### Raw VCF (before filtering)

| File | Description |
|------|-------------|
| `Child.bam.generator.V2.overlap.hashcount.fastq.bam.vcf` | Raw variant calls from RUFUS.interpret |

### Sorted and Filtered VCF (in `Intermediates/`)

| File | Description |
|------|-------------|
| `Intermediates/Child.bam.generator.V2.overlap.hashcount.fastq.bam.sorted.vcf` | Sorted and filtered VCF (before deduplication) |

### Directories

| Directory | Description |
|-----------|-------------|
| `TempOverlap/` | Temporary files from the overlap assembly stage |
| `Intermediates/` | Intermediate results including reference hashes and per-sample hash comparisons |

---

## Cleaning Up

To remove all intermediate files and keep only the final VCF output:

```bash
cd testRun
bash clean.sh
```

This removes all `*.generator*` files, `mer_counts_merged.jf`, and the
`TempOverlap/` and `Intermediates/` directories.

---

## CI Validation

The CI workflow (`.github/workflows/ci.yml`) validates the test run by:

1. Building the Docker image
2. Running `bash runTest.sh` inside the container
3. Checking that `Child.bam.generator.V2.overlap.hashcount.fastq.bam.vcf` exists
4. Verifying the VCF contains at least one variant call
5. Checking for the expected de novo variant at `5:177630000` position `12896` with ID `X-DeNovo`

If your local test run does not produce these results, something went wrong
during installation. See the main [README](../README.md) for troubleshooting.
