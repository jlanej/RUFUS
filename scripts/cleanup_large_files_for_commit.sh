#!/usr/bin/env bash
#
# cleanup_large_files_for_commit.sh
#
# Removes large files from the integration test data directory that should
# not be committed to Git. Run this before committing the test data.
#
# Usage:
#   bash scripts/cleanup_large_files_for_commit.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
DATA_DIR="${REPO_DIR}/integration_test/denovo_data"

echo "Cleaning up large files from integration test data..."

# Remove full reference genome (keep only the region-specific one)
if [[ -d "${DATA_DIR}/reference" ]]; then
    echo "  Removing large reference files..."
    rm -f "${DATA_DIR}/reference/GRCh38_full.fa"
    rm -f "${DATA_DIR}/reference/GRCh38_full.fa.gz"
    rm -f "${DATA_DIR}/reference/GRCh38_full.fa.fai"
fi

# Remove VCF directory if it exists (from discovery mode)
if [[ -d "${DATA_DIR}/vcfs" ]]; then
    echo "  Removing VCF directory..."
    rm -rf "${DATA_DIR}/vcfs"
fi

# Remove any temp files
rm -f "${DATA_DIR}"/*.tmp 2>/dev/null || true
rm -f "${DATA_DIR}/denovo_candidates.vcf.gz"* 2>/dev/null || true

echo ""
echo "Cleanup complete. Remaining files:"
find "${DATA_DIR}" -type f -exec ls -lh {} \; 2>/dev/null

echo ""
echo "Total size:"
du -sh "${DATA_DIR}"

echo ""
echo "You can now commit the integration test data:"
echo "  git add integration_test/denovo_data/"
echo "  git commit -m 'Add integration test data for de novo mutation detection'"

