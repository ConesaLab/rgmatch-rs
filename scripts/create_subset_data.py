#!/usr/bin/env python3
"""
Create subset test data for CI integration tests.
Takes 100 peaks per chromosome and extracts overlapping GTF annotations.
"""

import os
import sys
from collections import defaultdict

# Configuration
PEAKS_PER_CHROM = 100
BUFFER_SIZE = 50000  # 50kb buffer for GTF extraction
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# Paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(PROJECT_DIR, "tests", "data")

BED_INPUT = os.path.join(DATA_DIR, "full_peaks.bed")
GTF_INPUT = os.path.join(DATA_DIR, "full_genome.gtf")
BED_OUTPUT = os.path.join(DATA_DIR, "subset_peaks.bed")
GTF_OUTPUT = os.path.join(DATA_DIR, "subset_genome.gtf")


def create_subset_bed():
    """Create subset BED file with 100 peaks per chromosome."""
    print(f"Creating subset BED file: {BED_OUTPUT}")

    chrom_counts = defaultdict(int)
    selected_peaks = []

    with open(BED_INPUT, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue

            chrom = fields[0]
            if chrom in CHROMOSOMES and chrom_counts[chrom] < PEAKS_PER_CHROM:
                selected_peaks.append(line)
                chrom_counts[chrom] += 1

    # Sort peaks by chromosome order then by position
    def sort_key(line):
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = int(fields[1])
        # Sort by chromosome order in CHROMOSOMES list
        chrom_idx = CHROMOSOMES.index(chrom) if chrom in CHROMOSOMES else 999
        return (chrom_idx, start)

    selected_peaks.sort(key=sort_key)

    with open(BED_OUTPUT, 'w') as f:
        for peak in selected_peaks:
            f.write(peak)

    print(f"  Selected {len(selected_peaks)} peaks across {len(chrom_counts)} chromosomes")
    for chrom in CHROMOSOMES:
        if chrom in chrom_counts:
            print(f"    {chrom}: {chrom_counts[chrom]} peaks")

    return selected_peaks


def get_peak_ranges(peaks):
    """Get coordinate ranges for each chromosome from selected peaks."""
    chrom_ranges = {}

    for line in peaks:
        fields = line.strip().split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])

        if chrom not in chrom_ranges:
            chrom_ranges[chrom] = (start, end)
        else:
            cur_start, cur_end = chrom_ranges[chrom]
            chrom_ranges[chrom] = (min(cur_start, start), max(cur_end, end))

    # Add buffer
    for chrom in chrom_ranges:
        start, end = chrom_ranges[chrom]
        chrom_ranges[chrom] = (max(0, start - BUFFER_SIZE), end + BUFFER_SIZE)

    return chrom_ranges


def create_subset_gtf(peak_ranges):
    """Create subset GTF file with annotations overlapping peak ranges."""
    print(f"Creating subset GTF file: {GTF_OUTPUT}")

    header_lines = []
    selected_lines = []

    with open(GTF_INPUT, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line)
                continue

            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue

            chrom = fields[0]
            if chrom not in peak_ranges:
                continue

            try:
                gtf_start = int(fields[3])
                gtf_end = int(fields[4])
            except ValueError:
                continue

            range_start, range_end = peak_ranges[chrom]

            # Check for overlap
            if gtf_end >= range_start and gtf_start <= range_end:
                selected_lines.append(line)

    with open(GTF_OUTPUT, 'w') as f:
        for line in header_lines:
            f.write(line)
        for line in selected_lines:
            f.write(line)

    # Count by chromosome
    chrom_counts = defaultdict(int)
    for line in selected_lines:
        fields = line.strip().split('\t')
        chrom_counts[fields[0]] += 1

    print(f"  Selected {len(selected_lines)} GTF entries across {len(chrom_counts)} chromosomes")
    for chrom in CHROMOSOMES:
        if chrom in chrom_counts:
            print(f"    {chrom}: {chrom_counts[chrom]} entries")


def main():
    print("=" * 60)
    print("Creating subset test data for CI integration tests")
    print("=" * 60)

    # Check inputs exist
    if not os.path.exists(BED_INPUT):
        print(f"ERROR: BED input not found: {BED_INPUT}")
        sys.exit(1)
    if not os.path.exists(GTF_INPUT):
        print(f"ERROR: GTF input not found: {GTF_INPUT}")
        sys.exit(1)

    # Create subset BED
    selected_peaks = create_subset_bed()

    # Get peak ranges per chromosome
    peak_ranges = get_peak_ranges(selected_peaks)
    print("\nPeak ranges per chromosome (with 50kb buffer):")
    for chrom in CHROMOSOMES:
        if chrom in peak_ranges:
            start, end = peak_ranges[chrom]
            print(f"  {chrom}: {start:,} - {end:,}")

    # Create subset GTF
    print()
    create_subset_gtf(peak_ranges)

    # Print file sizes
    print("\n" + "=" * 60)
    print("Output file sizes:")
    bed_size = os.path.getsize(BED_OUTPUT)
    gtf_size = os.path.getsize(GTF_OUTPUT)
    print(f"  {BED_OUTPUT}: {bed_size:,} bytes ({bed_size/1024:.1f} KB)")
    print(f"  {GTF_OUTPUT}: {gtf_size:,} bytes ({gtf_size/1024/1024:.1f} MB)")
    print("=" * 60)


if __name__ == "__main__":
    main()
