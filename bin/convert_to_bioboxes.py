#!/usr/bin/env python3
"""
Convert TAXPASTA standardized TSV output to CAMI Bioboxes profiling format.

TAXPASTA produces TSV with columns: taxonomy_id, sample1_count, sample2_count, ...
This script converts those counts to relative abundances in the CAMI Bioboxes
profiling format required by OPAL for evaluation.

Requires taxonkit or the NCBI taxonomy dump for lineage resolution.
Uses the names.dmp and nodes.dmp from NCBI taxonomy.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict


VALID_RANKS = [
    "superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"
]


def load_taxonomy(taxdump_dir):
    """Load NCBI taxonomy from names.dmp and nodes.dmp files."""
    names = {}
    nodes = {}  # taxid -> (parent_taxid, rank)

    nodes_file = os.path.join(taxdump_dir, "nodes.dmp")
    names_file = os.path.join(taxdump_dir, "names.dmp")

    if not os.path.exists(nodes_file) or not os.path.exists(names_file):
        print(f"ERROR: Cannot find nodes.dmp or names.dmp in {taxdump_dir}", file=sys.stderr)
        sys.exit(1)

    # Load nodes
    with open(nodes_file, "r") as f:
        for line in f:
            parts = line.split("|")
            taxid = int(parts[0].strip())
            parent = int(parts[1].strip())
            rank = parts[2].strip()
            nodes[taxid] = (parent, rank)

    # Load names (scientific names only)
    with open(names_file, "r") as f:
        for line in f:
            parts = line.split("|")
            taxid = int(parts[0].strip())
            name = parts[1].strip()
            name_class = parts[3].strip()
            if name_class == "scientific name":
                names[taxid] = name

    return names, nodes


def get_lineage(taxid, nodes, names):
    """Get the full lineage for a taxid as lists of (taxid, rank, name)."""
    lineage = []
    current = taxid
    visited = set()

    while current in nodes and current not in visited:
        visited.add(current)
        parent, rank = nodes[current]
        name = names.get(current, "")
        lineage.append((current, rank, name))
        if current == parent:  # root
            break
        current = parent

    lineage.reverse()
    return lineage


def build_taxpath(lineage):
    """Build TAXPATH and TAXPATHSN strings from lineage."""
    taxpath_ids = []
    taxpath_names = []

    for taxid, rank, name in lineage:
        if rank in VALID_RANKS:
            taxpath_ids.append(str(taxid))
            taxpath_names.append(name)

    return "|".join(taxpath_ids), "|".join(taxpath_names)


def get_rank_for_taxid(taxid, nodes):
    """Get the rank for a given taxid."""
    if taxid in nodes:
        return nodes[taxid][1]
    return "no rank"


def convert_taxpasta_to_bioboxes(input_tsv, sample_id, taxdump_dir, output_file):
    """Convert TAXPASTA TSV to CAMI Bioboxes profiling format."""
    names, nodes = load_taxonomy(taxdump_dir)

    # Read TAXPASTA output
    taxa_counts = {}
    total = 0

    with open(input_tsv, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)

        # Find the sample column (first non-taxonomy_id column, or specific sample)
        sample_col_idx = 1  # default: second column
        for i, col in enumerate(header):
            if col == sample_id or (i > 0 and sample_col_idx == 1):
                sample_col_idx = i
                if col == sample_id:
                    break

        for row in reader:
            if len(row) <= sample_col_idx:
                continue
            try:
                taxid = int(row[0])
                count = float(row[sample_col_idx])
            except (ValueError, IndexError):
                continue

            if count > 0:
                taxa_counts[taxid] = count
                total += count

    if total == 0:
        # Write empty profile
        with open(output_file, "w") as out:
            out.write(f"@Version:0.9.1\n")
            out.write(f"@SampleID:{sample_id}\n")
            out.write(f"@Ranks:{'|'.join(VALID_RANKS)}\n")
            out.write(f"\n@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")
        return

    # Build profiles with percentages, aggregating at each rank
    rank_profiles = defaultdict(float)  # (taxid, rank) -> percentage

    for taxid, count in taxa_counts.items():
        percentage = (count / total) * 100.0
        lineage = get_lineage(taxid, nodes, names)

        # Add percentage at each rank level in the lineage
        for lin_taxid, lin_rank, _ in lineage:
            if lin_rank in VALID_RANKS:
                rank_profiles[(lin_taxid, lin_rank)] += percentage

    # Write CAMI Bioboxes format
    with open(output_file, "w") as out:
        out.write(f"@Version:0.9.1\n")
        out.write(f"@SampleID:{sample_id}\n")
        out.write(f"@Ranks:{'|'.join(VALID_RANKS)}\n")
        out.write(f"\n@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")

        # Sort by rank order, then by taxid
        rank_order = {r: i for i, r in enumerate(VALID_RANKS)}
        sorted_entries = sorted(
            rank_profiles.items(),
            key=lambda x: (rank_order.get(x[0][1], 99), x[0][0])
        )

        for (taxid, rank), percentage in sorted_entries:
            lineage = get_lineage(taxid, nodes, names)
            taxpath, taxpathsn = build_taxpath(lineage)

            # Only include up to current rank level
            rank_idx = VALID_RANKS.index(rank) if rank in VALID_RANKS else -1
            if rank_idx >= 0:
                lineage_filtered = [
                    (t, r, n) for t, r, n in lineage if r in VALID_RANKS[:rank_idx + 1]
                ]
                taxpath = "|".join(str(t) for t, _, _ in lineage_filtered)
                taxpathsn = "|".join(n for _, _, n in lineage_filtered)

            out.write(f"{taxid}\t{rank}\t{taxpath}\t{taxpathsn}\t{percentage:.5f}\n")

    print(f"Converted {len(taxa_counts)} taxa to Bioboxes format for sample {sample_id}")


def main():
    parser = argparse.ArgumentParser(
        description="Convert TAXPASTA TSV to CAMI Bioboxes profiling format"
    )
    parser.add_argument(
        "--input", required=True,
        help="Path to TAXPASTA standardized TSV file"
    )
    parser.add_argument(
        "--sample-id", required=True,
        help="Sample ID for the profile header"
    )
    parser.add_argument(
        "--taxdump", required=True,
        help="Path to NCBI taxonomy dump directory (containing names.dmp, nodes.dmp)"
    )
    parser.add_argument(
        "--output", default="profile.bioboxes",
        help="Output file path in CAMI Bioboxes format"
    )
    args = parser.parse_args()

    convert_taxpasta_to_bioboxes(args.input, args.sample_id, args.taxdump, args.output)


if __name__ == "__main__":
    main()
