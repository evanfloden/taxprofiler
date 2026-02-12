#!/usr/bin/env python3
"""
Evaluate a taxonomic profile against a gold standard.

Computes precision, recall, F1, and other metrics per taxonomic rank
by comparing a predicted profile against a truth profile, both in
CAMI Bioboxes format.

Produces a metrics.json file for the stimulus optimization framework.
"""

import argparse
import json
import sys
from collections import defaultdict


def parse_bioboxes_profile(filepath):
    """Parse a CAMI Bioboxes profiling format file.

    Returns:
        dict: {rank: {taxid: percentage}}
    """
    profiles = defaultdict(dict)  # rank -> {taxid -> percentage}
    sample_id = None

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()

            # Parse header
            if line.startswith("@SampleID:"):
                sample_id = line.split(":", 1)[1].strip()
                continue
            if line.startswith("@") or line.startswith("#"):
                continue
            if line.startswith("@@"):
                continue
            if not line:
                continue

            # Parse data line
            parts = line.split("\t")
            if len(parts) < 5:
                continue

            try:
                taxid = int(parts[0])
                rank = parts[1]
                percentage = float(parts[4])
            except (ValueError, IndexError):
                continue

            if percentage > 0:
                profiles[rank][taxid] = percentage

    return profiles, sample_id


def compute_metrics(truth_profile, pred_profile, rank):
    """Compute precision, recall, F1 for a given taxonomic rank.

    Uses presence/absence metrics (taxa with >0 abundance).
    """
    truth_taxa = set(truth_profile.get(rank, {}).keys())
    pred_taxa = set(pred_profile.get(rank, {}).keys())

    tp = len(truth_taxa & pred_taxa)
    fp = len(pred_taxa - truth_taxa)
    fn = len(truth_taxa - pred_taxa)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    # Abundance-based L1 norm
    all_taxa = truth_taxa | pred_taxa
    l1_norm = 0.0
    for taxid in all_taxa:
        truth_pct = truth_profile.get(rank, {}).get(taxid, 0.0)
        pred_pct = pred_profile.get(rank, {}).get(taxid, 0.0)
        l1_norm += abs(truth_pct - pred_pct)

    # Bray-Curtis
    sum_truth = sum(truth_profile.get(rank, {}).values())
    sum_pred = sum(pred_profile.get(rank, {}).values())
    bray_curtis = l1_norm / (sum_truth + sum_pred) if (sum_truth + sum_pred) > 0 else 1.0

    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "truth_total": len(truth_taxa),
        "pred_total": len(pred_taxa),
        "precision": round(precision, 6),
        "recall": round(recall, 6),
        "f1_score": round(f1, 6),
        "l1_norm": round(l1_norm, 6),
        "bray_curtis": round(bray_curtis, 6),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate taxonomic profile against gold standard"
    )
    parser.add_argument(
        "--prediction", required=True,
        help="Path to predicted profile in CAMI Bioboxes format"
    )
    parser.add_argument(
        "--truth", required=True,
        help="Path to gold standard profile in CAMI Bioboxes format"
    )
    parser.add_argument(
        "--sample", required=True,
        help="Sample name for output"
    )
    parser.add_argument(
        "--tools", required=True,
        help="Comma-separated list of profiling tools used"
    )
    parser.add_argument(
        "--output", default="metrics.json",
        help="Output JSON file path"
    )
    args = parser.parse_args()

    try:
        truth_profile, truth_sample = parse_bioboxes_profile(args.truth)
        pred_profile, pred_sample = parse_bioboxes_profile(args.prediction)

        # Compute metrics per rank
        ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"]
        per_rank = {}
        for rank in ranks:
            if rank in truth_profile or rank in pred_profile:
                per_rank[rank] = compute_metrics(truth_profile, pred_profile, rank)

        # Determine objective: species-level F1, falling back to genus
        objective_rank = "species"
        f1_val = 0.0
        for rank in ["species", "genus", "family", "order", "class", "phylum"]:
            if rank in per_rank and per_rank[rank]["f1_score"] > 0:
                f1_val = per_rank[rank]["f1_score"]
                objective_rank = rank
                break
            elif rank in per_rank:
                f1_val = per_rank[rank]["f1_score"]
                objective_rank = rank
                break

        species = per_rank.get("species", per_rank.get(objective_rank, {}))
        precision_val = species.get("precision", 0.0)
        recall_val = species.get("recall", 0.0)

        metrics = {
            "sample": args.sample,
            "tool_combination": args.tools,
            "per_rank_metrics": per_rank,
            "summary": {
                "weighted_f1": round(f1_val, 6),
                "species_precision": round(precision_val, 6),
                "species_recall": round(recall_val, 6),
                "species_f1": round(f1_val, 6),
                "objective_rank": objective_rank,
            },
        }

    except Exception as e:
        metrics = {
            "sample": args.sample,
            "tool_combination": args.tools,
            "per_rank_metrics": {},
            "summary": {
                "weighted_f1": 0.0,
                "species_precision": 0.0,
                "species_recall": 0.0,
                "species_f1": 0.0,
                "objective_rank": "error",
            },
            "error": str(e),
        }
        print(f"WARNING: Failed to evaluate profile: {e}", file=sys.stderr)

    with open(args.output, "w") as f:
        json.dump(metrics, f, indent=2)

    print(
        f"Metrics: weighted_f1={metrics['summary']['weighted_f1']}, "
        f"precision={metrics['summary']['species_precision']}, "
        f"recall={metrics['summary']['species_recall']}, "
        f"objective_rank={metrics['summary']['objective_rank']}"
    )


if __name__ == "__main__":
    main()
