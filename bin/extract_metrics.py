#!/usr/bin/env python3
"""
Extract benchmarking metrics from OPAL evaluation output.

Parses the OPAL results.tsv file and produces a JSON metrics file
with per-rank breakdown and a weighted F1 score for use with the
stimulus optimization framework.
"""

import argparse
import csv
import json
import os
import sys


def parse_opal_results(results_dir):
    """Parse OPAL output directory for metrics.

    OPAL produces results.tsv with columns including:
    tool, rank, metric, value
    """
    results_tsv = os.path.join(results_dir, "results.tsv")

    if not os.path.exists(results_tsv):
        # Try alternative location
        for f in os.listdir(results_dir):
            if f.endswith("results.tsv"):
                results_tsv = os.path.join(results_dir, f)
                break

    metrics = {}

    if not os.path.exists(results_tsv):
        return metrics

    with open(results_tsv, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            tool = row.get("tool", "")
            rank = row.get("rank", "")
            metric = row.get("metric", "")
            value_str = row.get("value", "")

            try:
                value = float(value_str) if value_str and value_str != "nan" else 0.0
            except ValueError:
                value = 0.0

            if rank not in metrics:
                metrics[rank] = {}
            metrics[rank][metric] = value

    return metrics


def parse_simple_metrics(metrics_file):
    """Parse a simple TSV metrics file with columns: metric, value."""
    metrics = {}
    with open(metrics_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                try:
                    metrics[parts[0]] = float(parts[1])
                except ValueError:
                    metrics[parts[0]] = parts[1]
    return metrics


def compute_objective(metrics):
    """Compute the objective metric from OPAL results.

    Strategy: Use species-level F1 score as the primary objective.
    If species F1 is not available, fall back to genus-level, then overall average.
    """
    # Prefer species-level F1
    for rank in ["species", "genus", "family", "order", "class", "phylum"]:
        if rank in metrics:
            rank_metrics = metrics[rank]
            f1 = rank_metrics.get("f1_score", rank_metrics.get("F1 score", 0.0))
            if f1 and f1 > 0:
                return f1, rank

    # Fallback: average F1 across all ranks
    f1_values = []
    for rank, rank_metrics in metrics.items():
        f1 = rank_metrics.get("f1_score", rank_metrics.get("F1 score", 0.0))
        if f1 and isinstance(f1, (int, float)) and f1 > 0:
            f1_values.append(f1)

    if f1_values:
        return sum(f1_values) / len(f1_values), "average"

    return 0.0, "none"


def build_metrics(sample_id, tool_combination, opal_metrics, precision_val,
                  recall_val, f1_val, objective_rank):
    """Build the full metrics JSON structure."""

    # Extract per-rank metrics
    per_rank = {}
    for rank, rank_metrics in opal_metrics.items():
        per_rank[rank] = {
            "precision": rank_metrics.get("purity (precision)", rank_metrics.get("precision", 0.0)),
            "recall": rank_metrics.get("completeness (recall)", rank_metrics.get("recall", 0.0)),
            "f1_score": rank_metrics.get("f1_score", rank_metrics.get("F1 score", 0.0)),
            "l1_norm": rank_metrics.get("L1 norm error", rank_metrics.get("l1_norm", 0.0)),
            "weighted_unifrac": rank_metrics.get("Weighted UniFrac error",
                                                  rank_metrics.get("weighted_unifrac", 0.0)),
        }

    return {
        "sample": sample_id,
        "tool_combination": tool_combination,
        "per_rank_metrics": per_rank,
        "summary": {
            "weighted_f1": round(f1_val, 6),
            "species_precision": round(precision_val, 6),
            "species_recall": round(recall_val, 6),
            "species_f1": round(f1_val, 6),
            "objective_rank": objective_rank,
        },
    }


def main():
    parser = argparse.ArgumentParser(
        description="Extract benchmarking metrics from OPAL evaluation output"
    )
    parser.add_argument(
        "--opal-results", required=False,
        help="Path to OPAL results directory"
    )
    parser.add_argument(
        "--metrics-tsv", required=False,
        help="Path to simple metrics TSV file (alternative to OPAL)"
    )
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument(
        "--tools", required=True,
        help="Comma-separated list of profiling tools used"
    )
    parser.add_argument(
        "--output", default="metrics.json", help="Output JSON file path"
    )
    args = parser.parse_args()

    try:
        if args.opal_results:
            opal_metrics = parse_opal_results(args.opal_results)
        elif args.metrics_tsv:
            raw = parse_simple_metrics(args.metrics_tsv)
            opal_metrics = {"species": raw}
        else:
            opal_metrics = {}

        f1_val, objective_rank = compute_objective(opal_metrics)

        # Get species-level precision/recall if available
        species_metrics = opal_metrics.get("species", opal_metrics.get(objective_rank, {}))
        precision_val = species_metrics.get("purity (precision)",
                                             species_metrics.get("precision", 0.0))
        recall_val = species_metrics.get("completeness (recall)",
                                          species_metrics.get("recall", 0.0))

        if not isinstance(precision_val, (int, float)):
            precision_val = 0.0
        if not isinstance(recall_val, (int, float)):
            recall_val = 0.0

        metrics = build_metrics(
            args.sample, args.tools, opal_metrics,
            precision_val, recall_val, f1_val, objective_rank
        )

    except Exception as e:
        # On failure, produce degraded metrics so optimizer gets a valid but bad score
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
        print(f"WARNING: Failed to extract metrics: {e}", file=sys.stderr)

    with open(args.output, "w") as f:
        json.dump(metrics, f, indent=2)

    print(
        f"Metrics: weighted_f1={metrics['summary']['weighted_f1']}, "
        f"species_f1={metrics['summary']['species_f1']}, "
        f"precision={metrics['summary']['species_precision']}, "
        f"recall={metrics['summary']['species_recall']}"
    )


if __name__ == "__main__":
    main()
