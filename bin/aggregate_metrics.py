#!/usr/bin/env python3
"""
Aggregate per-tool metrics into a single summary for Stimulus optimization.

Takes multiple metrics.json files (one per profiling tool/database) and
produces a single aggregated metrics.json. Uses the best species-level F1
across all tools as the objective metric.
"""

import argparse
import json
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate per-tool benchmark metrics"
    )
    parser.add_argument(
        "metrics_files", nargs="+",
        help="One or more metrics.json files to aggregate"
    )
    parser.add_argument(
        "--output", default="metrics.json",
        help="Output aggregated metrics JSON file"
    )
    args = parser.parse_args()

    all_metrics = []
    for fpath in args.metrics_files:
        try:
            with open(fpath) as f:
                all_metrics.append(json.load(f))
        except Exception as e:
            print(f"WARNING: Could not read {fpath}: {e}", file=sys.stderr)

    if not all_metrics:
        aggregated = {
            "per_tool_metrics": [],
            "summary": {
                "weighted_f1": 0.0,
                "species_precision": 0.0,
                "species_recall": 0.0,
                "species_f1": 0.0,
                "objective_rank": "error",
                "aggregation": "none",
            },
            "error": "No valid metrics files found",
        }
    else:
        # Find best species F1 across tools
        best = max(all_metrics, key=lambda m: m.get("summary", {}).get("species_f1", 0.0))
        best_summary = best.get("summary", {})

        # Compute average metrics
        n = len(all_metrics)
        avg_f1 = sum(m.get("summary", {}).get("species_f1", 0.0) for m in all_metrics) / n
        avg_precision = sum(m.get("summary", {}).get("species_precision", 0.0) for m in all_metrics) / n
        avg_recall = sum(m.get("summary", {}).get("species_recall", 0.0) for m in all_metrics) / n

        tool_names = [m.get("sample", "unknown") for m in all_metrics]
        tool_combinations = list(set(m.get("tool_combination", "") for m in all_metrics))

        aggregated = {
            "tools_evaluated": tool_names,
            "tool_combination": ",".join(tool_combinations),
            "num_profiles": n,
            "per_tool_metrics": [
                {
                    "sample": m.get("sample", "unknown"),
                    "species_f1": m.get("summary", {}).get("species_f1", 0.0),
                    "species_precision": m.get("summary", {}).get("species_precision", 0.0),
                    "species_recall": m.get("summary", {}).get("species_recall", 0.0),
                }
                for m in all_metrics
            ],
            "summary": {
                "weighted_f1": round(best_summary.get("species_f1", 0.0), 6),
                "species_precision": round(best_summary.get("species_precision", 0.0), 6),
                "species_recall": round(best_summary.get("species_recall", 0.0), 6),
                "species_f1": round(best_summary.get("species_f1", 0.0), 6),
                "avg_f1": round(avg_f1, 6),
                "avg_precision": round(avg_precision, 6),
                "avg_recall": round(avg_recall, 6),
                "objective_rank": best_summary.get("objective_rank", "species"),
                "aggregation": "best_of_n",
                "best_tool": best.get("sample", "unknown"),
            },
        }

    with open(args.output, "w") as f:
        json.dump(aggregated, f, indent=2)

    s = aggregated["summary"]
    print(
        f"Aggregated metrics: weighted_f1={s['weighted_f1']}, "
        f"n_tools={aggregated.get('num_profiles', 0)}"
    )


if __name__ == "__main__":
    main()
