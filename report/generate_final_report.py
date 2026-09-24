import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt

def main():
    print("Compiling final reproduction deliverables...")
    # Load summaries
    with open("local_reference/full_dataset15_galaxy_summary.json") as f:
        gal_full_summary = json.load(f)
    with open("local_reference/full_dataset15_acta2_summary.json") as f:
        loc_full_summary = json.load(f)
    with open("local_reference/galaxy_smoke_summary.json") as f:
        gal_tutorial_summary = json.load(f)
    with open("local_reference/smoke_acta2_sklearn15_summary.json") as f:
        loc_tutorial_summary = json.load(f)

    # Comparison metrics table
    metrics_data = [
        {
            "Dataset": "Dataset 15 (Tutorial subset)",
            "Environment": "Local Reference",
            "Empirical Coverage": f"{loc_tutorial_summary['empirical_coverage']:.4f}",
            "Expected Coverage": f"{loc_tutorial_summary['expected_coverage']:.2f}",
            "MAE": f"{loc_tutorial_summary['mae']:.4f}",
            "Pearson r": f"{loc_tutorial_summary['pearson_r']:.4f}",
            "Mean PI Width (log1p)": f"{loc_tutorial_summary['mean_pi_width_log1p']:.4f}",
            "P95 PI Width (log1p)": f"{loc_tutorial_summary['p95_pi_width_log1p']:.4f}"
        },
        {
            "Dataset": "Dataset 15 (Tutorial subset)",
            "Environment": "Galaxy UDT",
            "Empirical Coverage": f"{gal_tutorial_summary['empirical_coverage']:.4f}",
            "Expected Coverage": f"{gal_tutorial_summary['expected_coverage']:.2f}",
            "MAE": f"{gal_tutorial_summary['mae']:.4f}",
            "Pearson r": f"{gal_tutorial_summary['pearson_r']:.4f}",
            "Mean PI Width (log1p)": f"{gal_tutorial_summary['mean_pi_width_log1p']:.4f}",
            "P95 PI Width (log1p)": f"{gal_tutorial_summary['p95_pi_width_log1p']:.4f}"
        },
        {
            "Dataset": "Full Dataset 15 (Full scRNA)",
            "Environment": "Local Reference",
            "Empirical Coverage": f"{loc_full_summary['empirical_coverage']:.4f}",
            "Expected Coverage": f"{loc_full_summary['expected_coverage']:.2f}",
            "MAE": f"{loc_full_summary['mae']:.4f}",
            "Pearson r": f"{loc_full_summary['pearson_r']:.4f}",
            "Mean PI Width (log1p)": f"{loc_full_summary['mean_pi_width_log1p']:.4f}",
            "P95 PI Width (log1p)": f"{loc_full_summary['p95_pi_width_log1p']:.4f}"
        },
        {
            "Dataset": "Full Dataset 15 (Full scRNA)",
            "Environment": "Galaxy UDT",
            "Empirical Coverage": f"{gal_full_summary['empirical_coverage']:.4f}",
            "Expected Coverage": f"{gal_full_summary['expected_coverage']:.2f}",
            "MAE": f"{gal_full_summary['mae']:.4f}",
            "Pearson r": f"{gal_full_summary['pearson_r']:.4f}",
            "Mean PI Width (log1p)": f"{gal_full_summary['mean_pi_width_log1p']:.4f}",
            "P95 PI Width (log1p)": f"{gal_full_summary['p95_pi_width_log1p']:.4f}"
        }
    ]
    df_metrics = pd.DataFrame(metrics_data)
    df_metrics.to_csv("report/metrics_comparison.csv", index=False)
    print("Saved report/metrics_comparison.csv")

if __name__ == "__main__":
    main()
