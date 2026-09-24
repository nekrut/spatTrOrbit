import subprocess
import json
import pandas as pd
import numpy as np

seeds = [1, 2, 3, 4, 5]
results = []

for s in seeds:
    print(f"Running replicate with seed {s}...")
    cmd = [
        ".loom/env/bin/python", "scripts/run_tissue_wrapper.py",
        "--spatial-counts", "data/benchmark_tutorial/Spatial_count.txt",
        "--locations", "data/benchmark_tutorial/Locations.txt",
        "--scrna-counts", "data/benchmark_tutorial/scRNA_count.txt",
        "--target-gene", "acta2",
        "--seed", str(s),
        "--out-prefix", f"local_reference/rep_tutorial_seed{s}"
    ]
    subprocess.run(cmd, check=True)
    with open(f"local_reference/rep_tutorial_seed{s}_summary.json") as f:
        res = json.load(f)
    results.append(res)

summary_df = pd.DataFrame(results)
summary_df.to_csv("local_reference/five_runs_summary.csv", index=False)
print("Five runs completed!")
print(summary_df[["random_seed", "empirical_coverage", "mae", "pearson_r", "mean_pi_width_log1p", "p95_pi_width_log1p"]])
