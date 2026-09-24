#!/usr/bin/env python3
"""
TISSUE reproduction pipeline for Acta2 prediction and prediction intervals.
Executes data loading, preprocessing, SpaGE spatial prediction, spatial graph construction,
conformal uncertainty calibration, prediction interval extraction, and figure generation.
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad

def run_pipeline(
    spatial_count_path: str,
    locations_path: str,
    scrna_count_path: str,
    target_gene: str = "acta2",
    prediction_method: str = "spage",
    n_folds: int = 10,
    n_pv: int = 10,
    n_neighbors: int = 15,
    grouping_method: str = "kmeans_gene_cell",
    k: int = 4,
    k2: int = 2,
    alpha_level: float = 0.23,
    seed: int = 42,
    output_prefix: str = "tissue_acta2"
):
    np.random.seed(seed)
    
    # Import TISSUE
    import tissue.main
    
    print(f"Loading paired datasets...")
    adata, RNAseq_adata = tissue.main.load_paired_datasets(
        spatial_count_path, locations_path, scrna_count_path
    )
    
    # Case normalization
    adata.var_names = [x.lower() for x in adata.var_names]
    RNAseq_adata.var_names = [x.lower() for x in RNAseq_adata.var_names]
    
    # Preprocess scRNA-seq
    tissue.main.preprocess_data(RNAseq_adata, standardize=False, normalize=True)
    
    # Shared genes
    gene_names = np.intersect1d(adata.var_names, RNAseq_adata.var_names)
    adata = adata[:, gene_names].copy()
    
    target_gene_lower = target_gene.lower()
    if target_gene_lower not in adata.var_names:
        raise ValueError(f"Target gene '{target_gene}' not found in spatial dataset.")
    
    # Hold out target gene
    target_expn = adata[:, target_gene_lower].X.copy()
    adata = adata[:, [g for g in gene_names if g != target_gene_lower]].copy()
    
    print(f"Predicting gene expression for {target_gene_lower} using {prediction_method}...")
    tissue.main.predict_gene_expression(
        adata, RNAseq_adata, [target_gene_lower],
        method=prediction_method, n_folds=n_folds, n_pv=n_pv
    )
    
    print("Building spatial graph...")
    tissue.main.build_spatial_graph(adata, method="fixed_radius", n_neighbors=n_neighbors)
    
    print("Conformalizing spatial uncertainty...")
    pred_key = f"{prediction_method}_predicted_expression"
    tissue.main.conformalize_spatial_uncertainty(
        adata, pred_key, calib_genes=adata.var_names,
        grouping_method=grouping_method, k=k, k2=k2
    )
    
    print(f"Conformalizing prediction interval (alpha={alpha_level})...")
    tissue.main.conformalize_prediction_interval(
        adata, pred_key, calib_genes=adata.var_names, alpha_level=alpha_level
    )
    
    # Compute PI width and log1p transform
    lo_vals = adata.obsm[f"{pred_key}_lo"][target_gene_lower].values
    hi_vals = adata.obsm[f"{pred_key}_hi"][target_gene_lower].values
    pi_width_raw = hi_vals - lo_vals
    pi_width_clean = pi_width_raw.copy()
    pi_width_clean[pi_width_clean < 0] = 0
    pi_width_log = np.log1p(pi_width_clean)
    p95 = np.percentile(pi_width_log, 95)
    
    # Plot matching Extended Data Figure 1C
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.axis("off")
    cmap_plot = pi_width_log.copy()
    cmap_plot[cmap_plot > p95] = p95
    im = ax.scatter(
        adata.obsm["spatial"][:, 0], adata.obsm["spatial"][:, 1],
        s=1, c=cmap_plot, rasterized=True
    )
    ax.set_title(f"PI Width {target_gene_lower}", fontsize=12)
    cbar = fig.colorbar(im, ax=ax)
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel("Log Expression", rotation=270)
    fig_path = f"{output_prefix}_pi_width.png"
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close()
    
    # Prediction vs Actual plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    ax1.axis("off")
    act_log = np.log1p(np.maximum(0, target_expn.flatten()))
    act_p95 = np.percentile(act_log, 95)
    act_log[act_log > act_p95] = act_p95
    im1 = ax1.scatter(adata.obsm["spatial"][:, 0], adata.obsm["spatial"][:, 1], s=1, c=act_log, rasterized=True)
    ax1.set_title(f"Actual {target_gene_lower}")
    cbar1 = fig.colorbar(im1, ax=ax1)
    cbar1.ax.get_yaxis().labelpad = 15
    cbar1.ax.set_ylabel("Log Expression", rotation=270)
    
    ax2.axis("off")
    pred_vals = adata.obsm[pred_key][target_gene_lower].values
    pred_log = np.log1p(np.maximum(0, pred_vals))
    pred_p95 = np.percentile(pred_log, 95)
    pred_log[pred_log > pred_p95] = pred_p95
    im2 = ax2.scatter(adata.obsm["spatial"][:, 0], adata.obsm["spatial"][:, 1], s=1, c=pred_log, rasterized=True)
    ax2.set_title(f"Predicted {target_gene_lower}")
    cbar2 = fig.colorbar(im2, ax=ax2)
    cbar2.ax.get_yaxis().labelpad = 15
    cbar2.ax.set_ylabel("Log Expression", rotation=270)
    pred_fig_path = f"{output_prefix}_prediction_vs_actual.png"
    plt.savefig(pred_fig_path, dpi=300, bbox_inches="tight")
    plt.close()
    
    # Save per-cell table
    res_df = pd.DataFrame({
        "cell_id": range(len(adata.obs)),
        "x_coord": adata.obsm["spatial"][:, 0],
        "y_coord": adata.obsm["spatial"][:, 1],
        "actual_expression": target_expn.flatten(),
        "predicted_expression": pred_vals,
        "pi_lower": lo_vals,
        "pi_upper": hi_vals,
        "pi_width_raw": pi_width_raw,
        "pi_width_log1p": pi_width_log
    })
    csv_path = f"{output_prefix}_per_cell_results.csv"
    res_df.to_csv(csv_path, index=False)
    
    # Save evaluation summary
    in_interval = (target_expn.flatten() >= lo_vals) & (target_expn.flatten() <= hi_vals)
    coverage = np.mean(in_interval)
    mae = np.mean(np.abs(target_expn.flatten() - pred_vals))
    corr = np.corrcoef(target_expn.flatten(), pred_vals)[0, 1]
    
    summary = {
        "target_gene": target_gene_lower,
        "prediction_method": prediction_method,
        "alpha_level": alpha_level,
        "expected_coverage": 1.0 - alpha_level,
        "empirical_coverage": float(coverage),
        "mae": float(mae),
        "pearson_r": float(corr),
        "mean_pi_width_log1p": float(np.mean(pi_width_log)),
        "min_pi_width_log1p": float(np.min(pi_width_log)),
        "max_pi_width_log1p": float(np.max(pi_width_log)),
        "p95_pi_width_log1p": float(p95),
        "num_cells": len(res_df),
        "random_seed": seed
    }
    
    import json
    with open(f"{output_prefix}_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
        
    print("Done! Outputs:")
    print(f"  Figure 1C: {fig_path}")
    print(f"  Prediction vs actual: {pred_fig_path}")
    print(f"  Per-cell table: {csv_path}")
    print(f"  Summary: {output_prefix}_summary.json")
    print(f"Empirical coverage: {coverage:.4f} (target: {1-alpha_level:.4f}), MAE: {mae:.4f}, r: {corr:.4f}")
    return summary

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--spatial-counts", required=True)
    parser.add_argument("--locations", required=True)
    parser.add_argument("--scrna-counts", required=True)
    parser.add_argument("--target-gene", default="acta2")
    parser.add_argument("--prediction-method", default="spage")
    parser.add_argument("--n-folds", type=int, default=10)
    parser.add_argument("--n-pv", type=int, default=10)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--grouping-method", default="kmeans_gene_cell")
    parser.add_argument("--k", type=int, default=4)
    parser.add_argument("--k2", type=int, default=2)
    parser.add_argument("--alpha-level", type=float, default=0.23)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out-prefix", default="tissue_acta2")
    args = parser.parse_args()
    
    run_pipeline(
        spatial_count_path=args.spatial_counts,
        locations_path=args.locations,
        scrna_count_path=args.scrna_counts,
        target_gene=args.target_gene,
        prediction_method=args.prediction_method,
        n_folds=args.n_folds,
        n_pv=args.n_pv,
        n_neighbors=args.n_neighbors,
        grouping_method=args.grouping_method,
        k=args.k,
        k2=args.k2,
        alpha_level=args.alpha_level,
        seed=args.seed,
        output_prefix=args.out_prefix
    )
