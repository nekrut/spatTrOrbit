#!/usr/bin/env python3
import sys
import traceback

def main():
    try:
        import argparse
        import os
        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        import anndata as ad
        import json

        parser = argparse.ArgumentParser(description="TISSUE Spatial Transcriptomics Analysis")
        parser.add_argument("--spatial-counts", required=True, help="Spatial count matrix")
        parser.add_argument("--locations", required=True, help="Spatial coordinates (X, Y)")
        parser.add_argument("--scrna-counts", required=True, help="scRNA-seq count matrix")
        parser.add_argument("--target-gene", default="acta2", help="Target gene to predict")
        parser.add_argument("--prediction-method", default="spage", help="Prediction method (spage)")
        parser.add_argument("--n-folds", type=int, default=10, help="CV folds")
        parser.add_argument("--n-pv", type=int, default=10, help="Number of principal vectors")
        parser.add_argument("--n-neighbors", type=int, default=15, help="Neighbors for spatial graph")
        parser.add_argument("--grouping-method", default="kmeans_gene_cell", help="Stratified grouping method")
        parser.add_argument("--k", type=int, default=4, help="Gene groups")
        parser.add_argument("--k2", type=int, default=2, help="Cell groups")
        parser.add_argument("--alpha-level", type=float, default=0.23, help="Alpha level for prediction interval")
        parser.add_argument("--seed", type=int, default=42, help="Random seed")
        parser.add_argument("--out-prefix", default="galaxy_out", help="Prefix for outputs")
        args = parser.parse_args()

        np.random.seed(args.seed)

        import tissue.main

        print("1. Loading paired datasets...")
        adata, RNAseq_adata = tissue.main.load_paired_datasets(
            args.spatial_counts, args.locations, args.scrna_counts
        )

        adata.var_names = [x.lower() for x in adata.var_names]
        RNAseq_adata.var_names = [x.lower() for x in RNAseq_adata.var_names]

        print("2. Preprocessing scRNA...")
        tissue.main.preprocess_data(RNAseq_adata, standardize=False, normalize=True)

        gene_names = np.intersect1d(adata.var_names, RNAseq_adata.var_names)
        adata = adata[:, gene_names].copy()

        target_gene_lower = args.target_gene.lower()
        if target_gene_lower not in adata.var_names:
            raise ValueError(f"Target gene '{args.target_gene}' not found in spatial data.")

        target_expn = adata[:, target_gene_lower].X.copy()
        adata = adata[:, [g for g in gene_names if g != target_gene_lower]].copy()

        print(f"3. Predicting gene expression for {target_gene_lower} using {args.prediction_method}...")
        tissue.main.predict_gene_expression(
            adata, RNAseq_adata, [target_gene_lower],
            method=args.prediction_method, n_folds=args.n_folds, n_pv=args.n_pv
        )

        print("4. Building spatial graph...")
        tissue.main.build_spatial_graph(adata, method="fixed_radius", n_neighbors=args.n_neighbors)

        print("5. Conformalizing spatial uncertainty...")
        pred_key = f"{args.prediction_method}_predicted_expression"
        tissue.main.conformalize_spatial_uncertainty(
            adata, pred_key, calib_genes=adata.var_names,
            grouping_method=args.grouping_method, k=args.k, k2=args.k2
        )

        print(f"6. Conformalizing prediction interval (alpha={args.alpha_level})...")
        tissue.main.conformalize_prediction_interval(
            adata, pred_key, calib_genes=adata.var_names, alpha_level=args.alpha_level
        )

        print("7. Extracting intervals and generating plots...")
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
        fig_path = f"{args.out_prefix}_pi_width.png"
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
        pred_fig_path = f"{args.out_prefix}_prediction_vs_actual.png"
        plt.savefig(pred_fig_path, dpi=300, bbox_inches="tight")
        plt.close()

        print("8. Saving CSV and JSON...")
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
        csv_path = f"{args.out_prefix}_per_cell_results.csv"
        res_df.to_csv(csv_path, index=False)

        in_interval = (target_expn.flatten() >= lo_vals) & (target_expn.flatten() <= hi_vals)
        coverage = np.mean(in_interval)
        mae = np.mean(np.abs(target_expn.flatten() - pred_vals))
        corr = np.corrcoef(target_expn.flatten(), pred_vals)[0, 1]

        summary = {
            "target_gene": target_gene_lower,
            "prediction_method": args.prediction_method,
            "alpha_level": args.alpha_level,
            "expected_coverage": 1.0 - args.alpha_level,
            "empirical_coverage": float(coverage),
            "mae": float(mae),
            "pearson_r": float(corr),
            "mean_pi_width_log1p": float(np.mean(pi_width_log)),
            "min_pi_width_log1p": float(np.min(pi_width_log)),
            "max_pi_width_log1p": float(np.max(pi_width_log)),
            "p95_pi_width_log1p": float(p95),
            "num_cells": len(res_df),
            "random_seed": args.seed
        }

        with open(f"{args.out_prefix}_summary.json", "w") as f:
            json.dump(summary, f, indent=2)

        print("9. Pipeline finished successfully!")

    except Exception as e:
        err_msg = traceback.format_exc()
        sys.stderr.write(f"FATAL ERROR:\n{err_msg}\n")
        # Also write to summary.json so we can see it in Galaxy dataset!
        with open("galaxy_out_summary.json", "w") as f:
            f.write(json.dumps({"error": str(e), "traceback": err_msg}, indent=2))
        sys.exit(1)

if __name__ == "__main__":
    main()
