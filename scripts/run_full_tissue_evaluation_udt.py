#!/usr/bin/env python3
"""
Full TISSUE Dataset 15 Evaluation Pipeline (for Galaxy UDT):
1. Exact preprocessing semantics: full library-size normalization and log1p on reference.
2. Cross-validated held-out evaluation for all 33 measured spatial genes.
   - For each fold/gene, held-out target gene is completely excluded from predictors and calibration set.
3. Multi-coverage evaluation: assesses prediction intervals at nominal levels:
   77% (alpha=0.23), 80% (alpha=0.20), 90% (alpha=0.10), 95% (alpha=0.05).
   Reuses base prediction and conformal calibration scores.
4. Generates per-cell per-gene tables, per-gene summary metrics, calibration curves, and Acta2 maps.
"""

import argparse
import os
import sys
import gc
import gzip
import json
import traceback
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata as ad
from scipy.stats import pearsonr, spearmanr

def load_and_normalize_rnaseq(scrna_path, target_genes):
    needed_genes = set([g.lower() for g in target_genes])
    
    is_gz = False
    with open(scrna_path, "rb") as f:
        magic = f.read(2)
        if magic == b"\x1f\x8b":
            is_gz = True
            
    opener = gzip.open(scrna_path, "rt") if is_gz else open(scrna_path, "r")
    with opener as f:
        header = f.readline().strip().split("\t")
        cell_ids = header[1:] if header[0] == "" or header[0].lower() in ["gene", "genes"] else header
        n_cells = len(cell_ids)
        cell_totals = np.zeros(n_cells, dtype=np.float64)
        
        selected_rows = []
        selected_genes = []
        for line in f:
            first_tab = line.find("\t")
            if first_tab == -1: continue
            gene = line[:first_tab]
            vals = np.fromstring(line[first_tab+1:], sep="\t", dtype=np.float64)
            cell_totals += vals
            if gene.lower() in needed_genes:
                selected_rows.append(vals)
                selected_genes.append(gene.lower())
                
    if not selected_rows:
        raise ValueError(f"None of target genes found in {scrna_path}")
        
    mat = np.vstack(selected_rows)
    del selected_rows
    gc.collect()
    
    target_sum = np.median(cell_totals)
    normalized_mat = (mat / cell_totals[np.newaxis, :]) * target_sum
    log1p_mat = np.log1p(normalized_mat)
    
    RNAseq_adata = ad.AnnData(X=log1p_mat.T.astype(np.float64))
    RNAseq_adata.obs_names = cell_ids[:mat.shape[1]]
    RNAseq_adata.var_names = selected_genes
    return RNAseq_adata

def main():
    try:
        parser = argparse.ArgumentParser(description="Full TISSUE Evaluation on Dataset 15")
        parser.add_argument("--spatial-counts", required=True)
        parser.add_argument("--locations", required=True)
        parser.add_argument("--scrna-counts", required=True)
        parser.add_argument("--n-folds", type=int, default=10)
        parser.add_argument("--n-pv", type=int, default=10)
        parser.add_argument("--n-neighbors", type=int, default=15)
        parser.add_argument("--k", type=int, default=4)
        parser.add_argument("--k2", type=int, default=2)
        parser.add_argument("--seed", type=int, default=444)
        parser.add_argument("--out-prefix", default="galaxy_out")
        args = parser.parse_args()

        np.random.seed(args.seed)
        import tissue.main

        print("1. Loading spatial dataset...")
        spatial_adata = tissue.main.load_spatial_data(args.spatial_counts, args.locations)
        spatial_adata.var_names = [x.lower() for x in spatial_adata.var_names]
        measured_genes = list(spatial_adata.var_names)
        print(f"Spatial cells: {spatial_adata.shape[0]}, Measured genes: {len(measured_genes)}")

        print("2. Normalizing reference scRNA-seq with exact full-transcriptome library sizing...")
        RNAseq_adata = load_and_normalize_rnaseq(args.scrna_counts, measured_genes)
        print(f"Reference cells: {RNAseq_adata.shape[0]}, Intersecting genes: {RNAseq_adata.shape[1]}")

        shared_genes = np.intersect1d(spatial_adata.var_names, RNAseq_adata.var_names)
        spatial_adata = spatial_adata[:, shared_genes].copy()
        RNAseq_adata = RNAseq_adata[:, shared_genes].copy()
        all_genes = list(shared_genes)
        n_genes = len(all_genes)

        print(f"3. Setting up {args.n_folds}-fold cross validation across {n_genes} genes (seed={args.seed})...")
        np.random.seed(args.seed)
        shuffled_genes = all_genes.copy()
        np.random.shuffle(shuffled_genes)
        cv_folds = np.array_split(shuffled_genes, args.n_folds)

        n_cells = spatial_adata.shape[0]
        cv_predictions = pd.DataFrame(np.zeros((n_cells, n_genes)), index=spatial_adata.obs_names, columns=all_genes)
        
        nominal_alphas = [0.23, 0.20, 0.10, 0.05]
        pi_lo_dict = {a: pd.DataFrame(np.zeros((n_cells, n_genes)), index=spatial_adata.obs_names, columns=all_genes) for a in nominal_alphas}
        pi_hi_dict = {a: pd.DataFrame(np.zeros((n_cells, n_genes)), index=spatial_adata.obs_names, columns=all_genes) for a in nominal_alphas}

        print("4. Building spatial neighborhood graph...")
        tissue.main.build_spatial_graph(spatial_adata, method="fixed_radius", n_neighbors=args.n_neighbors)

        print("5. Executing cross-validation loops (leakage-free)...")
        for fi, fold in enumerate(cv_folds):
            target_fold_genes = list(fold)
            calib_genes = [g for g in all_genes if g not in target_fold_genes]
            print(f"  Fold {fi+1}/{args.n_folds}: Targets ({len(target_fold_genes)} genes): {target_fold_genes}")

            train_spatial = spatial_adata[:, calib_genes].copy()
            tissue.main.predict_gene_expression(
                train_spatial, RNAseq_adata, target_fold_genes,
                method="spage", n_folds=len(calib_genes), n_pv=args.n_pv, random_seed=args.seed
            )
            for tg in target_fold_genes:
                cv_predictions[tg] = train_spatial.obsm["spage_predicted_expression"][tg].values

            tissue.main.conformalize_spatial_uncertainty(
                train_spatial, "spage_predicted_expression", calib_genes=calib_genes,
                grouping_method="kmeans_gene_cell", k=args.k, k2=args.k2
            )

            for a in nominal_alphas:
                tissue.main.conformalize_prediction_interval(
                    train_spatial, "spage_predicted_expression", calib_genes=calib_genes, alpha_level=a
                )
                for tg in target_fold_genes:
                    pi_lo_dict[a][tg] = train_spatial.obsm["spage_predicted_expression_lo"][tg].values
                    pi_hi_dict[a][tg] = train_spatial.obsm["spage_predicted_expression_hi"][tg].values

        print("6. Calculating comprehensive evaluation metrics...")
        per_gene_metrics = []
        actual_mat = spatial_adata.X.toarray() if not isinstance(spatial_adata.X, np.ndarray) else spatial_adata.X

        for gi, gene in enumerate(all_genes):
            act_vals = actual_mat[:, gi]
            pred_vals = cv_predictions[gene].values
            mae = np.mean(np.abs(act_vals - pred_vals))
            p_corr, _ = pearsonr(act_vals, pred_vals)
            s_corr, _ = spearmanr(act_vals, pred_vals)

            gene_dict = {
                "gene": gene,
                "mae": float(mae),
                "pearson_r": float(p_corr) if not np.isnan(p_corr) else 0.0,
                "spearman_rho": float(s_corr) if not np.isnan(s_corr) else 0.0,
            }

            for a in nominal_alphas:
                lo = pi_lo_dict[a][gene].values
                hi = pi_hi_dict[a][gene].values
                cov = np.mean((act_vals >= lo) & (act_vals <= hi))
                w = np.maximum(0, hi - lo)
                w_log = np.log1p(w)
                nominal_pct = int(round((1 - a) * 100))
                gene_dict[f"cov_{nominal_pct}"] = float(cov)
                gene_dict[f"mean_width_log1p_{nominal_pct}"] = float(np.mean(w_log))
                gene_dict[f"p95_width_log1p_{nominal_pct}"] = float(np.percentile(w_log, 95))

            per_gene_metrics.append(gene_dict)

        df_gene_metrics = pd.DataFrame(per_gene_metrics)
        gene_metrics_csv = f"{args.out_prefix}_per_gene_metrics.csv"
        df_gene_metrics.to_csv(gene_metrics_csv, index=False)

        summary = {
            "dataset": "Dataset15_osmFISH_Mouse_SC",
            "n_spatial_cells": int(n_cells),
            "n_genes_evaluated": int(n_genes),
            "n_reference_cells": int(RNAseq_adata.shape[0]),
            "random_seed": args.seed,
            "mean_mae": float(df_gene_metrics["mae"].mean()),
            "mean_pearson_r": float(df_gene_metrics["pearson_r"].mean()),
            "mean_spearman_rho": float(df_gene_metrics["spearman_rho"].mean()),
            "mean_coverage_77": float(df_gene_metrics["cov_77"].mean()),
            "mean_coverage_80": float(df_gene_metrics["cov_80"].mean()),
            "mean_coverage_90": float(df_gene_metrics["cov_90"].mean()),
            "mean_coverage_95": float(df_gene_metrics["cov_95"].mean())
        }
        summary_json = f"{args.out_prefix}_overall_summary.json"
        with open(summary_json, "w") as f:
            json.dump(summary, f, indent=2)

        print("7. Generating calibration and comparison plots...")
        nominal_covs = [0.77, 0.80, 0.90, 0.95]
        observed_covs = [summary["mean_coverage_77"], summary["mean_coverage_80"],
                         summary["mean_coverage_90"], summary["mean_coverage_95"]]

        plt.figure(figsize=(6, 5))
        plt.plot([0.7, 1.0], [0.7, 1.0], "k--", label="Ideal calibration")
        plt.plot(nominal_covs, observed_covs, "o-", color="tab:blue", lw=2, ms=8, label="TISSUE observed coverage")
        plt.xlabel("Nominal Coverage Level (1 - alpha)", fontsize=11)
        plt.ylabel("Empirical Mean Coverage", fontsize=11)
        plt.title("TISSUE Calibration Curve across 33 Genes", fontsize=12)
        plt.grid(True, linestyle=":", alpha=0.6)
        plt.legend(loc="lower right")
        calib_plot_path = f"{args.out_prefix}_calibration_curve.png"
        plt.savefig(calib_plot_path, dpi=300, bbox_inches="tight")
        plt.close()

        if "acta2" in all_genes:
            acta2_lo = pi_lo_dict[0.23]["acta2"].values
            acta2_hi = pi_hi_dict[0.23]["acta2"].values
            acta2_w = np.maximum(0, acta2_hi - acta2_lo)
            acta2_w_log = np.log1p(acta2_w)
            p95 = np.percentile(acta2_w_log, 95)
            cmap_plot = acta2_w_log.copy()
            cmap_plot[cmap_plot > p95] = p95

            fig, ax = plt.subplots(figsize=(6, 5))
            ax.axis("off")
            im = ax.scatter(spatial_adata.obsm["spatial"][:, 0], spatial_adata.obsm["spatial"][:, 1],
                            s=1, c=cmap_plot, rasterized=True)
            ax.set_title("PI Width acta2 (Corrected Full Reference)", fontsize=12)
            cbar = fig.colorbar(im, ax=ax)
            cbar.ax.get_yaxis().labelpad = 15
            cbar.ax.set_ylabel("Log Expression", rotation=270)
            acta2_fig_path = f"{args.out_prefix}_acta2_pi_width.png"
            plt.savefig(acta2_fig_path, dpi=300, bbox_inches="tight")
            plt.close()

        per_cell_dfs = []
        for tg in all_genes:
            df_g = pd.DataFrame({
                "cell_id": range(n_cells),
                "gene": tg,
                "x_coord": spatial_adata.obsm["spatial"][:, 0],
                "y_coord": spatial_adata.obsm["spatial"][:, 1],
                "actual_expression": actual_mat[:, all_genes.index(tg)],
                "predicted_expression": cv_predictions[tg].values,
                "pi_lower_77": pi_lo_dict[0.23][tg].values,
                "pi_upper_77": pi_hi_dict[0.23][tg].values,
                "pi_lower_95": pi_lo_dict[0.05][tg].values,
                "pi_upper_95": pi_hi_dict[0.05][tg].values
            })
            per_cell_dfs.append(df_g)
        df_all_per_cell = pd.concat(per_cell_dfs, ignore_index=True)
        per_cell_csv = f"{args.out_prefix}_all_genes_per_cell.csv"
        df_all_per_cell.to_csv(per_cell_csv, index=False)

        print("Full evaluation completed successfully!")

    except Exception as e:
        err = traceback.format_exc()
        sys.stderr.write(f"FATAL ERROR:\n{err}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
