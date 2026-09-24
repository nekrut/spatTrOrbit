# spatTrOrbit

Artifacts from an **Orbit** session. In that session we reproduced the TISSUE spatial transcriptomics analysis (Miao et al., *Nature* 2026, Extended Data Fig. 1C) on Galaxy (usegalaxy.org). TISSUE is the uncertainty-calibrated method from Sun et al., *Nat Methods* 2024. The run predicts spatial *Acta2* expression and computes conformal prediction intervals for it.

The full analysis record, with methods, provenance, Galaxy history/workflow IDs and results, is in [`notebook.md`](notebook.md).

## Layout

- `notebook.md`: session lab notebook (main record)
- `report/`: final reproduction report (`.md`/`.pdf`), comparison figures, generating scripts
- `workflow/`: Galaxy workflows (`.ga`) and the scripts that build them
- `scripts/`: pipeline runners, replicate and full-evaluation drivers
- `local_reference/`: local reference runs (per-cell results, summaries, plots)
- `paper_assets/`: figures and supplementary info from the original paper
- `tissue_bundle/`: bundled TISSUE package used in the Galaxy tool
- `udt_test/`: Galaxy User-Defined Tool probes
- `data/`: input data (Dataset 15, tutorial benchmark, Zenodo record 8259942)
- `sources/`: upstream repos as submodules (`TISSUE`, `tissue-figures-and-analyses`, `tissue_mcp`)

## Notes

- Clone with `git clone --recurse-submodules` to get `sources/`.
- Two uncompressed `scRNA_count.txt` files (404 MB) exceed GitHub's file limit and are not committed. Use `data/full_dataset15/scRNA_count.txt.gz` or re-download them from Zenodo record 8259942.
