import json
import uuid

workflow_dict = {
    "a_galaxy_workflow": "true",
    "format-version": "0.1",
    "name": "TISSUE Full Cross-Validated Evaluation (Dataset 15 All Measured Genes)",
    "annotation": "Full cross-validated held-out spatial gene expression prediction and conformal uncertainty interval evaluation across all 33 measured genes on mouse somatosensory cortex with exact full reference normalization.",
    "tags": ["spatial-transcriptomics", "TISSUE", "calibration", "cross-validation", "Dataset-15"],
    "version": 1,
    "steps": {
        "0": {
            "annotation": "Spatial count matrix (Spatial_count.txt)",
            "content_id": None,
            "errors": None,
            "id": 0,
            "input_connections": {},
            "inputs": [{"description": "Spatial count matrix", "name": "Spatial counts"}],
            "label": "Spatial counts",
            "name": "Input dataset",
            "outputs": [],
            "position": {"left": 10, "top": 10},
            "tool_id": None,
            "tool_state": "{\"optional\": false, \"tag\": null}",
            "tool_version": None,
            "type": "data_input",
            "uuid": str(uuid.uuid4())
        },
        "1": {
            "annotation": "Spatial locations file (Locations.txt)",
            "content_id": None,
            "errors": None,
            "id": 1,
            "input_connections": {},
            "inputs": [{"description": "Spatial cell coordinates", "name": "Spatial locations"}],
            "label": "Spatial locations",
            "name": "Input dataset",
            "outputs": [],
            "position": {"left": 10, "top": 150},
            "tool_id": None,
            "tool_state": "{\"optional\": false, \"tag\": null}",
            "tool_version": None,
            "type": "data_input",
            "uuid": str(uuid.uuid4())
        },
        "2": {
            "annotation": "Full reference scRNA-seq counts (scRNA_count.txt.gz)",
            "content_id": None,
            "errors": None,
            "id": 2,
            "input_connections": {},
            "inputs": [{"description": "Reference scRNA-seq counts", "name": "scRNA-seq counts"}],
            "label": "scRNA-seq counts",
            "name": "Input dataset",
            "outputs": [],
            "position": {"left": 10, "top": 290},
            "tool_id": None,
            "tool_state": "{\"optional\": false, \"tag\": null}",
            "tool_version": None,
            "type": "data_input",
            "uuid": str(uuid.uuid4())
        },
        "3": {
            "annotation": "TISSUE Python package archive (tissue_pkg.tar.gz)",
            "content_id": None,
            "errors": None,
            "id": 3,
            "input_connections": {},
            "inputs": [{"description": "TISSUE core library bundle", "name": "TISSUE code package"}],
            "label": "TISSUE code package",
            "name": "Input dataset",
            "outputs": [],
            "position": {"left": 10, "top": 430},
            "tool_id": None,
            "tool_state": "{\"optional\": false, \"tag\": null}",
            "tool_version": None,
            "type": "data_input",
            "uuid": str(uuid.uuid4())
        },
        "4": {
            "annotation": "Full evaluation runner script (run_full_tissue_evaluation_udt.py)",
            "content_id": None,
            "errors": None,
            "id": 4,
            "input_connections": {},
            "inputs": [{"description": "Evaluation runner script", "name": "Evaluation runner script"}],
            "label": "Evaluation runner script",
            "name": "Input dataset",
            "outputs": [],
            "position": {"left": 10, "top": 570},
            "tool_id": None,
            "tool_state": "{\"optional\": false, \"tag\": null}",
            "tool_version": None,
            "type": "data_input",
            "uuid": str(uuid.uuid4())
        },
        "5": {
            "annotation": "Runs full cross-validated SpaGE predictions and TISSUE multi-coverage uncertainty evaluations across all measured genes.",
            "content_id": "0ce0b774-f5b9-40bc-b006-2d6cab6b939f",
            "errors": None,
            "id": 5,
            "input_connections": {
                "spatial_counts": {"id": 0, "output_name": "output"},
                "locations": {"id": 1, "output_name": "output"},
                "scrna_counts": {"id": 2, "output_name": "output"},
                "code_bundle": {"id": 3, "output_name": "output"},
                "wrapper_script": {"id": 4, "output_name": "output"}
            },
            "inputs": [],
            "label": "TISSUE Full Evaluation Tool",
            "name": "TISSUE Full Cross-Validated Evaluation",
            "outputs": [
                {"name": "per_gene_metrics", "type": "csv"},
                {"name": "overall_summary", "type": "json"},
                {"name": "calibration_curve", "type": "png"},
                {"name": "acta2_pi_map", "type": "png"},
                {"name": "all_genes_per_cell", "type": "csv"}
            ],
            "position": {"left": 350, "top": 200},
            "post_job_actions": {},
            "tool_id": "0ce0b774-f5b9-40bc-b006-2d6cab6b939f",
            "tool_state": json.dumps({
                "__page__": None,
                "__rerun_remap_job_id__": None,
                "spatial_counts": {"__class__": "ConnectedValue"},
                "locations": {"__class__": "ConnectedValue"},
                "scrna_counts": {"__class__": "ConnectedValue"},
                "code_bundle": {"__class__": "ConnectedValue"},
                "wrapper_script": {"__class__": "ConnectedValue"},
                "n_folds": 10,
                "n_pv": 10,
                "seed": 444
            }),
            "tool_version": "1.0.0",
            "tool_uuid": "0ce0b774-f5b9-40bc-b006-2d6cab6b939f",
            "type": "tool",
            "uuid": str(uuid.uuid4())
        }
    }
}

with open("workflow/tissue_full_evaluation_workflow.ga", "w") as f:
    json.dump(workflow_dict, f, indent=2)

print("Saved workflow/tissue_full_evaluation_workflow.ga")
