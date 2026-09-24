import json
import uuid

workflow_dict = {
    "a_galaxy_workflow": "true",
    "format-version": "0.1",
    "name": "TISSUE Spatial Gene Expression and Prediction Intervals (Nature 2026 Acta2)",
    "annotation": "End-to-end TISSUE pipeline reproducing Nature 2026 Extended Data Fig. 1C for Acta2 spatial gene expression prediction and conformal uncertainty prediction intervals on mouse somatosensory cortex data.",
    "tags": ["spatial-transcriptomics", "TISSUE", "uncertainty-calibration", "Nature-2026"],
    "version": 1,
    "steps": {
        "0": {
            "annotation": "Spatial counts matrix (e.g. Spatial_count.txt)",
            "content_id": None,
            "errors": None,
            "id": 0,
            "input_connections": {},
            "inputs": [
                {
                    "description": "Spatial counts matrix (cells in rows, genes in columns)",
                    "name": "Spatial counts"
                }
            ],
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
            "annotation": "Spatial coordinates matrix (Locations.txt)",
            "content_id": None,
            "errors": None,
            "id": 1,
            "input_connections": {},
            "inputs": [
                {
                    "description": "Spatial cell coordinates (X, Y)",
                    "name": "Spatial locations"
                }
            ],
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
            "annotation": "Reference scRNA-seq counts matrix (scRNA_count.txt or scRNA_count.txt.gz)",
            "content_id": None,
            "errors": None,
            "id": 2,
            "input_connections": {},
            "inputs": [
                {
                    "description": "scRNA-seq count matrix (cells in columns, genes in rows)",
                    "name": "scRNA-seq counts"
                }
            ],
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
            "inputs": [
                {
                    "description": "TISSUE core library tar.gz",
                    "name": "TISSUE code package"
                }
            ],
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
            "annotation": "Python execution wrapper script (run_tissue_wrapper.py)",
            "content_id": None,
            "errors": None,
            "id": 4,
            "input_connections": {},
            "inputs": [
                {
                    "description": "Wrapper script",
                    "name": "TISSUE execution wrapper"
                }
            ],
            "label": "TISSUE execution wrapper",
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
            "annotation": "Runs SpaGE prediction, spatial graph construction, conformal uncertainty calibration, and prediction interval extraction.",
            "content_id": "tissue_spatial_prediction",
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
            "label": "TISSUE Analysis",
            "name": "TISSUE Spatial Gene Expression and Prediction Intervals",
            "outputs": [
                {"name": "pi_plot", "type": "png"},
                {"name": "prediction_plot", "type": "png"},
                {"name": "per_cell_table", "type": "csv"},
                {"name": "summary_metrics", "type": "json"}
            ],
            "position": {"left": 350, "top": 200},
            "post_job_actions": {},
            "tool_id": "tissue_spatial_prediction",
            "tool_state": json.dumps({
                "__page__": None,
                "__rerun_remap_job_id__": None,
                "spatial_counts": {"__class__": "ConnectedValue"},
                "locations": {"__class__": "ConnectedValue"},
                "scrna_counts": {"__class__": "ConnectedValue"},
                "code_bundle": {"__class__": "ConnectedValue"},
                "wrapper_script": {"__class__": "ConnectedValue"},
                "target_gene": "acta2",
                "prediction_method": "spage",
                "alpha_level": 0.23,
                "seed": 42
            }),
            "tool_version": "1.2.0",
            "type": "tool",
            "uuid": str(uuid.uuid4())
        }
    }
}

with open("workflow/tissue_spatial_prediction_workflow.ga", "w") as f:
    json.dump(workflow_dict, f, indent=2)

print("Saved workflow/tissue_spatial_prediction_workflow.ga with valid UUIDs")
