import argparse
import json
import subprocess
import sys
from pathlib import Path

from rcsb_receptor_utils import normalize_gene_name
from rcsb_receptor_utils import resolve_default_output_dir


STAGE_COMMANDS = {
    "fetch-candidates": ["fetch_receptor_by_gene_name.py"],
    "select-structures": ["select_holo_structures.py"],
    "extract-protomers": ["extract_protomers.py"],
    "align-protomers": ["align_protomers.py"],
    "cluster-sites": ["cluster_binding_sites.py"],
    "build-box": ["build_docking_box.py"],
    "visualize": ["generate_visual_artifacts.py"],
}

RUN_ALL_SEQUENCE = [
    "fetch-candidates",
    "select-structures",
    "extract-protomers",
    "align-protomers",
    "cluster-sites",
    "build-box",
    "visualize",
]


def run_stage(stage_name: str, gene_name: str, output_dir: str | Path) -> int:
    script_name = STAGE_COMMANDS[stage_name][0]
    command = [sys.executable, str(Path(__file__).resolve().parent / script_name), gene_name, "--output-dir", str(output_dir)]
    completed = subprocess.run(command, check=False)
    return completed.returncode


def load_json_if_exists(path: Path) -> dict | list | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def build_site_pipeline_summary(gene_name: str, output_dir: str | Path) -> dict:
    gene_root = Path(output_dir) / gene_name

    selection_index = load_json_if_exists(gene_root / "01_selection" / "index_selected_structures.json")
    protomer_index = load_json_if_exists(gene_root / "02_protomers" / "index_protomers.json")
    alignment_qc = load_json_if_exists(gene_root / "03_alignment" / "alignment_qc.json")
    site_qc = load_json_if_exists(gene_root / "04_sites" / "site_cluster_qc.json")
    final_box = load_json_if_exists(gene_root / "05_box" / "final_box.json")

    summary = {
        "gene_name": gene_name,
        "pipeline_status": "completed",
        "artifacts": {
            "selection_index": str(gene_root / "01_selection" / "index_selected_structures.json"),
            "protomer_index": str(gene_root / "02_protomers" / "index_protomers.json"),
            "alignment_qc": str(gene_root / "03_alignment" / "alignment_qc.json"),
            "site_cluster_qc": str(gene_root / "04_sites" / "site_cluster_qc.json"),
            "final_box": str(gene_root / "05_box" / "final_box.json"),
            "visual_review": str(gene_root / "06_visuals" / "overlay_review.html"),
            "vina_box": str(gene_root / "07_exports" / "docking_box_for_vina.txt"),
            "visualize_ligand_cloud_py3dmol": str(gene_root / "07_exports" / "visualize_ligand_cloud_py3dmol.py"),
            "visualize_box_comparison_py3dmol": str(gene_root / "07_exports" / "visualize_box_comparison_py3dmol.py"),
        },
    }

    if isinstance(selection_index, dict):
        summary["selection"] = {
            "input_structure_count": selection_index.get("input_structure_count"),
            "accepted_count": selection_index.get("accepted_count"),
            "uncertain_count": selection_index.get("uncertain_count"),
            "rejected_count": selection_index.get("rejected_count"),
            "accepted_structure_ids": selection_index.get("accepted_structure_ids"),
        }

    if isinstance(protomer_index, dict):
        summary["protomers"] = {
            "accepted_structure_count": protomer_index.get("accepted_structure_count"),
            "accepted_protomer_count": protomer_index.get("accepted_protomer_count"),
            "accepted_ligand_instance_count": protomer_index.get("accepted_ligand_instance_count"),
            "canonical_reference_candidate_id": protomer_index.get("canonical_reference_candidate_id"),
        }

    if isinstance(alignment_qc, dict):
        summary["alignment"] = {
            "reference_protomer_instance_id": alignment_qc.get("reference_protomer_instance_id"),
            "accepted_alignment_count": alignment_qc.get("accepted_alignment_count"),
            "requires_review_alignment_count": alignment_qc.get("requires_review_alignment_count"),
        }

    if isinstance(site_qc, dict):
        summary["sites"] = {
            "aligned_ligand_instance_count": site_qc.get("aligned_ligand_instance_count"),
            "cluster_count": site_qc.get("cluster_count"),
            "primary_site_cluster_id": site_qc.get("primary_site_cluster_id"),
            "primary_site_supporting_structure_ids": site_qc.get("primary_site_supporting_structure_ids"),
        }

    if isinstance(final_box, dict):
        summary["box"] = {
            "reference_protomer_instance_id": final_box.get("reference_protomer_instance_id"),
            "site_cluster_id": final_box.get("site_cluster_id"),
            "ligand_instances_used": final_box.get("ligand_instances_used"),
            "box_center": final_box.get("box_center"),
            "box_size": final_box.get("box_size"),
        }

    return summary


def write_site_pipeline_summary(gene_name: str, output_dir: str | Path) -> Path:
    gene_root = Path(output_dir) / gene_name
    exports_dir = gene_root / "07_exports"
    exports_dir.mkdir(parents=True, exist_ok=True)
    summary_path = exports_dir / "site_pipeline_summary.json"
    summary_path.write_text(
        json.dumps(build_site_pipeline_summary(gene_name, output_dir), ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary_path


def print_run_all_outputs(gene_name: str, output_dir: str | Path) -> None:
    gene_root = Path(output_dir) / gene_name
    key_outputs = [
        gene_root / "01_selection" / "index_selected_structures.json",
        gene_root / "02_protomers" / "index_protomers.json",
        gene_root / "03_alignment" / "alignment_qc.json",
        gene_root / "04_sites" / "site_cluster_qc.json",
        gene_root / "05_box" / "final_box.json",
        gene_root / "05_box" / "final_box_qc.json",
        gene_root / "06_visuals" / "overlay_review.html",
        gene_root / "07_exports" / "docking_box_for_vina.txt",
        gene_root / "07_exports" / "site_pipeline_summary.json",
        gene_root / "07_exports" / "visualize_ligand_cloud_py3dmol.py",
        gene_root / "07_exports" / "visualize_box_comparison_py3dmol.py",
    ]

    optional_outputs = [
        gene_root / "06_visuals_ref" / "overlay_review_ref.html",
    ]

    print("\nRun-all completed. Key outputs:")
    for path in key_outputs:
        if path.exists():
            print(path)

    for path in optional_outputs:
        if path.exists():
            print(path)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run the staged receptor-site pipeline step by step or end to end."
    )
    parser.add_argument(
        "command",
        choices=[*STAGE_COMMANDS.keys(), "run-all"],
        help="Pipeline stage to run, or run-all for the full sequence.",
    )
    parser.add_argument("gene_name", help="Gene symbol, for example ACHE or CNR1")
    parser.add_argument(
        "--output-dir",
        default=str(resolve_default_output_dir()),
        help="Base output directory for all stage artifacts.",
    )
    args = parser.parse_args()
    gene_name = normalize_gene_name(args.gene_name)

    if args.command == "run-all":
        for stage_name in RUN_ALL_SEQUENCE:
            return_code = run_stage(stage_name, gene_name, args.output_dir)
            if return_code != 0:
                return return_code
        write_site_pipeline_summary(gene_name, args.output_dir)
        print_run_all_outputs(gene_name, args.output_dir)
        return 0

    return run_stage(args.command, gene_name, args.output_dir)


if __name__ == "__main__":
    raise SystemExit(main())
