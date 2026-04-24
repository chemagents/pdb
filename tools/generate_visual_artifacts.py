import argparse
import subprocess
import sys
from pathlib import Path

from rcsb_receptor_utils import normalize_gene_name
from rcsb_receptor_utils import resolve_default_output_dir


def run_script(script_name: str, gene_name: str, output_dir: str | Path) -> int:
    script_path = Path(__file__).resolve().parent / script_name
    command = [sys.executable, str(script_path), gene_name, "--output-dir", str(output_dir)]
    completed = subprocess.run(command, check=False)
    return completed.returncode


def referent_exists_for_visual_ref(gene_name: str, output_dir: str | Path) -> bool:
    output_root = Path(output_dir) / gene_name
    reference_path = output_root / "03_alignment" / "reference_protomer.json"
    if not reference_path.exists():
        return False

    import json

    reference = json.loads(reference_path.read_text(encoding="utf-8"))
    structure_id = str(reference.get("structure_id", "")).lower()
    if not structure_id:
        return False
    referent_path = Path(__file__).resolve().parent.parent / "referents" / f"{structure_id}_config.txt"
    return referent_path.exists()


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Generate mandatory visualization artifacts: 2D overlays, optional 2D box comparison, and 3D py3Dmol script."
        )
    )
    parser.add_argument("gene_name", help="Gene symbol, for example ACHE or CNR1")
    parser.add_argument(
        "--output-dir",
        default=str(resolve_default_output_dir()),
        help="Base output directory for pipeline artifacts.",
    )
    args = parser.parse_args()

    gene_name = normalize_gene_name(args.gene_name)

    for script_name in [
        "render_site_overlay.py",
        "visualize_ligand_cloud_py3dmol.py",
        "visualize_box_comparison_py3dmol.py",
    ]:
        return_code = run_script(script_name, gene_name, args.output_dir)
        if return_code != 0:
            return return_code

    if referent_exists_for_visual_ref(gene_name, args.output_dir):
        return_code = run_script("render_site_overlay_ref.py", gene_name, args.output_dir)
        if return_code != 0:
            return return_code

    print(Path(args.output_dir) / gene_name / "06_visuals" / "overlay_review.html")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
