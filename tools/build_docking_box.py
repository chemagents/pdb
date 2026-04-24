import argparse
import json
import sys
from pathlib import Path

from rcsb_receptor_utils import normalize_gene_name
from rcsb_receptor_utils import resolve_default_output_dir
from rcsb_receptor_utils import write_json


BOX_MARGIN_ANGSTROM = 4.0


def load_json(path: str | Path) -> dict | list:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def load_primary_cluster(gene_name: str, output_dir: str | Path) -> dict:
    path = Path(output_dir) / gene_name / "04_sites" / "primary_site_cluster.json"
    if not path.exists():
        raise RuntimeError(f"Missing primary site cluster at {path}")
    data = load_json(path)
    if not isinstance(data, dict):
        raise RuntimeError(f"Primary site cluster file at {path} must contain a JSON object")
    return data


def load_site_qc(gene_name: str, output_dir: str | Path) -> dict:
    path = Path(output_dir) / gene_name / "04_sites" / "site_cluster_qc.json"
    if not path.exists():
        raise RuntimeError(f"Missing site cluster QC file at {path}")
    data = load_json(path)
    if not isinstance(data, dict):
        raise RuntimeError(f"Site cluster QC file at {path} must contain a JSON object")
    return data


def load_aligned_ligands(gene_name: str, output_dir: str | Path) -> list[dict]:
    path = Path(output_dir) / gene_name / "04_sites" / "accepted_ligand_instances_aligned.json"
    if not path.exists():
        raise RuntimeError(f"Missing aligned ligand file at {path}")
    data = load_json(path)
    if not isinstance(data, list):
        raise RuntimeError(f"Aligned ligand file at {path} must contain a JSON list")
    return data


def centroid(points: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    count = len(points)
    return (
        sum(point[0] for point in points) / count,
        sum(point[1] for point in points) / count,
        sum(point[2] for point in points) / count,
    )


def build_box(points: list[tuple[float, float, float]], margin: float) -> dict:
    min_x = min(point[0] for point in points)
    min_y = min(point[1] for point in points)
    min_z = min(point[2] for point in points)
    max_x = max(point[0] for point in points)
    max_y = max(point[1] for point in points)
    max_z = max(point[2] for point in points)

    padded_min_x = min_x - margin
    padded_min_y = min_y - margin
    padded_min_z = min_z - margin
    padded_max_x = max_x + margin
    padded_max_y = max_y + margin
    padded_max_z = max_z + margin

    center_x = (padded_min_x + padded_max_x) / 2.0
    center_y = (padded_min_y + padded_max_y) / 2.0
    center_z = (padded_min_z + padded_max_z) / 2.0

    size_x = padded_max_x - padded_min_x
    size_y = padded_max_y - padded_min_y
    size_z = padded_max_z - padded_min_z

    return {
        "raw_min": [min_x, min_y, min_z],
        "raw_max": [max_x, max_y, max_z],
        "padded_min": [padded_min_x, padded_min_y, padded_min_z],
        "padded_max": [padded_max_x, padded_max_y, padded_max_z],
        "box_center": [center_x, center_y, center_z],
        "box_size": [size_x, size_y, size_z],
    }


def build_box_inputs(primary_cluster: dict, aligned_entries: list[dict]) -> dict:
    selected_ids = set(primary_cluster["member_ligand_instance_ids"])
    selected_entries = [entry for entry in aligned_entries if entry["ligand_instance_id"] in selected_ids]
    if not selected_entries:
        raise RuntimeError("Primary site cluster did not match any aligned ligand entries")

    point_cloud = []
    centroid_points = []
    for entry in selected_entries:
        centroid_points.append(tuple(entry["aligned_centroid"]))
        for point in entry["aligned_atom_coordinates"]:
            point_cloud.append((point[0], point[1], point[2]))

    return {
        "selected_entries": selected_entries,
        "point_cloud": point_cloud,
        "centroid_points": centroid_points,
    }


def build_final_box_record(gene_name: str, site_qc: dict, primary_cluster: dict, box_geometry: dict) -> dict:
    return {
        "gene_name": gene_name,
        "reference_protomer_instance_id": site_qc["reference_protomer_instance_id"],
        "site_cluster_id": primary_cluster["site_cluster_id"],
        "supporting_structure_ids": primary_cluster["supporting_structure_ids"],
        "ligand_instances_used": primary_cluster["member_count"],
        "box_center": box_geometry["box_center"],
        "box_size": box_geometry["box_size"],
        "box_margin_angstrom": BOX_MARGIN_ANGSTROM,
        "notes": [
            f"Built from primary site cluster {primary_cluster['site_cluster_id']}",
            f"Reference protomer: {site_qc['reference_protomer_instance_id']}",
            f"Padding margin: {BOX_MARGIN_ANGSTROM} A",
        ],
    }


def build_box_qc(gene_name: str, site_qc: dict, primary_cluster: dict, box_geometry: dict, selected_entries: list[dict]) -> dict:
    size_x, size_y, size_z = box_geometry["box_size"]
    centroid_of_ligand_centroids = centroid([tuple(entry["aligned_centroid"]) for entry in selected_entries])
    center_x, center_y, center_z = box_geometry["box_center"]

    return {
        "gene_name": gene_name,
        "pipeline_stage": 5,
        "stage_name": "docking_box_construction",
        "reference_protomer_instance_id": site_qc["reference_protomer_instance_id"],
        "site_cluster_id": primary_cluster["site_cluster_id"],
        "supporting_structure_count": len(primary_cluster["supporting_structure_ids"]),
        "ligand_instance_count": len(selected_entries),
        "raw_min": box_geometry["raw_min"],
        "raw_max": box_geometry["raw_max"],
        "padded_min": box_geometry["padded_min"],
        "padded_max": box_geometry["padded_max"],
        "box_center": box_geometry["box_center"],
        "box_size": box_geometry["box_size"],
        "box_margin_angstrom": BOX_MARGIN_ANGSTROM,
        "centroid_of_ligand_centroids": [
            centroid_of_ligand_centroids[0],
            centroid_of_ligand_centroids[1],
            centroid_of_ligand_centroids[2],
        ],
        "box_center_minus_ligand_centroid": [
            center_x - centroid_of_ligand_centroids[0],
            center_y - centroid_of_ligand_centroids[1],
            center_z - centroid_of_ligand_centroids[2],
        ],
        "vina_large_box_warning": any(size > 30.0 for size in (size_x, size_y, size_z)),
    }


def build_vina_text(final_box_record: dict) -> str:
    center_x, center_y, center_z = final_box_record["box_center"]
    size_x, size_y, size_z = final_box_record["box_size"]
    lines = [
        f"center_x = {center_x:.4f}",
        f"center_y = {center_y:.4f}",
        f"center_z = {center_z:.4f}",
        f"size_x = {size_x:.4f}",
        f"size_y = {size_y:.4f}",
        f"size_z = {size_z:.4f}",
    ]
    return "\n".join(lines) + "\n"


def write_box_report(final_box_record: dict, qc_record: dict, output_path: str | Path) -> None:
    lines = [
        f"# Final Box Summary: {final_box_record['gene_name']}",
        "",
        f"- Reference protomer: `{final_box_record['reference_protomer_instance_id']}`",
        f"- Site cluster: `{final_box_record['site_cluster_id']}`",
        f"- Supporting structures: {len(final_box_record['supporting_structure_ids'])}",
        f"- Ligand instances used: {final_box_record['ligand_instances_used']}",
        f"- Box center: {final_box_record['box_center']}",
        f"- Box size: {final_box_record['box_size']}",
        f"- Margin: {final_box_record['box_margin_angstrom']} A",
        f"- Vina large box warning: {qc_record['vina_large_box_warning']}",
    ]
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Build a docking box from the primary aligned ligand site cluster."
        )
    )
    parser.add_argument("gene_name", help="Gene symbol, for example ACHE")
    parser.add_argument(
        "--output-dir",
        default=str(resolve_default_output_dir()),
        help=(
            "Base output directory. Inputs are read from output/<gene_name>/04_sites/ and "
            "box artifacts are written into output/<gene_name>/05_box/."
        ),
    )
    args = parser.parse_args()

    try:
        gene_name = normalize_gene_name(args.gene_name)
        primary_cluster = load_primary_cluster(gene_name, args.output_dir)
        site_qc = load_site_qc(gene_name, args.output_dir)
        aligned_entries = load_aligned_ligands(gene_name, args.output_dir)

        selected_ids = set(primary_cluster["member_ligand_instance_ids"])
        selected_entries = [entry for entry in aligned_entries if entry["ligand_instance_id"] in selected_ids]
        if not selected_entries:
            raise RuntimeError("Primary site cluster did not match any aligned ligand entries")

        point_cloud = []
        for entry in selected_entries:
            for point in entry["aligned_atom_coordinates"]:
                point_cloud.append((point[0], point[1], point[2]))

        box_geometry = build_box(point_cloud, BOX_MARGIN_ANGSTROM)
        final_box_record = build_final_box_record(gene_name, site_qc, primary_cluster, box_geometry)
        qc_record = build_box_qc(gene_name, site_qc, primary_cluster, box_geometry, selected_entries)
        box_inputs = {
            "reference_protomer_instance_id": site_qc["reference_protomer_instance_id"],
            "site_cluster_id": primary_cluster["site_cluster_id"],
            "member_ligand_instance_ids": primary_cluster["member_ligand_instance_ids"],
            "supporting_structure_ids": primary_cluster["supporting_structure_ids"],
            "raw_atom_point_count": len(point_cloud),
            "box_margin_angstrom": BOX_MARGIN_ANGSTROM,
        }

        stage_dir = Path(args.output_dir) / gene_name / "05_box"
        reports_dir = stage_dir / "reports"
        exports_dir = Path(args.output_dir) / gene_name / "07_exports"
        stage_dir.mkdir(parents=True, exist_ok=True)
        reports_dir.mkdir(parents=True, exist_ok=True)
        exports_dir.mkdir(parents=True, exist_ok=True)

        write_json(final_box_record, stage_dir / "final_box.json")
        write_json(qc_record, stage_dir / "final_box_qc.json")
        write_json(box_inputs, stage_dir / "box_inputs.json")
        write_box_report(final_box_record, qc_record, reports_dir / "final_box_summary.md")
        (exports_dir / "docking_box_for_vina.txt").write_text(
            build_vina_text(final_box_record),
            encoding="utf-8",
        )
        write_json(final_box_record, exports_dir / "docking_box_for_vina.json")

        print(stage_dir / "final_box.json")
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"Failed to build docking box for {args.gene_name}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
