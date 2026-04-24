import argparse
import json
import shlex
import sys
from pathlib import Path

from rcsb_receptor_utils import normalize_gene_name
from rcsb_receptor_utils import resolve_default_output_dir
from rcsb_receptor_utils import write_json


SVG_WIDTH = 900
SVG_HEIGHT = 300
PANEL_WIDTH = 260
PANEL_HEIGHT = 240
PANEL_GAP = 20
MARGIN = 20


def load_json(path: str | Path) -> dict | list:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def parse_loop_row(line: str) -> list[str]:
    lexer = shlex.shlex(line, posix=True)
    lexer.whitespace_split = True
    lexer.commenters = ""
    return list(lexer)


def iter_atom_site_rows(mmcif_path: str | Path):
    path = Path(mmcif_path)
    lines = path.read_text(encoding="utf-8").splitlines()

    atom_site_fields = []
    field_index = {}
    in_atom_site_loop = False

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue

        if not in_atom_site_loop:
            if line == "loop_":
                candidate_fields = []
                j = i + 1
                while j < len(lines):
                    candidate_line = lines[j].strip()
                    if candidate_line.startswith("_atom_site."):
                        candidate_fields.append(candidate_line)
                        j += 1
                        continue
                    break
                if candidate_fields:
                    atom_site_fields = candidate_fields
                    field_index = {name: idx for idx, name in enumerate(atom_site_fields)}
                    in_atom_site_loop = True
                    i = j
                    continue
            i += 1
            continue

        if line == "#" or line == "loop_" or line.startswith("_"):
            break

        row = parse_loop_row(line)
        if len(row) == len(atom_site_fields):
            yield row, field_index
        i += 1


def extract_reference_ca_points(mmcif_path: str | Path, auth_asym_id: str) -> list[tuple[float, float, float]]:
    points = []
    for row, field_index in iter_atom_site_rows(mmcif_path):
        if row[field_index["_atom_site.auth_asym_id"]] != auth_asym_id:
            continue
        if row[field_index["_atom_site.label_atom_id"]] != "CA":
            continue
        try:
            x = float(row[field_index["_atom_site.Cartn_x"]])
            y = float(row[field_index["_atom_site.Cartn_y"]])
            z = float(row[field_index["_atom_site.Cartn_z"]])
        except ValueError:
            continue
        points.append((x, y, z))
    if not points:
        raise RuntimeError(f"No CA points found in {mmcif_path} for chain {auth_asym_id}")
    return points


def project_point(point: tuple[float, float, float], plane: str) -> tuple[float, float]:
    x, y, z = point
    if plane == "xy":
        return (x, y)
    if plane == "xz":
        return (x, z)
    if plane == "yz":
        return (y, z)
    raise ValueError(f"Unsupported plane {plane}")


def scale_projected_points(
    points: list[tuple[float, float]],
    width: int,
    height: int,
    margin: int,
) -> tuple[list[tuple[float, float]], tuple[float, float, float, float, float]]:
    min_x = min(point[0] for point in points)
    max_x = max(point[0] for point in points)
    min_y = min(point[1] for point in points)
    max_y = max(point[1] for point in points)

    span_x = max(max_x - min_x, 1.0)
    span_y = max(max_y - min_y, 1.0)
    usable_width = width - 2 * margin
    usable_height = height - 2 * margin
    scale = min(usable_width / span_x, usable_height / span_y)

    scaled = []
    for x, y in points:
        sx = margin + (x - min_x) * scale
        sy = height - (margin + (y - min_y) * scale)
        scaled.append((sx, sy))
    return scaled, (min_x, max_x, min_y, max_y, scale)


def scale_with_bounds(
    points: list[tuple[float, float]],
    bounds: tuple[float, float, float, float, float],
    height: int,
    margin: int,
) -> list[tuple[float, float]]:
    min_x, _max_x, min_y, _max_y, scale = bounds
    scaled = []
    for x, y in points:
        sx = margin + (x - min_x) * scale
        sy = height - (margin + (y - min_y) * scale)
        scaled.append((sx, sy))
    return scaled


def box_edges(box: dict, plane: str) -> list[tuple[float, float]]:
    min_x, min_y, min_z = box["padded_min"]
    max_x, max_y, max_z = box["padded_max"]
    if plane == "xy":
        return [(min_x, min_y), (max_x, min_y), (max_x, max_y), (min_x, max_y), (min_x, min_y)]
    if plane == "xz":
        return [(min_x, min_z), (max_x, min_z), (max_x, max_z), (min_x, max_z), (min_x, min_z)]
    if plane == "yz":
        return [(min_y, min_z), (max_y, min_z), (max_y, max_z), (min_y, max_z), (min_y, min_z)]
    raise ValueError(f"Unsupported plane {plane}")


def center_size_to_box(center: list[float], size: list[float]) -> dict:
    cx, cy, cz = center
    sx, sy, sz = size
    return {
        "padded_min": [cx - sx / 2.0, cy - sy / 2.0, cz - sz / 2.0],
        "padded_max": [cx + sx / 2.0, cy + sy / 2.0, cz + sz / 2.0],
    }


def parse_ref_box(config_path: str | Path) -> dict:
    values = {}
    for raw_line in Path(config_path).read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = [part.strip() for part in line.split("=", 1)]
        values[key] = value

    required = ["center_x", "center_y", "center_z", "size_x", "size_y", "size_z"]
    missing = [key for key in required if key not in values]
    if missing:
        raise RuntimeError(f"Reference config {config_path} is missing keys: {missing}")

    center = [float(values["center_x"]), float(values["center_y"]), float(values["center_z"])]
    size = [float(values["size_x"]), float(values["size_y"]), float(values["size_z"])]
    return {
        "config_path": str(config_path),
        "box_center": center,
        "box_size": size,
        **center_size_to_box(center, size),
    }


def build_panel_svg(
    plane: str,
    reference_points: list[tuple[float, float, float]],
    ligand_entries: list[dict],
    pipeline_box: dict,
    ref_box: dict,
    title: str,
) -> str:
    projected_reference = [project_point(point, plane) for point in reference_points]
    projected_ligands = []
    for entry in ligand_entries:
        for point in entry["aligned_atom_coordinates"]:
            projected_ligands.append(project_point((point[0], point[1], point[2]), plane))

    combined = projected_reference + projected_ligands + box_edges(pipeline_box, plane) + box_edges(ref_box, plane)
    scaled_reference, bounds = scale_projected_points(combined, PANEL_WIDTH, PANEL_HEIGHT, MARGIN)
    reference_count = len(projected_reference)
    ligand_start = reference_count
    ligand_end = reference_count + len(projected_ligands)
    scaled_ligands = scaled_reference[ligand_start:ligand_end]
    scaled_pipeline_box = scale_with_bounds(box_edges(pipeline_box, plane), bounds, PANEL_HEIGHT, MARGIN)
    scaled_ref_box = scale_with_bounds(box_edges(ref_box, plane), bounds, PANEL_HEIGHT, MARGIN)
    scaled_reference_only = scaled_reference[:reference_count]

    path_d = " ".join(
        [f"M {scaled_reference_only[0][0]:.2f} {scaled_reference_only[0][1]:.2f}"]
        + [f"L {point[0]:.2f} {point[1]:.2f}" for point in scaled_reference_only[1:]]
    )
    pipeline_box_d = " ".join(
        [f"M {scaled_pipeline_box[0][0]:.2f} {scaled_pipeline_box[0][1]:.2f}"]
        + [f"L {point[0]:.2f} {point[1]:.2f}" for point in scaled_pipeline_box[1:]]
    )
    ref_box_d = " ".join(
        [f"M {scaled_ref_box[0][0]:.2f} {scaled_ref_box[0][1]:.2f}"]
        + [f"L {point[0]:.2f} {point[1]:.2f}" for point in scaled_ref_box[1:]]
    )

    circles = []
    for point in scaled_ligands:
        circles.append(
            f'<circle cx="{point[0]:.2f}" cy="{point[1]:.2f}" r="2.4" fill="#b91c1c" fill-opacity="0.65" />'
        )

    return (
        f'<g>'
        f'<rect x="0" y="0" width="{PANEL_WIDTH}" height="{PANEL_HEIGHT}" fill="#ffffff" stroke="#d1d5db" />'
        f'<text x="12" y="18" font-size="13" font-family="Arial" fill="#111827">{title}</text>'
        f'<path d="{path_d}" stroke="#9ca3af" stroke-width="1.2" fill="none" />'
        f'<path d="{pipeline_box_d}" stroke="#2563eb" stroke-width="1.6" fill="none" stroke-dasharray="6 4" />'
        f'<path d="{ref_box_d}" stroke="#7c3aed" stroke-width="1.6" fill="none" stroke-dasharray="3 5" />'
        + "".join(circles)
        + '</g>'
    )


def build_multi_panel_svg(reference_points: list[tuple[float, float, float]], ligand_entries: list[dict], pipeline_box: dict, ref_box: dict, title_prefix: str) -> str:
    panels = []
    for index, plane in enumerate(["xy", "xz", "yz"]):
        panel_svg = build_panel_svg(plane, reference_points, ligand_entries, pipeline_box, ref_box, f"{title_prefix} ({plane.upper()})")
        x_offset = index * (PANEL_WIDTH + PANEL_GAP)
        panels.append(f'<g transform="translate({x_offset},0)">{panel_svg}</g>')
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{SVG_WIDTH}" height="{SVG_HEIGHT}" viewBox="0 0 {SVG_WIDTH} {SVG_HEIGHT}">'
        f'<rect x="0" y="0" width="{SVG_WIDTH}" height="{SVG_HEIGHT}" fill="#f9fafb" />'
        + "".join(panels)
        + '</svg>'
    )


def write_text(path: str | Path, text: str) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(text, encoding="utf-8")


def build_html_review(reference_protomer: dict, final_box: dict, ref_box: dict, primary_entries: list[dict], outlier_entries: list[dict]) -> str:
    return f"""<!DOCTYPE html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\" />
  <title>Site Overlay Review (ref): {reference_protomer['gene_name']}</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 24px; color: #111827; background: #f9fafb; }}
    h1, h2 {{ margin-bottom: 8px; }}
    .meta {{ margin-bottom: 18px; }}
    .card {{ background: white; border: 1px solid #d1d5db; padding: 16px; margin-bottom: 20px; }}
    code {{ background: #f3f4f6; padding: 2px 4px; }}
  </style>
</head>
<body>
  <h1>Site Overlay Review (ref): {reference_protomer['gene_name']}</h1>
  <div class=\"meta\">
    <div>Reference protomer: <code>{reference_protomer['protomer_instance_id']}</code></div>
    <div>Primary site cluster: <code>{final_box['site_cluster_id']}</code></div>
    <div>Ligand instances used: {final_box['ligand_instances_used']}</div>
    <div>Outliers: {len(outlier_entries)}</div>
    <div>Blue dashed box: pipeline box</div>
    <div>Purple dashed box: human reference box from <code>{Path(ref_box['config_path']).name}</code></div>
  </div>
  <div class=\"card\">
    <h2>Accepted Overlay</h2>
    <p>Reference backbone in gray, aligned ligand atoms in red, pipeline box in blue dashed outline, human box in purple dashed outline.</p>
    <img src=\"overlay_all_accepted_ref.svg\" alt=\"Accepted overlay ref\" />
  </div>
  <div class=\"card\">
    <h2>Outlier Overlay</h2>
    <p>{'No outlier ligand poses were identified for this run.' if not outlier_entries else 'Outlier poses are shown in the same frame for manual review.'}</p>
    <img src=\"overlay_outliers_ref.svg\" alt=\"Outlier overlay ref\" />
  </div>
  <div class=\"card\">
    <h2>Final Box Overlay</h2>
    <p>Primary site cluster with both pipeline and human reference boxes shown together.</p>
    <img src=\"final_box_overlay_ref.svg\" alt=\"Final box overlay ref\" />
  </div>
</body>
</html>
"""


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Render SVG/HTML visualization artifacts with an overlaid human reference box."
    )
    parser.add_argument("gene_name", help="Gene symbol, for example ACHE")
    parser.add_argument(
        "--output-dir",
        default=str(resolve_default_output_dir()),
        help=(
            "Base output directory. Inputs are read from output/<gene_name>/03_alignment/, "
            "04_sites/, 05_box/ and referents/. Outputs are written into output/<gene_name>/06_visuals_ref/."
        ),
    )
    parser.add_argument(
        "--referents-dir",
        default=r"C:\CODE\pdb\referents",
        help="Directory containing human-annotated *_config.txt box files.",
    )
    args = parser.parse_args()

    try:
        gene_name = normalize_gene_name(args.gene_name)
        reference_protomer = load_json(Path(args.output_dir) / gene_name / "03_alignment" / "reference_protomer.json")
        aligned_entries = load_json(Path(args.output_dir) / gene_name / "04_sites" / "accepted_ligand_instances_aligned.json")
        outlier_entries = load_json(Path(args.output_dir) / gene_name / "04_sites" / "rejected_ligand_instances.json")
        final_box = load_json(Path(args.output_dir) / gene_name / "05_box" / "final_box.json")
        final_box_qc = load_json(Path(args.output_dir) / gene_name / "05_box" / "final_box_qc.json")

        if not isinstance(reference_protomer, dict) or not isinstance(final_box, dict) or not isinstance(final_box_qc, dict):
            raise RuntimeError("Reference protomer and final box inputs must be JSON objects")
        if not isinstance(aligned_entries, list) or not isinstance(outlier_entries, list):
            raise RuntimeError("Aligned and outlier ligand inputs must be JSON lists")

        reference_mmcif = Path(args.output_dir) / reference_protomer["source_mmcif_file"]
        reference_points = extract_reference_ca_points(reference_mmcif, reference_protomer["auth_asym_id"])

        primary_entries = [entry for entry in aligned_entries if entry.get("cluster_role") == "primary"]
        if not primary_entries:
            raise RuntimeError(f"No primary aligned ligand entries found for {gene_name}")

        reference_structure_id = reference_protomer["structure_id"].lower()
        ref_config_path = Path(args.referents_dir) / f"{reference_structure_id}_config.txt"
        if not ref_config_path.exists():
            raise RuntimeError(f"Missing human reference config for {reference_structure_id} at {ref_config_path}")
        ref_box = parse_ref_box(ref_config_path)

        visuals_dir = Path(args.output_dir) / gene_name / "06_visuals_ref"
        manifests_dir = visuals_dir / "manifests"
        visuals_dir.mkdir(parents=True, exist_ok=True)
        manifests_dir.mkdir(parents=True, exist_ok=True)

        accepted_svg = build_multi_panel_svg(reference_points, primary_entries, final_box_qc, ref_box, "Accepted overlay")
        outlier_svg = build_multi_panel_svg(reference_points, outlier_entries, final_box_qc, ref_box, "Outlier overlay")
        box_svg = build_multi_panel_svg(reference_points, primary_entries, final_box_qc, ref_box, "Final box overlay")

        write_text(visuals_dir / "overlay_all_accepted_ref.svg", accepted_svg)
        write_text(visuals_dir / "overlay_outliers_ref.svg", outlier_svg)
        write_text(visuals_dir / "final_box_overlay_ref.svg", box_svg)
        write_text(visuals_dir / "overlay_review_ref.html", build_html_review(reference_protomer, final_box, ref_box, primary_entries, outlier_entries))

        visuals_manifest = {
            "gene_name": gene_name,
            "pipeline_stage": "ref_review",
            "reference_protomer_instance_id": reference_protomer["protomer_instance_id"],
            "reference_config_file": str(ref_config_path),
            "accepted_overlay_file": f"output/{gene_name}/06_visuals_ref/overlay_all_accepted_ref.svg",
            "outlier_overlay_file": f"output/{gene_name}/06_visuals_ref/overlay_outliers_ref.svg",
            "final_box_overlay_file": f"output/{gene_name}/06_visuals_ref/final_box_overlay_ref.svg",
            "html_review_file": f"output/{gene_name}/06_visuals_ref/overlay_review_ref.html",
        }
        write_json(visuals_manifest, manifests_dir / "visuals_ref.json")

        print(visuals_dir / "overlay_review_ref.html")
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"Failed to render reference overlay for {args.gene_name}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
