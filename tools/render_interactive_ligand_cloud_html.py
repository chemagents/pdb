import argparse
import json
import shlex
from pathlib import Path

from rcsb_receptor_utils import normalize_gene_name
from rcsb_receptor_utils import resolve_default_output_dir


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


def load_stage_artifacts(gene_name: str, output_dir: str | Path) -> tuple[dict, list[dict], dict, dict]:
    output_root = Path(output_dir) / gene_name
    reference_protomer = load_json(output_root / "03_alignment" / "reference_protomer.json")
    aligned_ligands = load_json(output_root / "04_sites" / "accepted_ligand_instances_aligned.json")
    primary_cluster = load_json(output_root / "04_sites" / "primary_site_cluster.json")
    final_box = load_json(output_root / "05_box" / "final_box_qc.json")
    if not isinstance(reference_protomer, dict) or not isinstance(primary_cluster, dict) or not isinstance(final_box, dict):
        raise RuntimeError("Expected JSON objects for reference protomer, primary cluster and final box")
    if not isinstance(aligned_ligands, list):
        raise RuntimeError("Expected aligned ligand JSON list")
    return reference_protomer, aligned_ligands, primary_cluster, final_box


def select_ligands(aligned_ligands: list[dict], primary_cluster: dict, include_alternative: bool) -> list[dict]:
    if include_alternative:
        return aligned_ligands
    selected_ids = set(primary_cluster["member_ligand_instance_ids"])
    return [entry for entry in aligned_ligands if entry["ligand_instance_id"] in selected_ids]


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


def build_ca_trace_pdb(reference_protomer: dict, output_root: Path) -> str:
    mmcif_path = output_root / reference_protomer["source_mmcif_file"]
    chain_id = reference_protomer["auth_asym_id"]
    lines = []
    serial = 1
    for row, field_index in iter_atom_site_rows(mmcif_path):
        if row[field_index["_atom_site.auth_asym_id"]] != chain_id:
            continue
        if row[field_index["_atom_site.label_atom_id"]] != "CA":
            continue
        auth_comp_id = row[field_index["_atom_site.auth_comp_id"]]
        auth_seq_id = row[field_index["_atom_site.auth_seq_id"]]
        try:
            x = float(row[field_index["_atom_site.Cartn_x"]])
            y = float(row[field_index["_atom_site.Cartn_y"]])
            z = float(row[field_index["_atom_site.Cartn_z"]])
            resseq = int(auth_seq_id)
        except ValueError:
            continue
        lines.append(
            f"ATOM  {serial:5d}  CA  {auth_comp_id[:3].rjust(3)} A{resseq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C"
        )
        serial += 1
    lines.append("END")
    return "\n".join(lines) + "\n"


def build_html_document(
    gene_name: str,
    reference_protomer: dict,
    ca_trace_pdb: str,
    ligand_entries: list[dict],
    final_box: dict,
    include_alternative: bool,
) -> str:
    box_center = final_box["box_center"]
    box_size = final_box["box_size"]

    ligand_payload = []
    palette = [
        "#16a34a", "#0891b2", "#c026d3", "#ca8a04",
        "#ea580c", "#7c3aed", "#2563eb", "#dc2626",
    ]
    for index, entry in enumerate(ligand_entries):
        ligand_payload.append(
            {
                "ligand_instance_id": entry["ligand_instance_id"],
                "ligand_code": entry["ligand_code"],
                "aligned_centroid": entry["aligned_centroid"],
                "aligned_atom_coordinates": entry["aligned_atom_coordinates"],
                "cluster_role": entry.get("cluster_role"),
                "site_cluster_id": entry.get("site_cluster_id"),
                "color": palette[index % len(palette)],
            }
        )

    return f"""<!DOCTYPE html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\" />
  <title>Ligand Cloud Viewer: {gene_name}</title>
  <script src=\"https://3Dmol.org/build/3Dmol-min.js\"></script>
  <style>
    body {{ margin: 0; font-family: Arial, sans-serif; background: #f8fafc; color: #111827; }}
    #app {{ display: flex; height: 100vh; }}
    #viewer {{ flex: 1; min-width: 0; }}
    #sidebar {{ width: 360px; overflow: auto; background: white; border-left: 1px solid #d1d5db; padding: 16px; box-sizing: border-box; }}
    h1 {{ font-size: 20px; margin: 0 0 12px 0; }}
    h2 {{ font-size: 16px; margin: 20px 0 8px 0; }}
    .meta {{ font-size: 13px; line-height: 1.45; }}
    .control {{ margin: 10px 0; }}
    .legend-item {{ display: flex; align-items: center; gap: 8px; font-size: 12px; margin: 4px 0; }}
    .swatch {{ width: 12px; height: 12px; border-radius: 999px; flex: 0 0 12px; }}
    code {{ background: #f3f4f6; padding: 2px 4px; }}
  </style>
</head>
<body>
  <div id=\"app\">
    <div id=\"viewer\"></div>
    <div id=\"sidebar\">
      <h1>Ligand Cloud Viewer: {gene_name}</h1>
      <div class=\"meta\">
        <div>Reference protomer: <code>{reference_protomer['protomer_instance_id']}</code></div>
        <div>Include alternative clusters: <code>{include_alternative}</code></div>
        <div>Ligands rendered: <code>{len(ligand_payload)}</code></div>
        <div>Box center: <code>{box_center}</code></div>
        <div>Box size: <code>{box_size}</code></div>
      </div>

      <h2>Controls</h2>
      <div class=\"control\"><label><input type=\"checkbox\" id=\"toggleProtein\" checked /> Protein backbone</label></div>
      <div class=\"control\"><label><input type=\"checkbox\" id=\"toggleAtomBouquet\" checked /> Ligand atom bouquet</label></div>
      <div class=\"control\"><label><input type=\"checkbox\" id=\"toggleCentroids\" checked /> Ligand centroid spheres</label></div>
      <div class=\"control\"><label><input type=\"checkbox\" id=\"toggleLabels\" checked /> Ligand labels</label></div>
      <div class=\"control\"><label><input type=\"checkbox\" id=\"toggleBox\" checked /> Docking box</label></div>

      <h2>Legend</h2>
      <div id=\"legend\"></div>
    </div>
  </div>

  <script>
    const proteinPdb = {json.dumps(ca_trace_pdb)};
    const ligandPayload = {json.dumps(ligand_payload)};
    const boxCenter = {json.dumps(box_center)};
    const boxSize = {json.dumps(box_size)};

    const viewer = $3Dmol.createViewer('viewer', {{ backgroundColor: 'white' }});
    viewer.addModel(proteinPdb, 'pdb');
    viewer.setStyle({{model: 0}}, {{line: {{color: '#9ca3af', linewidth: 1.2}}}});

    let proteinVisible = true;
    let bouquetVisible = true;
    let centroidsVisible = true;
    let labelsVisible = true;
    let boxVisible = true;

    const centroidIds = [];
    const atomSphereIds = [];
    const labelIds = [];
    let boxShape = null;

    function addLigands() {{
      ligandPayload.forEach((ligand) => {{
        ligand.aligned_atom_coordinates.forEach((coord) => {{
          const shape = viewer.addSphere({{
            center: {{x: coord[0], y: coord[1], z: coord[2]}},
            radius: 0.55,
            color: ligand.color,
            opacity: 0.82
          }});
          atomSphereIds.push(shape);
        }});

        const centroid = ligand.aligned_centroid;
        const centroidShape = viewer.addSphere({{
          center: {{x: centroid[0], y: centroid[1], z: centroid[2]}},
          radius: 1.45,
          color: ligand.color,
          opacity: 0.96
        }});
        centroidIds.push(centroidShape);

        const label = viewer.addLabel(ligand.ligand_code, {{
          position: {{x: centroid[0], y: centroid[1], z: centroid[2]}},
          backgroundColor: 'white',
          fontColor: 'black',
          fontSize: 12,
          showBackground: true,
          inFront: true
        }});
        labelIds.push(label);
      }});
    }}

    function addBox() {{
      boxShape = viewer.addBox({{
        center: {{x: boxCenter[0], y: boxCenter[1], z: boxCenter[2]}},
        dimensions: {{w: boxSize[0], h: boxSize[1], d: boxSize[2]}},
        color: 'cyan',
        opacity: 0.22
      }});
    }}

    function updateVisibility() {{
      viewer.setStyle({{model: 0}}, proteinVisible ? {{line: {{color: '#9ca3af', linewidth: 1.2}}}} : {{}});
      atomSphereIds.forEach((shape) => shape.updateStyle({{hidden: !bouquetVisible}}));
      centroidIds.forEach((shape) => shape.updateStyle({{hidden: !centroidsVisible}}));
      labelIds.forEach((label) => label.setHidden(!labelsVisible));
      if (boxShape) boxShape.updateStyle({{hidden: !boxVisible}});
      viewer.render();
    }}

    function buildLegend() {{
      const legend = document.getElementById('legend');
      legend.innerHTML = '';
      ligandPayload.forEach((ligand) => {{
        const item = document.createElement('div');
        item.className = 'legend-item';
        const swatch = document.createElement('div');
        swatch.className = 'swatch';
        swatch.style.background = ligand.color;
        const text = document.createElement('div');
        text.textContent = `${{ligand.ligand_code}} (${{ligand.ligand_instance_id.split(':')[1]}})`;
        item.appendChild(swatch);
        item.appendChild(text);
        legend.appendChild(item);
      }});
    }}

    addLigands();
    addBox();
    buildLegend();

    viewer.zoomTo();
    viewer.zoom(1.35);
    viewer.render();

    document.getElementById('toggleProtein').addEventListener('change', (e) => {{ proteinVisible = e.target.checked; updateVisibility(); }});
    document.getElementById('toggleAtomBouquet').addEventListener('change', (e) => {{ bouquetVisible = e.target.checked; updateVisibility(); }});
    document.getElementById('toggleCentroids').addEventListener('change', (e) => {{ centroidsVisible = e.target.checked; updateVisibility(); }});
    document.getElementById('toggleLabels').addEventListener('change', (e) => {{ labelsVisible = e.target.checked; updateVisibility(); }});
    document.getElementById('toggleBox').addEventListener('change', (e) => {{ boxVisible = e.target.checked; updateVisibility(); }});
  </script>
</body>
</html>
"""


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Render standalone interactive HTML for ligand bouquet review."
    )
    parser.add_argument("gene_name", help="Gene symbol, for example ACHE")
    parser.add_argument(
        "--output-dir",
        default=str(Path(resolve_default_output_dir())),
        help="Base output directory containing pipeline artifacts.",
    )
    parser.add_argument(
        "--include-alternative",
        action="store_true",
        help="Include ligands from alternative site clusters as well.",
    )
    parser.add_argument(
        "--write-path",
        default=None,
        help="Optional output path. Defaults to output/<gene>/07_exports/ligand_cloud_viewer.html.",
    )
    args = parser.parse_args()

    gene_name = normalize_gene_name(args.gene_name)
    output_root = Path(args.output_dir)
    reference_protomer, aligned_ligands, primary_cluster, final_box = load_stage_artifacts(gene_name, output_root)
    selected_ligands = select_ligands(aligned_ligands, primary_cluster, args.include_alternative)
    ca_trace_pdb = build_ca_trace_pdb(reference_protomer, output_root)

    write_path = Path(args.write_path) if args.write_path else output_root / gene_name / "07_exports" / "ligand_cloud_viewer.html"
    write_path.parent.mkdir(parents=True, exist_ok=True)
    write_path.write_text(
        build_html_document(
            gene_name,
            reference_protomer,
            ca_trace_pdb,
            selected_ligands,
            final_box,
            args.include_alternative,
        ),
        encoding="utf-8",
    )
    print(write_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
