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


def find_ligand_atom_metadata(entry: dict, output_root: Path) -> list[dict]:
    mmcif_path = output_root / entry["source_mmcif_file"]
    ligand_code = entry["ligand_code"]
    locator = entry["source_ligand_locator"]
    auth_asym_id = locator["auth_asym_id"]
    auth_seq_id = locator["auth_seq_id"]

    atoms = []
    for row, field_index in iter_atom_site_rows(mmcif_path):
        if row[field_index["_atom_site.auth_asym_id"]] != auth_asym_id:
            continue
        if row[field_index["_atom_site.auth_seq_id"]] != auth_seq_id:
            continue
        auth_comp_id = row[field_index["_atom_site.auth_comp_id"]]
        label_comp_id = row[field_index["_atom_site.label_comp_id"]]
        if ligand_code not in {auth_comp_id, label_comp_id}:
            continue
        atom_name = row[field_index["_atom_site.label_atom_id"]]
        element = row[field_index["_atom_site.type_symbol"]]
        atoms.append(
            {
                "atom_name": atom_name,
                "element": element,
            }
        )

    if not atoms:
        raise RuntimeError(
            f"Could not recover ligand atom metadata for {entry['ligand_instance_id']} from {mmcif_path}"
        )
    return atoms


def build_ligand_atom_payload(entry: dict, output_root: Path, model_index: int) -> dict:
    atom_metadata = find_ligand_atom_metadata(entry, output_root)
    coordinates = entry["aligned_atom_coordinates"]
    if len(atom_metadata) != len(coordinates):
        raise RuntimeError(
            f"Atom metadata count mismatch for {entry['ligand_instance_id']}: "
            f"metadata={len(atom_metadata)}, coords={len(coordinates)}"
        )

    atoms = []
    for atom_info, coord in zip(atom_metadata, coordinates):
        if atom_info["element"].upper() == "H":
            continue
        atoms.append(
            {
                "atom_name": atom_info["atom_name"],
                "element": atom_info["element"],
                "x": coord[0],
                "y": coord[1],
                "z": coord[2],
            }
        )

    return {
        "ligand_instance_id": entry["ligand_instance_id"],
        "ligand_code": entry["ligand_code"],
        "model_index": model_index + 1,
        "atom_count": len(atoms),
        "atoms": atoms,
        "centroid": [
            sum(atom[axis] for atom in atoms) / len(atoms)
            for axis in ["x", "y", "z"]
        ],
    }


def build_anchor_pdb_block(ligand_payloads: list[dict]) -> str:
    lines = []
    serial = 1
    residue_serial = 1
    for payload in ligand_payloads:
        for atom in payload["atoms"]:
            atom_name = atom["atom_name"][:4].rjust(4)
            res_name = payload["ligand_code"][:3].rjust(3)
            element = atom["element"][:2].rjust(2)
            lines.append(
                f"HETATM{serial:5d} {atom_name} {res_name} Z{residue_serial:4d}    "
                f"{atom['x']:8.3f}{atom['y']:8.3f}{atom['z']:8.3f}  1.00 20.00          {element}"
            )
            serial += 1
        residue_serial += 1
    lines.append("END")
    return "\n".join(lines) + "\n"


def generate_palette(count: int) -> list[str]:
    base = [
        "#16a34a",
        "#0891b2",
        "#c026d3",
        "#ca8a04",
        "#ea580c",
        "#7c3aed",
        "#2563eb",
        "#dc2626",
    ]
    colors = []
    for index in range(count):
        colors.append(base[index % len(base)])
    return colors


def build_py3dmol_script(gene_name: str, output_dir: str | Path, protein_cif_path: str | Path | None, include_alternative: bool) -> str:
    output_root = Path(output_dir)
    reference_protomer, aligned_ligands, primary_cluster, final_box = load_stage_artifacts(gene_name, output_root)
    selected_ligands = select_ligands(aligned_ligands, primary_cluster, include_alternative)

    default_protein_cif = output_root / reference_protomer["source_mmcif_file"]
    protein_cif = Path(protein_cif_path) if protein_cif_path else default_protein_cif
    palette = generate_palette(len(selected_ligands))

    ligand_payloads = []
    ligand_meta = []
    for index, entry in enumerate(selected_ligands):
        payload = build_ligand_atom_payload(entry, output_root, index)
        ligand_payloads.append(payload)
        ligand_meta.append(
            {
                "ligand_instance_id": entry["ligand_instance_id"],
                "ligand_code": entry["ligand_code"],
                "model_index": index + 1,
                "color": palette[index],
                "site_cluster_id": entry.get("site_cluster_id"),
                "cluster_role": entry.get("cluster_role"),
                "atom_count": payload["atom_count"],
            }
        )

    center_x, center_y, center_z = final_box["box_center"]
    size_x, size_y, size_z = final_box["box_size"]
    anchor_pdb_block = build_anchor_pdb_block(ligand_payloads)

    script_lines = [
        "from pathlib import Path",
        "import py3Dmol",
        "from IPython.display import HTML",
        "",
        f"GENE_NAME = {gene_name!r}",
        f"PROTEIN_CIF_PATH = Path({str(protein_cif)!r})",
        f"REFERENCE_PROTOMER_INSTANCE_ID = {reference_protomer['protomer_instance_id']!r}",
        f"PRIMARY_SITE_CLUSTER_ID = {primary_cluster['site_cluster_id']!r}",
        f"INCLUDE_ALTERNATIVE = {include_alternative!r}",
        f"EXPECTED_LIGAND_COUNT = {len(selected_ligands)}",
        f"LIGAND_META = {json.dumps(ligand_meta, ensure_ascii=False, indent=2)}",
        f"LIGAND_PAYLOADS = {json.dumps(ligand_payloads, ensure_ascii=False, indent=2)}",
        f"ANCHOR_PDB_BLOCK = {json.dumps(anchor_pdb_block, ensure_ascii=False)}",
        f"BOX_CENTER = {json.dumps([center_x, center_y, center_z])}",
        f"BOX_SIZE = {json.dumps([size_x, size_y, size_z])}",
        "",
        "def _base_radius(element):",
        "    element = element.upper()",
        "    if element in {'CL', 'BR', 'I'}:",
        "        return 1.2",
        "    if element in {'O', 'N', 'S', 'P', 'F'}:",
        "        return 0.95",
        "    return 0.85",
        "",
        "def build_view(protein_cif_path=None, sphere_scale=1.0, ligand_opacity=0.92, protein_opacity=0.22, box_opacity=0.25, show_labels=True, width=1100, height=800):",
        "    protein_path = Path(protein_cif_path) if protein_cif_path else PROTEIN_CIF_PATH",
        "    protein_data = protein_path.read_text(encoding='utf-8')",
        "    view = py3Dmol.view(width=width, height=height)",
        "    view.setBackgroundColor('white')",
        "    view.addModel(protein_data, 'mmcif')",
        "    view.setStyle({'model': 0, 'protein': True}, {'cartoon': {'color': 'lightgray', 'opacity': protein_opacity}})",
        "    view.setStyle({'model': 0, 'hetflag': True}, {'stick': {'hidden': True}})",
        "    view.addModel(ANCHOR_PDB_BLOCK, 'pdb')",
        "    view.setStyle({'model': 1}, {'sphere': {'radius': 0.05, 'opacity': 0.0}})",
        "    rendered_count = 0",
        "    for payload, meta in zip(LIGAND_PAYLOADS, LIGAND_META):",
        "        rendered_count += 1",
        "        for atom in payload['atoms']:",
        "            view.addSphere({",
        "                'center': {'x': atom['x'], 'y': atom['y'], 'z': atom['z']},",
        "                'radius': _base_radius(atom['element']) * sphere_scale,",
        "                'color': meta['color'],",
        "                'opacity': ligand_opacity",
        "            })",
        "        if show_labels and payload['atoms']:",
        "            first_atom = payload['atoms'][0]",
        "            view.addLabel(meta['ligand_code'], {",
        "                'position': {'x': first_atom['x'], 'y': first_atom['y'], 'z': first_atom['z']},",
        "                'backgroundColor': 'white',",
        "                'fontColor': 'black',",
        "                'fontSize': 12,",
        "                'showBackground': True,",
        "                'inFront': True",
        "            })",
        "    if rendered_count != EXPECTED_LIGAND_COUNT:",
        "        raise RuntimeError(f'Rendered {rendered_count} ligands but expected {EXPECTED_LIGAND_COUNT}')",
        "    view.addBox({",
        "        'center': {'x': BOX_CENTER[0], 'y': BOX_CENTER[1], 'z': BOX_CENTER[2]},",
        "        'dimensions': {'w': BOX_SIZE[0], 'h': BOX_SIZE[1], 'd': BOX_SIZE[2]},",
        "        'color': 'cyan',",
        "        'opacity': box_opacity",
        "    })",
        "    view.zoomTo({'model': 1})",
        "    view.zoom(1.8)",
        "    return view",
        "",
        "def show_view(**kwargs):",
        "    view = build_view(**kwargs)",
        "    return HTML(view._make_html())",
    ]
    return "\n".join(script_lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate a py3Dmol script for bouquet-style ligand cloud visualization."
    )
    parser.add_argument("gene_name", help="Gene symbol, for example ACHE")
    parser.add_argument(
        "--output-dir",
        default=str(Path(resolve_default_output_dir())),
        help="Base output directory containing pipeline artifacts.",
    )
    parser.add_argument(
        "--protein-cif",
        default=None,
        help=(
            "Optional path to a CIF file to use as the background protein model. "
            "This CIF must already be in the same aligned frame as the reference protomer."
        ),
    )
    parser.add_argument(
        "--include-alternative",
        action="store_true",
        help="Include ligands from alternative clusters in addition to the primary site cluster.",
    )
    parser.add_argument(
        "--write-path",
        default=None,
        help="Optional path for writing the generated Python script. Defaults to output/<gene>/07_exports/visualize_ligand_cloud_py3dmol.py.",
    )
    args = parser.parse_args()

    gene_name = normalize_gene_name(args.gene_name)
    default_write_path = Path(args.output_dir) / gene_name / "07_exports" / "visualize_ligand_cloud_py3dmol.py"
    write_path = Path(args.write_path) if args.write_path else default_write_path

    script_text = build_py3dmol_script(
        gene_name,
        args.output_dir,
        args.protein_cif,
        args.include_alternative,
    )
    write_path.parent.mkdir(parents=True, exist_ok=True)
    write_path.write_text(script_text, encoding="utf-8")
    print(write_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
