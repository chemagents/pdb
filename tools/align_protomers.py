import argparse
import json
import math
import shlex
import sys
from pathlib import Path

from rcsb_receptor_utils import ensure_mmcif_file
from rcsb_receptor_utils import normalize_gene_name
from rcsb_receptor_utils import resolve_default_output_dir
from rcsb_receptor_utils import write_json


def load_json(path: str | Path) -> dict | list:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def load_index(gene_name: str, output_dir: str | Path) -> dict:
    path = Path(output_dir) / gene_name / "02_protomers" / "index_protomers.json"
    if not path.exists():
        raise RuntimeError(
            f"Missing protomer index file at {path}. Alignment stage expects completed Stage 2 artifacts."
        )
    data = load_json(path)
    if not isinstance(data, dict):
        raise RuntimeError(f"Protomer index at {path} must contain a JSON object")
    return data


def load_protomer_records(gene_name: str, output_dir: str | Path) -> list[dict]:
    root = Path(output_dir) / gene_name / "02_protomers" / "per_structure"
    protomer_paths = sorted(root.glob("*/protomer_*.json"))
    protomer_paths = [path for path in protomer_paths if path.stem != "protomer_stage"]
    if not protomer_paths:
        raise RuntimeError(f"No protomer JSON files found under {root}")

    records = []
    for protomer_path in protomer_paths:
        data = load_json(protomer_path)
        if not isinstance(data, dict):
            raise RuntimeError(f"Protomer artifact at {protomer_path} must contain a JSON object")
        records.append(data)
    return records


def load_receptor_record(gene_name: str, structure_id: str, output_dir: str | Path) -> dict:
    path = Path(output_dir) / gene_name / f"receptor_structure_pdb_{structure_id}.json"
    if not path.exists():
        raise RuntimeError(f"Missing receptor_structure JSON for {structure_id} at {path}")
    data = load_json(path)
    if not isinstance(data, dict):
        raise RuntimeError(f"Receptor record at {path} must contain a JSON object")
    return data


def resolve_referents_dir() -> Path:
    return Path(__file__).resolve().parent.parent / "referents"


def discover_referent_backed_base_reference(gene_name: str, output_dir: str | Path) -> dict | None:
    referents_dir = resolve_referents_dir()
    if not referents_dir.exists():
        return None

    gene_dir = Path(output_dir) / gene_name
    receptor_paths = sorted(gene_dir.glob("receptor_structure_pdb_*.json"))
    if not receptor_paths:
        return None

    candidates = []
    for receptor_path in receptor_paths:
        record = load_json(receptor_path)
        if not isinstance(record, dict):
            continue

        structure_id = record.get("structure_id")
        chains = record.get("chains") or []
        oligomeric_state = record.get("oligomeric_state")
        referent_config = referents_dir / f"{str(structure_id).lower()}_config.txt"
        if not referent_config.exists():
            continue
        if oligomeric_state != "monomer" or not chains:
            continue

        candidates.append(
            {
                "structure_id": structure_id,
                "auth_asym_id": sorted(chains)[0],
                "resolution_angstrom": record.get("resolution_angstrom"),
                "referent_config": str(referent_config),
            }
        )

    if not candidates:
        return None

    candidates.sort(
        key=lambda item: (
            item["resolution_angstrom"] is None,
            item["resolution_angstrom"] or 999.0,
            item["structure_id"],
            item["auth_asym_id"],
        )
    )
    return candidates[0]


def build_external_reference_record(gene_name: str, output_dir: str | Path, structure_id: str, auth_asym_id: str) -> dict:
    record = load_receptor_record(gene_name, structure_id, output_dir)
    stage_dir = Path(output_dir) / gene_name / "02_protomers"
    source_structures_dir = stage_dir / "source_structures"
    source_structures_dir.mkdir(parents=True, exist_ok=True)
    mmcif_path = ensure_mmcif_file(structure_id, source_structures_dir / f"{structure_id}.cif")
    mmcif_relative_path = str(mmcif_path.relative_to(Path(output_dir))).replace("\\", "/")

    referent_config = resolve_referents_dir() / f"{str(structure_id).lower()}_config.txt"

    return {
        "protomer_instance_id": f"protomer_instance:{structure_id}:{auth_asym_id}",
        "structure_id": structure_id,
        "gene_name": gene_name,
        "organism": record.get("organism"),
        "uniprot_id": record.get("uniprot_id"),
        "assembly_type": record.get("assembly_type"),
        "oligomeric_state": record.get("oligomeric_state"),
        "stoichiometry": record.get("stoichiometry"),
        "auth_asym_id": auth_asym_id,
        "chain_rank_in_record": 0,
        "protomer_mode": "external_base_reference",
        "selection_status": "external_reference",
        "stage_status": "accepted_reference_candidate",
        "relevance": "fixed_base_reference",
        "primary_ligand_code": record.get("primary_ligand_code"),
        "bound_ligands": record.get("bound_ligands") or [],
        "ligand_state": record.get("ligand_state"),
        "experimental_method": record.get("experimental_method"),
        "resolution_angstrom": record.get("resolution_angstrom"),
        "source_receptor_record": f"receptor_structure_pdb_{structure_id}.json",
        "source_mmcif_file": mmcif_relative_path,
        "referent_config": str(referent_config) if referent_config.exists() else None,
        "notes": record.get("notes"),
    }


def resolve_reference_protomer(gene_name: str, index_data: dict, protomer_records: list[dict]) -> dict:
    base_protomer = discover_referent_backed_base_reference(gene_name, index_data["output_dir"])
    if base_protomer is not None:
        for record in protomer_records:
            if (
                record.get("structure_id") == base_protomer["structure_id"]
                and record.get("auth_asym_id") == base_protomer["auth_asym_id"]
            ):
                resolved_record = dict(record)
                resolved_record["referent_config"] = base_protomer["referent_config"]
                return resolved_record
        return build_external_reference_record(
            gene_name,
            index_data["output_dir"],
            base_protomer["structure_id"],
            base_protomer["auth_asym_id"],
        )

    reference_id = index_data.get("canonical_reference_candidate_id")
    if not reference_id:
        raise RuntimeError("Protomer index is missing canonical_reference_candidate_id")

    for record in protomer_records:
        if record.get("protomer_instance_id") == reference_id:
            return record
    raise RuntimeError(f"Unable to resolve reference protomer {reference_id}")


def parse_loop_row(line: str) -> list[str]:
    lexer = shlex.shlex(line, posix=True)
    lexer.whitespace_split = True
    lexer.commenters = ""
    return list(lexer)


def parse_mmcif_ca_coords(mmcif_path: str | Path, auth_asym_id: str) -> dict[int, tuple[float, float, float]]:
    path = Path(mmcif_path)
    lines = path.read_text(encoding="utf-8").splitlines()

    in_atom_site_loop = False
    atom_site_fields = []
    field_index = {}
    coords = {}

    required_fields = {
        "_atom_site.auth_asym_id",
        "_atom_site.auth_seq_id",
        "_atom_site.label_atom_id",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
    }

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
                    in_atom_site_loop = True
                    atom_site_fields = candidate_fields
                    field_index = {name: idx for idx, name in enumerate(atom_site_fields)}
                    if not required_fields.issubset(field_index):
                        raise RuntimeError(
                            f"mmCIF file {path} is missing required atom_site fields for alignment"
                        )
                    i = j
                    continue
            i += 1
            continue

        if line == "#":
            break

        if line == "loop_" or line.startswith("_"):
            break

        row = parse_loop_row(line)
        if len(row) != len(atom_site_fields):
            i += 1
            continue

        chain_value = row[field_index["_atom_site.auth_asym_id"]]
        atom_name = row[field_index["_atom_site.label_atom_id"]]
        auth_seq_id = row[field_index["_atom_site.auth_seq_id"]]
        if chain_value != auth_asym_id or atom_name != "CA" or auth_seq_id in {"?", "."}:
            i += 1
            continue

        try:
            seq_id = int(auth_seq_id)
            x = float(row[field_index["_atom_site.Cartn_x"]])
            y = float(row[field_index["_atom_site.Cartn_y"]])
            z = float(row[field_index["_atom_site.Cartn_z"]])
        except ValueError:
            i += 1
            continue

        coords[seq_id] = (x, y, z)
        i += 1

    if not coords:
        raise RuntimeError(f"No C-alpha coordinates found in {path} for chain {auth_asym_id}")
    return coords


def parse_referent_box(config_path: str | Path) -> dict | None:
    path = Path(config_path)
    if not path.exists():
        return None

    values = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = [part.strip() for part in line.split("=", 1)]
        values[key] = value

    required = {"center_x", "center_y", "center_z", "size_x", "size_y", "size_z"}
    if not required.issubset(values):
        return None

    return {
        "center_x": float(values["center_x"]),
        "center_y": float(values["center_y"]),
        "center_z": float(values["center_z"]),
        "size_x": float(values["size_x"]),
        "size_y": float(values["size_y"]),
        "size_z": float(values["size_z"]),
    }


def is_point_inside_expanded_box(point: tuple[float, float, float], box: dict, padding: float = 8.0) -> bool:
    x, y, z = point
    return (
        (box["center_x"] - box["size_x"] / 2.0 - padding) <= x <= (box["center_x"] + box["size_x"] / 2.0 + padding)
        and (box["center_y"] - box["size_y"] / 2.0 - padding) <= y <= (box["center_y"] + box["size_y"] / 2.0 + padding)
        and (box["center_z"] - box["size_z"] / 2.0 - padding) <= z <= (box["center_z"] + box["size_z"] / 2.0 + padding)
    )


def filter_reference_coords_by_referent_box(reference_coords: dict[int, tuple[float, float, float]], referent_box: dict | None) -> dict[int, tuple[float, float, float]]:
    if referent_box is None:
        return reference_coords

    filtered = {
        seq_id: point
        for seq_id, point in reference_coords.items()
        if is_point_inside_expanded_box(point, referent_box)
    }
    return filtered or reference_coords


def shared_ca_points(
    reference_coords: dict[int, tuple[float, float, float]],
    target_coords: dict[int, tuple[float, float, float]],
) -> tuple[list[tuple[float, float, float]], list[tuple[float, float, float]], list[int]]:
    shared_seq_ids = sorted(set(reference_coords) & set(target_coords))
    reference_points = [reference_coords[seq_id] for seq_id in shared_seq_ids]
    target_points = [target_coords[seq_id] for seq_id in shared_seq_ids]
    return reference_points, target_points, shared_seq_ids


def centroid(points: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    count = len(points)
    return (
        sum(point[0] for point in points) / count,
        sum(point[1] for point in points) / count,
        sum(point[2] for point in points) / count,
    )


def subtract(point: tuple[float, float, float], other: tuple[float, float, float]) -> tuple[float, float, float]:
    return (point[0] - other[0], point[1] - other[1], point[2] - other[2])


def add(point: tuple[float, float, float], other: tuple[float, float, float]) -> tuple[float, float, float]:
    return (point[0] + other[0], point[1] + other[1], point[2] + other[2])


def mat_vec_mul(matrix: list[list[float]], vector: tuple[float, float, float]) -> tuple[float, float, float]:
    return (
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2],
        matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2],
    )


def determinant_3x3(matrix: list[list[float]]) -> float:
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def quaternion_to_rotation_matrix(quaternion: tuple[float, float, float, float]) -> list[list[float]]:
    q0, q1, q2, q3 = quaternion
    return [
        [
            q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3,
            2.0 * (q1 * q2 - q0 * q3),
            2.0 * (q1 * q3 + q0 * q2),
        ],
        [
            2.0 * (q1 * q2 + q0 * q3),
            q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3,
            2.0 * (q2 * q3 - q0 * q1),
        ],
        [
            2.0 * (q1 * q3 - q0 * q2),
            2.0 * (q2 * q3 + q0 * q1),
            q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3,
        ],
    ]


def normalize_quaternion(quaternion: tuple[float, float, float, float]) -> tuple[float, float, float, float]:
    norm = math.sqrt(sum(component * component for component in quaternion))
    if norm == 0.0:
        return (1.0, 0.0, 0.0, 0.0)
    return (
        quaternion[0] / norm,
        quaternion[1] / norm,
        quaternion[2] / norm,
        quaternion[3] / norm,
    )


def dominant_quaternion_from_profile(profile: list[list[float]]) -> tuple[float, float, float, float]:
    vector = (1.0, 0.0, 0.0, 0.0)
    for _ in range(50):
        next_vector = []
        for row in profile:
            next_vector.append(sum(row[idx] * vector[idx] for idx in range(4)))
        vector = normalize_quaternion(tuple(next_vector))
    return vector


def kabsch_via_quaternion(
    reference_points: list[tuple[float, float, float]],
    target_points: list[tuple[float, float, float]],
) -> tuple[list[list[float]], tuple[float, float, float], float]:
    reference_centroid = centroid(reference_points)
    target_centroid = centroid(target_points)

    centered_reference = [subtract(point, reference_centroid) for point in reference_points]
    centered_target = [subtract(point, target_centroid) for point in target_points]

    sxx = syy = szz = sxy = sxz = syx = syz = szx = szy = 0.0
    for target_point, reference_point in zip(centered_target, centered_reference):
        tx, ty, tz = target_point
        rx, ry, rz = reference_point
        sxx += tx * rx
        sxy += tx * ry
        sxz += tx * rz
        syx += ty * rx
        syy += ty * ry
        syz += ty * rz
        szx += tz * rx
        szy += tz * ry
        szz += tz * rz

    profile = [
        [sxx + syy + szz, syz - szy, szx - sxz, sxy - syx],
        [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
        [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
        [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz],
    ]

    quaternion = dominant_quaternion_from_profile(profile)
    rotation = quaternion_to_rotation_matrix(quaternion)
    translation = subtract(reference_centroid, mat_vec_mul(rotation, target_centroid))

    squared_error_sum = 0.0
    for target_point, reference_point in zip(target_points, reference_points):
        transformed_point = add(mat_vec_mul(rotation, target_point), translation)
        dx = transformed_point[0] - reference_point[0]
        dy = transformed_point[1] - reference_point[1]
        dz = transformed_point[2] - reference_point[2]
        squared_error_sum += dx * dx + dy * dy + dz * dz

    rmsd = math.sqrt(squared_error_sum / len(target_points))
    return rotation, translation, rmsd


def build_alignment_record(
    protomer_record: dict,
    reference_record: dict,
    shared_seq_ids: list[int],
    rotation: list[list[float]],
    translation: tuple[float, float, float],
    rmsd: float,
) -> dict:
    determinant = determinant_3x3(rotation)
    status = "accepted" if len(shared_seq_ids) >= 30 and rmsd <= 5.0 else "requires_review"

    warnings = []
    if len(shared_seq_ids) < 30:
        warnings.append(f"low_alignment_coverage:{len(shared_seq_ids)}")
    if rmsd > 5.0:
        warnings.append(f"high_rmsd:{round(rmsd, 4)}")
    if abs(determinant - 1.0) > 0.05:
        warnings.append(f"rotation_determinant_suspicious:{round(determinant, 4)}")

    return {
        "protomer_instance_id": protomer_record["protomer_instance_id"],
        "structure_id": protomer_record["structure_id"],
        "auth_asym_id": protomer_record["auth_asym_id"],
        "reference_protomer_instance_id": reference_record["protomer_instance_id"],
        "reference_structure_id": reference_record["structure_id"],
        "reference_auth_asym_id": reference_record["auth_asym_id"],
        "source_mmcif_file": protomer_record["source_mmcif_file"],
        "reference_mmcif_file": reference_record["source_mmcif_file"],
        "fit_residue_count": len(shared_seq_ids),
        "fit_sequence_range": {
            "start": shared_seq_ids[0] if shared_seq_ids else None,
            "end": shared_seq_ids[-1] if shared_seq_ids else None,
        },
        "rmsd_backbone_ca": rmsd,
        "rotation_matrix": rotation,
        "translation_vector": [translation[0], translation[1], translation[2]],
        "rotation_determinant": determinant,
        "status": status,
        "warnings": warnings,
    }


def build_alignment_summary(gene_name: str, alignment_records: list[dict], reference_record: dict) -> dict:
    accepted = [item for item in alignment_records if item["status"] == "accepted"]
    requires_review = [item for item in alignment_records if item["status"] != "accepted"]

    rmsd_values = [item["rmsd_backbone_ca"] for item in alignment_records if item["rmsd_backbone_ca"] is not None]
    if rmsd_values:
        min_rmsd = min(rmsd_values)
        max_rmsd = max(rmsd_values)
        mean_rmsd = sum(rmsd_values) / len(rmsd_values)
    else:
        min_rmsd = max_rmsd = mean_rmsd = None

    return {
        "gene_name": gene_name,
        "pipeline_stage": 3,
        "stage_name": "canonical_alignment",
        "reference_protomer_instance_id": reference_record["protomer_instance_id"],
        "reference_structure_id": reference_record["structure_id"],
        "reference_auth_asym_id": reference_record["auth_asym_id"],
        "alignment_record_count": len(alignment_records),
        "accepted_alignment_count": len(accepted),
        "requires_review_alignment_count": len(requires_review),
        "requires_review_protomer_instance_ids": [item["protomer_instance_id"] for item in requires_review],
        "rmsd_backbone_ca_min": min_rmsd,
        "rmsd_backbone_ca_mean": mean_rmsd,
        "rmsd_backbone_ca_max": max_rmsd,
    }


def write_alignment_report(summary: dict, alignment_records: list[dict], output_path: str | Path) -> None:
    lines = [
        f"# Alignment Summary: {summary['gene_name']}",
        "",
        f"- Reference protomer: `{summary['reference_protomer_instance_id']}`",
        f"- Alignment records: {summary['alignment_record_count']}",
        f"- Accepted: {summary['accepted_alignment_count']}",
        f"- Requires review: {summary['requires_review_alignment_count']}",
        f"- RMSD min: {summary['rmsd_backbone_ca_min']}",
        f"- RMSD mean: {summary['rmsd_backbone_ca_mean']}",
        f"- RMSD max: {summary['rmsd_backbone_ca_max']}",
        "",
        "## Requires Review",
        "",
    ]

    for item in alignment_records:
        if item["status"] != "accepted":
            lines.append(f"- `{item['protomer_instance_id']}`: {', '.join(item['warnings'])}")

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Align protomer instances to a canonical reference protomer using local mmCIF coordinate artifacts."
        )
    )
    parser.add_argument("gene_name", help="Gene symbol, for example ACHE")
    parser.add_argument(
        "--output-dir",
        default=str(resolve_default_output_dir()),
        help=(
            "Base output directory. Inputs are read from output/<gene_name>/02_protomers/ and "
            "alignment artifacts are written into output/<gene_name>/03_alignment/."
        ),
    )
    args = parser.parse_args()

    try:
        gene_name = normalize_gene_name(args.gene_name)
        index_data = load_index(gene_name, args.output_dir)
        index_data["output_dir"] = args.output_dir
        protomer_records = load_protomer_records(gene_name, args.output_dir)
        reference_record = resolve_reference_protomer(gene_name, index_data, protomer_records)

        reference_mmcif_path = Path(args.output_dir) / reference_record["source_mmcif_file"]
        reference_coords = parse_mmcif_ca_coords(reference_mmcif_path, reference_record["auth_asym_id"])
        referent_config = reference_record.get("referent_config")
        referent_box = parse_referent_box(referent_config) if isinstance(referent_config, str) and referent_config else None
        reference_coords = filter_reference_coords_by_referent_box(reference_coords, referent_box)

        stage_dir = Path(args.output_dir) / gene_name / "03_alignment"
        per_structure_dir = stage_dir / "per_structure"
        reports_dir = stage_dir / "reports"
        stage_dir.mkdir(parents=True, exist_ok=True)
        per_structure_dir.mkdir(parents=True, exist_ok=True)
        reports_dir.mkdir(parents=True, exist_ok=True)

        alignment_records = []
        for protomer_record in protomer_records:
            target_mmcif_path = Path(args.output_dir) / protomer_record["source_mmcif_file"]
            target_coords = parse_mmcif_ca_coords(target_mmcif_path, protomer_record["auth_asym_id"])
            reference_points, target_points, shared_seq_ids = shared_ca_points(reference_coords, target_coords)
            if len(shared_seq_ids) < 3:
                alignment_record = {
                    "protomer_instance_id": protomer_record["protomer_instance_id"],
                    "structure_id": protomer_record["structure_id"],
                    "auth_asym_id": protomer_record["auth_asym_id"],
                    "reference_protomer_instance_id": reference_record["protomer_instance_id"],
                    "reference_structure_id": reference_record["structure_id"],
                    "reference_auth_asym_id": reference_record["auth_asym_id"],
                    "source_mmcif_file": protomer_record["source_mmcif_file"],
                    "reference_mmcif_file": reference_record["source_mmcif_file"],
                    "fit_residue_count": len(shared_seq_ids),
                    "fit_sequence_range": {
                        "start": shared_seq_ids[0] if shared_seq_ids else None,
                        "end": shared_seq_ids[-1] if shared_seq_ids else None,
                    },
                    "rmsd_backbone_ca": None,
                    "rotation_matrix": None,
                    "translation_vector": None,
                    "rotation_determinant": None,
                    "status": "requires_review",
                    "warnings": [f"insufficient_shared_residues:{len(shared_seq_ids)}"],
                }
            else:
                rotation, translation, rmsd = kabsch_via_quaternion(reference_points, target_points)
                alignment_record = build_alignment_record(
                    protomer_record,
                    reference_record,
                    shared_seq_ids,
                    rotation,
                    translation,
                    rmsd,
                )

            alignment_records.append(alignment_record)

            structure_dir = per_structure_dir / protomer_record["structure_id"]
            structure_dir.mkdir(parents=True, exist_ok=True)
            output_name = f"{protomer_record['auth_asym_id']}.align.json"
            write_json(alignment_record, structure_dir / output_name)

        summary = build_alignment_summary(gene_name, alignment_records, reference_record)
        transform_index = {
            record["protomer_instance_id"]: {
                "rotation_matrix": record["rotation_matrix"],
                "translation_vector": record["translation_vector"],
                "status": record["status"],
            }
            for record in alignment_records
        }

        write_json(reference_record, stage_dir / "reference_protomer.json")
        write_json(alignment_records, stage_dir / "alignment_records.json")
        write_json(transform_index, stage_dir / "transforms.json")
        write_json(summary, stage_dir / "alignment_qc.json")
        write_alignment_report(summary, alignment_records, reports_dir / "alignment_summary.md")

        print(stage_dir / "alignment_qc.json")
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"Failed to align protomers for {args.gene_name}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
