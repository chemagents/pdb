import argparse
import json
import math
import shlex
import sys
from pathlib import Path

from rcsb_receptor_utils import normalize_gene_name
from rcsb_receptor_utils import resolve_default_output_dir
from rcsb_receptor_utils import write_json


CLUSTER_DISTANCE_THRESHOLD = 6.0


def load_json(path: str | Path) -> dict | list:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def load_ligand_instances(gene_name: str, output_dir: str | Path) -> list[dict]:
    path = Path(output_dir) / gene_name / "02_protomers" / "accepted_ligand_instances.json"
    if not path.exists():
        raise RuntimeError(f"Missing Stage 2 ligand instance file at {path}")
    data = load_json(path)
    if not isinstance(data, list):
        raise RuntimeError(f"Ligand instance file at {path} must contain a JSON list")
    return data


def load_alignment_records(gene_name: str, output_dir: str | Path) -> list[dict]:
    path = Path(output_dir) / gene_name / "03_alignment" / "alignment_records.json"
    if not path.exists():
        raise RuntimeError(f"Missing Stage 3 alignment record file at {path}")
    data = load_json(path)
    if not isinstance(data, list):
        raise RuntimeError(f"Alignment record file at {path} must contain a JSON list")
    return data


def load_reference_protomer(gene_name: str, output_dir: str | Path) -> dict:
    path = Path(output_dir) / gene_name / "03_alignment" / "reference_protomer.json"
    if not path.exists():
        raise RuntimeError(f"Missing Stage 3 reference protomer file at {path}")
    data = load_json(path)
    if not isinstance(data, dict):
        raise RuntimeError(f"Reference protomer file at {path} must contain a JSON object")
    return data


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


def centroid(points: list[tuple[float, float, float]]) -> tuple[float, float, float]:
    count = len(points)
    return (
        sum(point[0] for point in points) / count,
        sum(point[1] for point in points) / count,
        sum(point[2] for point in points) / count,
    )


def euclidean_distance(point_a: tuple[float, float, float], point_b: tuple[float, float, float]) -> float:
    dx = point_a[0] - point_b[0]
    dy = point_a[1] - point_b[1]
    dz = point_a[2] - point_b[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def mat_vec_mul(matrix: list[list[float]], vector: tuple[float, float, float]) -> tuple[float, float, float]:
    return (
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2],
        matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2],
    )


def add(point: tuple[float, float, float], other: tuple[float, float, float]) -> tuple[float, float, float]:
    return (point[0] + other[0], point[1] + other[1], point[2] + other[2])


def extract_chain_ca_centroid(mmcif_path: str | Path, auth_asym_id: str) -> tuple[float, float, float]:
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
        raise RuntimeError(f"No C-alpha points found for chain {auth_asym_id} in {mmcif_path}")
    return centroid(points)


def extract_ligand_groups(mmcif_path: str | Path, ligand_code: str) -> list[dict]:
    groups = {}

    for row, field_index in iter_atom_site_rows(mmcif_path):
        auth_comp_id = row[field_index.get("_atom_site.auth_comp_id")]
        label_comp_id = row[field_index.get("_atom_site.label_comp_id")]
        if ligand_code not in {auth_comp_id, label_comp_id}:
            continue

        try:
            x = float(row[field_index["_atom_site.Cartn_x"]])
            y = float(row[field_index["_atom_site.Cartn_y"]])
            z = float(row[field_index["_atom_site.Cartn_z"]])
        except ValueError:
            continue

        auth_asym_id = row[field_index["_atom_site.auth_asym_id"]]
        auth_seq_id = row[field_index["_atom_site.auth_seq_id"]]
        atom_name = row[field_index["_atom_site.label_atom_id"]]
        key = (auth_asym_id, auth_seq_id)
        groups.setdefault(
            key,
            {
                "auth_asym_id": auth_asym_id,
                "auth_seq_id": auth_seq_id,
                "atom_names": [],
                "atom_coordinates": [],
            },
        )
        groups[key]["atom_names"].append(atom_name)
        groups[key]["atom_coordinates"].append((x, y, z))

    result = []
    for group in groups.values():
        group["centroid"] = centroid(group["atom_coordinates"])
        result.append(group)
    return result


def choose_best_ligand_group(
    ligand_groups: list[dict],
    chain_centroid: tuple[float, float, float],
) -> dict:
    if not ligand_groups:
        raise RuntimeError("No ligand groups available for selection")

    ranked = sorted(
        ligand_groups,
        key=lambda group: (
            euclidean_distance(group["centroid"], chain_centroid),
            -len(group["atom_coordinates"]),
            group["auth_asym_id"],
            group["auth_seq_id"],
        ),
    )
    return ranked[0]


def transform_points(
    points: list[tuple[float, float, float]],
    rotation_matrix: list[list[float]],
    translation_vector: list[float],
) -> list[tuple[float, float, float]]:
    translation = (
        translation_vector[0],
        translation_vector[1],
        translation_vector[2],
    )
    return [add(mat_vec_mul(rotation_matrix, point), translation) for point in points]


def build_aligned_ligand_entry(ligand_instance: dict, alignment_record: dict, ligand_group: dict) -> dict:
    aligned_points = transform_points(
        ligand_group["atom_coordinates"],
        alignment_record["rotation_matrix"],
        alignment_record["translation_vector"],
    )
    aligned_centroid = centroid(aligned_points)

    return {
        "ligand_instance_id": ligand_instance["ligand_instance_id"],
        "structure_id": ligand_instance["structure_id"],
        "ligand_code": ligand_instance["ligand_code"],
        "assigned_protomer_instance_id": ligand_instance["assigned_protomer_instance_id"],
        "alignment_status": alignment_record["status"],
        "source_mmcif_file": ligand_instance["source_mmcif_file"],
        "source_ligand_locator": {
            "auth_asym_id": ligand_group["auth_asym_id"],
            "auth_seq_id": ligand_group["auth_seq_id"],
        },
        "native_centroid": [
            ligand_group["centroid"][0],
            ligand_group["centroid"][1],
            ligand_group["centroid"][2],
        ],
        "aligned_centroid": [
            aligned_centroid[0],
            aligned_centroid[1],
            aligned_centroid[2],
        ],
        "atom_count": len(aligned_points),
        "aligned_atom_coordinates": [[point[0], point[1], point[2]] for point in aligned_points],
        "is_primary_ligand": ligand_instance.get("is_primary_ligand"),
        "ligand_state": ligand_instance.get("ligand_state"),
        "notes": ligand_instance.get("notes"),
    }


def cluster_aligned_ligands(aligned_entries: list[dict]) -> list[dict]:
    clusters = []
    for entry in aligned_entries:
        entry_centroid = tuple(entry["aligned_centroid"])
        assigned = False
        for cluster in clusters:
            if euclidean_distance(entry_centroid, tuple(cluster["cluster_centroid"])) <= CLUSTER_DISTANCE_THRESHOLD:
                cluster["members"].append(entry)
                points = [tuple(member["aligned_centroid"]) for member in cluster["members"]]
                cluster["cluster_centroid"] = list(centroid(points))
                assigned = True
                break
        if not assigned:
            clusters.append(
                {
                    "members": [entry],
                    "cluster_centroid": list(entry["aligned_centroid"]),
                }
            )

    normalized_clusters = []
    for index, cluster in enumerate(clusters, start=1):
        members = cluster["members"]
        structure_ids = sorted({member["structure_id"] for member in members})
        atom_points = []
        for member in members:
            atom_points.extend(tuple(point) for point in member["aligned_atom_coordinates"])
        normalized_clusters.append(
            {
                "site_cluster_id": f"site_cluster_{index}",
                "member_ligand_instance_ids": [member["ligand_instance_id"] for member in members],
                "supporting_structure_ids": structure_ids,
                "member_count": len(members),
                "cluster_centroid": cluster["cluster_centroid"],
                "atom_point_count": len(atom_points),
            }
        )

    normalized_clusters.sort(
        key=lambda cluster: (
            -cluster["member_count"],
            -len(cluster["supporting_structure_ids"]),
            cluster["site_cluster_id"],
        )
    )

    primary_cluster_id = normalized_clusters[0]["site_cluster_id"] if normalized_clusters else None
    for cluster in normalized_clusters:
        if cluster["site_cluster_id"] == primary_cluster_id:
            cluster["cluster_role"] = "primary"
        elif cluster["member_count"] == 1 and normalized_clusters and normalized_clusters[0]["member_count"] > 1:
            cluster["cluster_role"] = "outlier"
        else:
            cluster["cluster_role"] = "alternative"
    return normalized_clusters


def annotate_cluster_members(aligned_entries: list[dict], clusters: list[dict]) -> tuple[list[dict], list[dict]]:
    cluster_map = {}
    for cluster in clusters:
        for ligand_instance_id in cluster["member_ligand_instance_ids"]:
            cluster_map[ligand_instance_id] = cluster

    annotated = []
    outliers = []
    for entry in aligned_entries:
        cluster = cluster_map[entry["ligand_instance_id"]]
        annotated_entry = dict(entry)
        annotated_entry["site_cluster_id"] = cluster["site_cluster_id"]
        annotated_entry["cluster_role"] = cluster["cluster_role"]
        annotated.append(annotated_entry)
        if cluster["cluster_role"] == "outlier":
            outliers.append(annotated_entry)
    return annotated, outliers


def build_primary_cluster(clusters: list[dict]) -> dict | None:
    for cluster in clusters:
        if cluster["cluster_role"] == "primary":
            return cluster
    return None


def build_stage_summary(gene_name: str, aligned_entries: list[dict], clusters: list[dict], outliers: list[dict], reference_protomer: dict) -> dict:
    primary_cluster = build_primary_cluster(clusters)
    return {
        "gene_name": gene_name,
        "pipeline_stage": 4,
        "stage_name": "binding_site_clustering",
        "reference_protomer_instance_id": reference_protomer["protomer_instance_id"],
        "aligned_ligand_instance_count": len(aligned_entries),
        "cluster_count": len(clusters),
        "outlier_count": len(outliers),
        "primary_site_cluster_id": primary_cluster["site_cluster_id"] if primary_cluster else None,
        "primary_site_member_count": primary_cluster["member_count"] if primary_cluster else None,
        "primary_site_supporting_structure_ids": primary_cluster["supporting_structure_ids"] if primary_cluster else [],
    }


def write_site_report(summary: dict, clusters: list[dict], output_path: str | Path) -> None:
    lines = [
        f"# Site Cluster Summary: {summary['gene_name']}",
        "",
        f"- Reference protomer: `{summary['reference_protomer_instance_id']}`",
        f"- Aligned ligand instances: {summary['aligned_ligand_instance_count']}",
        f"- Clusters: {summary['cluster_count']}",
        f"- Outliers: {summary['outlier_count']}",
        f"- Primary site cluster: `{summary['primary_site_cluster_id']}`",
        "",
        "## Clusters",
        "",
    ]

    for cluster in clusters:
        lines.append(
            f"- `{cluster['site_cluster_id']}`: role={cluster['cluster_role']}, "
            f"members={cluster['member_count']}, structures={cluster['supporting_structure_ids']}"
        )

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Extract aligned ligand clouds and cluster binding-site candidates using Stage 2 and Stage 3 artifacts."
        )
    )
    parser.add_argument("gene_name", help="Gene symbol, for example ACHE")
    parser.add_argument(
        "--output-dir",
        default=str(resolve_default_output_dir()),
        help=(
            "Base output directory. Inputs are read from output/<gene_name>/02_protomers/ and "
            "output/<gene_name>/03_alignment/, and site artifacts are written into output/<gene_name>/04_sites/."
        ),
    )
    args = parser.parse_args()

    try:
        gene_name = normalize_gene_name(args.gene_name)
        ligand_instances = load_ligand_instances(gene_name, args.output_dir)
        alignment_records = load_alignment_records(gene_name, args.output_dir)
        reference_protomer = load_reference_protomer(gene_name, args.output_dir)

        alignment_by_protomer_id = {
            record["protomer_instance_id"]: record for record in alignment_records if record["status"] == "accepted"
        }

        aligned_entries = []
        for ligand_instance in ligand_instances:
            protomer_instance_id = ligand_instance["assigned_protomer_instance_id"]
            alignment_record = alignment_by_protomer_id.get(protomer_instance_id)
            if alignment_record is None:
                continue

            mmcif_path = Path(args.output_dir) / ligand_instance["source_mmcif_file"]
            chain_centroid = extract_chain_ca_centroid(mmcif_path, protomer_instance_id.split(":")[-1])
            ligand_groups = extract_ligand_groups(mmcif_path, ligand_instance["ligand_code"])
            if not ligand_groups:
                continue
            best_group = choose_best_ligand_group(ligand_groups, chain_centroid)
            aligned_entries.append(build_aligned_ligand_entry(ligand_instance, alignment_record, best_group))

        if not aligned_entries:
            raise RuntimeError(f"No aligned ligand entries could be extracted for {gene_name}")

        clusters = cluster_aligned_ligands(aligned_entries)
        annotated_entries, outliers = annotate_cluster_members(aligned_entries, clusters)
        primary_cluster = build_primary_cluster(clusters)
        summary = build_stage_summary(gene_name, annotated_entries, clusters, outliers, reference_protomer)

        stage_dir = Path(args.output_dir) / gene_name / "04_sites"
        reports_dir = stage_dir / "reports"
        stage_dir.mkdir(parents=True, exist_ok=True)
        reports_dir.mkdir(parents=True, exist_ok=True)

        write_json(annotated_entries, stage_dir / "accepted_ligand_instances_aligned.json")
        write_json(outliers, stage_dir / "rejected_ligand_instances.json")
        write_json(clusters, stage_dir / "site_clusters.json")
        if primary_cluster is not None:
            write_json(primary_cluster, stage_dir / "primary_site_cluster.json")
        write_json(summary, stage_dir / "site_cluster_qc.json")
        write_site_report(summary, clusters, reports_dir / "site_cluster_summary.md")

        print(stage_dir / "site_cluster_qc.json")
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"Failed to cluster binding sites for {args.gene_name}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
