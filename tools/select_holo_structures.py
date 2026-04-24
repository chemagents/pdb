import argparse
import json
import sys
from pathlib import Path

from rcsb_receptor_utils import normalize_gene_name
from rcsb_receptor_utils import resolve_default_output_dir
from rcsb_receptor_utils import write_json


def classify_structure(record: dict) -> tuple[str, list[str]]:
    reasons = []

    organism = record.get("organism")
    if organism not in {"Homo sapiens", "Mus musculus"}:
        reasons.append(f"unsupported_organism:{organism}")

    ligand_state = record.get("ligand_state")
    if ligand_state == "apo":
        reasons.append("apo_record")

    bound_ligands = record.get("bound_ligands") or []
    if not bound_ligands:
        reasons.append("no_meaningful_bound_ligands")

    notes = record.get("notes") or ""
    lowered_notes = notes.lower()
    if "peptide/protein-bound inhibitory context" in lowered_notes:
        reasons.append("peptide_or_protein_bound_context")

    if ligand_state == "unknown":
        reasons.append("unknown_ligand_state")

    if ligand_state == "mixed_complex":
        reasons.append("mixed_ligand_context")

    if len(bound_ligands) >= 6:
        reasons.append(f"many_bound_ligands:{len(bound_ligands)}")

    primary_ligand_code = record.get("primary_ligand_code")
    if bound_ligands and primary_ligand_code is None:
        reasons.append("no_primary_ligand_code")

    experimental_method = record.get("experimental_method")
    if experimental_method in {None, "OTHER"}:
        reasons.append("weak_experimental_method")

    resolution_angstrom = record.get("resolution_angstrom")
    if resolution_angstrom is None:
        reasons.append("missing_resolution")
    elif isinstance(resolution_angstrom, (int, float)) and resolution_angstrom > 3.5:
        reasons.append(f"low_resolution:{resolution_angstrom}")

    if any(reason.startswith("unsupported_organism") for reason in reasons):
        return "rejected", reasons

    if any(
        reason in {
            "apo_record",
            "no_meaningful_bound_ligands",
            "peptide_or_protein_bound_context",
        }
        for reason in reasons
    ):
        return "rejected", reasons

    if any(reason.startswith("low_resolution") for reason in reasons):
        return "uncertain", reasons

    if any(
        reason in {"missing_resolution", "weak_experimental_method"}
        for reason in reasons
    ):
        return "uncertain", reasons

    if any(reason.startswith("many_bound_ligands") for reason in reasons):
        return "uncertain", reasons

    if any(
        reason in {"no_primary_ligand_code", "unknown_ligand_state", "mixed_ligand_context"}
        for reason in reasons
    ) and bound_ligands:
        reasons.append("accepted_with_ligand_ambiguity")
        return "accepted", reasons

    return "accepted", ["usable_holo_structure"]


def build_selection_entry(record: dict) -> dict:
    status, reasons = classify_structure(record)
    return {
        "structure_id": record.get("structure_id"),
        "gene_name": record.get("gene_name"),
        "organism": record.get("organism"),
        "experimental_method": record.get("experimental_method"),
        "resolution_angstrom": record.get("resolution_angstrom"),
        "ligand_state": record.get("ligand_state"),
        "bound_ligands": record.get("bound_ligands"),
        "primary_ligand_code": record.get("primary_ligand_code"),
        "status": status,
        "reasons": reasons,
        "notes": record.get("notes"),
    }


def build_selection_summary(gene_name: str, decisions: list[dict]) -> dict:
    accepted = [entry["structure_id"] for entry in decisions if entry["status"] == "accepted"]
    uncertain = [entry["structure_id"] for entry in decisions if entry["status"] == "uncertain"]
    rejected = [entry["structure_id"] for entry in decisions if entry["status"] == "rejected"]

    reason_counts = {}
    for entry in decisions:
        for reason in entry["reasons"]:
            reason_counts[reason] = reason_counts.get(reason, 0) + 1

    return {
        "gene_name": gene_name,
        "pipeline_stage": 1,
        "stage_name": "structure_selection",
        "input_source": f"output/{gene_name}/receptor_structure_pdb_*.json",
        "input_structure_count": len(decisions),
        "accepted_structure_ids": accepted,
        "uncertain_structure_ids": uncertain,
        "rejected_structure_ids": rejected,
        "accepted_count": len(accepted),
        "uncertain_count": len(uncertain),
        "rejected_count": len(rejected),
        "reason_counts": reason_counts,
    }


def write_selection_report(summary: dict, decisions: list[dict], output_path: str | Path) -> None:
    accepted = [entry for entry in decisions if entry["status"] == "accepted"]
    uncertain = [entry for entry in decisions if entry["status"] == "uncertain"]
    rejected = [entry for entry in decisions if entry["status"] == "rejected"]

    lines = [
        f"# Structure Selection Summary: {summary['gene_name']}",
        "",
        f"- Input structures: {summary['input_structure_count']}",
        f"- Accepted: {summary['accepted_count']}",
        f"- Uncertain: {summary['uncertain_count']}",
        f"- Rejected: {summary['rejected_count']}",
        "",
        "## Reason Counts",
        "",
    ]

    for reason, count in sorted(summary["reason_counts"].items()):
        lines.append(f"- `{reason}`: {count}")

    lines.extend([
        "",
        "## Accepted",
        "",
    ])
    for entry in accepted:
        lines.append(
            f"- `{entry['structure_id']}`: {entry['ligand_state']}, "
            f"primary={entry['primary_ligand_code']}, ligands={entry['bound_ligands']}"
        )

    lines.extend([
        "",
        "## Uncertain",
        "",
    ])
    for entry in uncertain:
        lines.append(f"- `{entry['structure_id']}`: {', '.join(entry['reasons'])}")

    lines.extend([
        "",
        "## Rejected",
        "",
    ])
    for entry in rejected:
        lines.append(f"- `{entry['structure_id']}`: {', '.join(entry['reasons'])}")

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_candidate_records(gene_name: str, output_dir: str | Path) -> list[dict]:
    gene_dir = Path(output_dir) / gene_name
    record_paths = sorted(gene_dir.glob("receptor_structure_pdb_*.json"))
    if not record_paths:
        raise RuntimeError(
            f"No candidate receptor_structure JSON files found in {gene_dir}. "
            "Selection stage expects pre-fetched structure records."
        )

    records = []
    for record_path in record_paths:
        records.append(json.loads(record_path.read_text(encoding="utf-8")))
    return records


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Select usable holo receptor structures from pre-fetched receptor_structure JSON files."
        )
    )
    parser.add_argument("gene_name", help="Gene symbol, for example ACHE")
    parser.add_argument(
        "--output-dir",
        default=str(resolve_default_output_dir()),
        help=(
            "Base output directory. Input receptor_structure JSON files are read from "
            "output/<gene_name>/ and selection artifacts are written into "
            "output/<gene_name>/01_selection/."
        ),
    )
    args = parser.parse_args()

    try:
        gene_name = normalize_gene_name(args.gene_name)
        candidate_records = load_candidate_records(gene_name, args.output_dir)

        stage_dir = Path(args.output_dir) / gene_name / "01_selection"
        per_structure_dir = stage_dir / "per_structure"
        reports_dir = stage_dir / "reports"
        stage_dir.mkdir(parents=True, exist_ok=True)
        per_structure_dir.mkdir(parents=True, exist_ok=True)
        reports_dir.mkdir(parents=True, exist_ok=True)

        decisions = []
        for record in candidate_records:
            entry = build_selection_entry(record)
            decisions.append(entry)
            write_json(entry, per_structure_dir / f"{entry['structure_id']}.selection.json")

        summary = build_selection_summary(gene_name, decisions)

        write_json(decisions, stage_dir / "selection_decisions.json")
        write_json(summary, stage_dir / "index_selected_structures.json")
        write_selection_report(summary, decisions, reports_dir / "selection_summary.md")

        print(stage_dir / "index_selected_structures.json")
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"Failed to select structures for {args.gene_name}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
