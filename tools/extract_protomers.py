import argparse
import json
import sys
import urllib.error
from pathlib import Path

from rcsb_receptor_utils import ensure_mmcif_file
from rcsb_receptor_utils import normalize_gene_name
from rcsb_receptor_utils import resolve_default_output_dir
from rcsb_receptor_utils import write_json


def load_json(path: str | Path) -> dict | list:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def load_selection_decisions(gene_name: str, output_dir: str | Path) -> list[dict]:
    path = Path(output_dir) / gene_name / "01_selection" / "selection_decisions.json"
    if not path.exists():
        raise RuntimeError(
            f"Missing selection decisions file at {path}. "
            "Protomer extraction expects completed Stage 1 artifacts."
        )
    data = load_json(path)
    if not isinstance(data, list):
        raise RuntimeError(f"Selection decisions file at {path} must contain a JSON list")
    return data


def load_receptor_record(gene_name: str, structure_id: str, output_dir: str | Path) -> dict:
    path = Path(output_dir) / gene_name / f"receptor_structure_pdb_{structure_id}.json"
    if not path.exists():
        raise RuntimeError(f"Missing receptor_structure JSON for {structure_id} at {path}")
    data = load_json(path)
    if not isinstance(data, dict):
        raise RuntimeError(f"Receptor record at {path} must contain a JSON object")
    return data


def infer_protomer_mode(record: dict) -> tuple[str, str]:
    oligomeric_state = record.get("oligomeric_state")
    chains = record.get("chains") or []

    if oligomeric_state == "monomer" and len(chains) >= 1:
        if len(chains) == 1:
            return "single_chain_monomer", "Single chain present in receptor_structure record."
        return (
            "symmetry_equivalent_monomer_candidates",
            "Record is marked monomeric but contains multiple auth_asym_ids; first chain is used as the stage-2 canonical candidate.",
        )

    if oligomeric_state in {"dimer", "trimer", "tetramer", "oligomer"}:
        return (
            "oligomeric_chain_set",
            "Multiple chains may represent biologically relevant assembly context; downstream stages must verify chain equivalence and interface context.",
        )

    return (
        "unknown_chain_context",
        "Unable to infer robust protomer mode from receptor_structure record alone.",
    )


def build_protomer_record(
    record: dict,
    chain_id: str,
    chain_rank: int,
    protomer_mode: str,
    mmcif_relative_path: str,
) -> dict:
    structure_id = record.get("structure_id")
    gene_name = record.get("gene_name")
    bound_ligands = record.get("bound_ligands") or []

    if protomer_mode == "single_chain_monomer":
        stage_status = "accepted_reference_candidate"
        relevance = "primary_candidate"
    elif protomer_mode == "symmetry_equivalent_monomer_candidates":
        stage_status = "accepted_reference_candidate" if chain_rank == 0 else "symmetry_equivalent_candidate"
        relevance = "primary_candidate" if chain_rank == 0 else "symmetry_equivalent_candidate"
    elif protomer_mode == "oligomeric_chain_set":
        stage_status = "requires_interface_review"
        relevance = "ambiguous_oligomer_context"
    else:
        stage_status = "requires_review"
        relevance = "unknown"

    return {
        "protomer_instance_id": f"protomer_instance:{structure_id}:{chain_id}",
        "structure_id": structure_id,
        "gene_name": gene_name,
        "organism": record.get("organism"),
        "uniprot_id": record.get("uniprot_id"),
        "assembly_type": record.get("assembly_type"),
        "oligomeric_state": record.get("oligomeric_state"),
        "stoichiometry": record.get("stoichiometry"),
        "auth_asym_id": chain_id,
        "chain_rank_in_record": chain_rank,
        "protomer_mode": protomer_mode,
        "selection_status": "accepted",
        "stage_status": stage_status,
        "relevance": relevance,
        "primary_ligand_code": record.get("primary_ligand_code"),
        "bound_ligands": bound_ligands,
        "ligand_state": record.get("ligand_state"),
        "experimental_method": record.get("experimental_method"),
        "resolution_angstrom": record.get("resolution_angstrom"),
        "source_receptor_record": f"receptor_structure_pdb_{structure_id}.json",
        "source_mmcif_file": mmcif_relative_path,
        "notes": record.get("notes"),
    }


def build_ligand_instances(record: dict, protomer_records: list[dict]) -> list[dict]:
    ligand_instances = []
    structure_id = record.get("structure_id")
    bound_ligands = record.get("bound_ligands") or []

    if not protomer_records:
        return ligand_instances

    primary_protomer = protomer_records[0]
    for ligand_code in bound_ligands:
        ligand_instances.append(
            {
                "ligand_instance_id": f"ligand_instance:{structure_id}:{ligand_code}",
                "structure_id": structure_id,
                "ligand_code": ligand_code,
                "assigned_protomer_instance_id": primary_protomer["protomer_instance_id"],
                "assignment_mode": "stage2_primary_protomer_assignment",
                "assignment_status": "provisional",
                "is_primary_ligand": ligand_code == record.get("primary_ligand_code"),
                "ligand_state": record.get("ligand_state"),
                "source_mmcif_file": primary_protomer["source_mmcif_file"],
                "notes": record.get("notes"),
            }
        )
    return ligand_instances


def build_structure_stage_record(decision: dict, record: dict, mmcif_relative_path: str) -> dict:
    chains = record.get("chains") or []
    protomer_mode, protomer_mode_note = infer_protomer_mode(record)
    protomer_records = [
        build_protomer_record(record, chain_id, index, protomer_mode, mmcif_relative_path)
        for index, chain_id in enumerate(chains)
    ]
    ligand_instances = build_ligand_instances(record, protomer_records)

    return {
        "structure_id": record.get("structure_id"),
        "selection_status": decision.get("status"),
        "selection_reasons": decision.get("reasons"),
        "protomer_mode": protomer_mode,
        "protomer_mode_note": protomer_mode_note,
        "protomer_count": len(protomer_records),
        "ligand_instance_count": len(ligand_instances),
        "canonical_protomer_candidate_id": protomer_records[0]["protomer_instance_id"] if protomer_records else None,
        "status": "ready_for_alignment" if protomer_records else "requires_review",
    }


def write_protomer_report(gene_name: str, index_data: dict, output_path: str | Path) -> None:
    lines = [
        f"# Protomer Summary: {gene_name}",
        "",
        f"- Accepted structures from selection: {index_data['accepted_structure_count']}",
        f"- Extracted protomer instances: {index_data['accepted_protomer_count']}",
        f"- Extracted ligand instances: {index_data['accepted_ligand_instance_count']}",
        f"- Structures requiring review before alignment: {len(index_data['requires_review_structure_ids'])}",
        "",
        "## Requires Review",
        "",
    ]

    for structure_id in index_data["requires_review_structure_ids"]:
        lines.append(f"- `{structure_id}`")

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Extract protomer-stage artifacts from accepted structure-selection outputs and pre-fetched receptor_structure JSON files."
        )
    )
    parser.add_argument("gene_name", help="Gene symbol, for example ACHE")
    parser.add_argument(
        "--output-dir",
        default=str(resolve_default_output_dir()),
        help=(
            "Base output directory. Inputs are read from output/<gene_name>/ and "
            "output/<gene_name>/01_selection/, and protomer artifacts are written into "
            "output/<gene_name>/02_protomers/."
        ),
    )
    args = parser.parse_args()

    try:
        gene_name = normalize_gene_name(args.gene_name)
        decisions = load_selection_decisions(gene_name, args.output_dir)
        accepted_decisions = [entry for entry in decisions if entry.get("status") == "accepted"]
        if not accepted_decisions:
            raise RuntimeError(f"No accepted structures found in Stage 1 outputs for {gene_name}")

        stage_dir = Path(args.output_dir) / gene_name / "02_protomers"
        per_structure_dir = stage_dir / "per_structure"
        reports_dir = stage_dir / "reports"
        stage_dir.mkdir(parents=True, exist_ok=True)
        per_structure_dir.mkdir(parents=True, exist_ok=True)
        reports_dir.mkdir(parents=True, exist_ok=True)

        structure_stage_records = []
        all_protomers = []
        all_ligand_instances = []

        for decision in accepted_decisions:
            structure_id = decision["structure_id"]
            record = load_receptor_record(gene_name, structure_id, args.output_dir)
            chains = record.get("chains") or []
            protomer_mode, _protomer_mode_note = infer_protomer_mode(record)
            source_structures_dir = stage_dir / "source_structures"
            source_structures_dir.mkdir(parents=True, exist_ok=True)
            mmcif_path = ensure_mmcif_file(
                structure_id,
                source_structures_dir / f"{structure_id}.cif",
            )
            mmcif_relative_path = str(mmcif_path.relative_to(stage_dir.parent.parent)).replace("\\", "/")
            protomer_records = [
                build_protomer_record(record, chain_id, index, protomer_mode, mmcif_relative_path)
                for index, chain_id in enumerate(chains)
            ]
            ligand_instances = build_ligand_instances(record, protomer_records)
            structure_stage_record = build_structure_stage_record(decision, record, mmcif_relative_path)
            structure_stage_record["source_mmcif_file"] = mmcif_relative_path

            structure_dir = per_structure_dir / structure_id
            structure_dir.mkdir(parents=True, exist_ok=True)
            for protomer_record in protomer_records:
                write_json(
                    protomer_record,
                    structure_dir / f"protomer_{protomer_record['auth_asym_id']}.json",
                )
            write_json(ligand_instances, structure_dir / "ligand_instances.json")
            write_json(structure_stage_record, structure_dir / "protomer_stage.json")

            structure_stage_records.append(structure_stage_record)
            all_protomers.extend(protomer_records)
            all_ligand_instances.extend(ligand_instances)

        reference_candidates = [
            protomer
            for protomer in all_protomers
            if protomer["stage_status"] == "accepted_reference_candidate"
        ]
        reference_candidates.sort(
            key=lambda item: (
                item.get("resolution_angstrom") is None,
                item.get("resolution_angstrom") or 999.0,
                item.get("structure_id") or "",
                item.get("auth_asym_id") or "",
            )
        )
        canonical_reference_candidate = reference_candidates[0] if reference_candidates else None

        requires_review_structure_ids = [
            item["structure_id"]
            for item in structure_stage_records
            if item["status"] != "ready_for_alignment"
        ]

        index_data = {
            "gene_name": gene_name,
            "pipeline_stage": 2,
            "stage_name": "protomer_extraction",
            "input_selection_file": f"output/{gene_name}/01_selection/selection_decisions.json",
            "accepted_structure_count": len(accepted_decisions),
            "accepted_structure_ids": [item["structure_id"] for item in accepted_decisions],
            "accepted_protomer_count": len(all_protomers),
            "accepted_ligand_instance_count": len(all_ligand_instances),
            "canonical_reference_candidate_id": canonical_reference_candidate["protomer_instance_id"] if canonical_reference_candidate else None,
            "requires_review_structure_ids": requires_review_structure_ids,
        }

        write_json(index_data, stage_dir / "index_protomers.json")
        write_json(structure_stage_records, stage_dir / "structure_stage_records.json")
        write_json(all_ligand_instances, stage_dir / "accepted_ligand_instances.json")
        if canonical_reference_candidate is not None:
            write_json(canonical_reference_candidate, stage_dir / "canonical_reference_candidate.json")
        write_protomer_report(gene_name, index_data, reports_dir / "protomer_summary.md")

        print(stage_dir / "index_protomers.json")
        return 0
    except urllib.error.HTTPError as exc:
        print(f"HTTP error while extracting protomers for {args.gene_name}: {exc}", file=sys.stderr)
        return 1
    except urllib.error.URLError as exc:
        print(f"Network error while extracting protomers for {args.gene_name}: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:  # noqa: BLE001
        print(f"Failed to extract protomers for {args.gene_name}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
