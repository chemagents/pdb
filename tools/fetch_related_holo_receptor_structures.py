import argparse
import sys
import urllib.error
from pathlib import Path

from rcsb_receptor_utils import TARGET_ORGANISMS
from rcsb_receptor_utils import detect_query_kind
from rcsb_receptor_utils import extract_identity_from_structure
from rcsb_receptor_utils import fetch_receptor_structure_record
from rcsb_receptor_utils import is_holo_record
from rcsb_receptor_utils import normalize_gene_name
from rcsb_receptor_utils import normalize_structure_id
from rcsb_receptor_utils import resolve_default_output_dir
from rcsb_receptor_utils import resolve_gene_scoped_output_dir
from rcsb_receptor_utils import search_related_structure_ids
from rcsb_receptor_utils import write_record_json


def resolve_candidate_structure_ids(query: str) -> tuple[str, str | None, list[str]]:
    query_kind = detect_query_kind(query)
    if query_kind == "pdb_id":
        structure_id = normalize_structure_id(query)
        uniprot_id, gene_name, _protein_name = extract_identity_from_structure(structure_id)
        candidate_ids = search_related_structure_ids(uniprot_id, gene_name)
        normalized_gene_name = normalize_gene_name(gene_name) if gene_name else None
        return query_kind, normalized_gene_name, candidate_ids

    gene_name = normalize_gene_name(query)
    return query_kind, gene_name, search_related_structure_ids(None, gene_name)


def resolve_output_dir(base_output_dir: str | Path, query_kind: str, gene_name: str | None) -> Path:
    if query_kind == "gene_name" and gene_name:
        return resolve_gene_scoped_output_dir(base_output_dir, gene_name)
    return Path(base_output_dir)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Find human and mouse holo structures for a receptor using either a PDB structure ID "
            "or a gene name, then save normalized receptor_structure JSON files."
        )
    )
    parser.add_argument("query", help="Reference PDB code like 5U09 or gene symbol like ACHE")
    parser.add_argument(
        "--output-dir",
        default=str(resolve_default_output_dir()),
        help=(
            "Base directory where normalized JSON files will be written. "
            "For gene-name queries, files are saved into a gene-specific subdirectory."
        ),
    )
    args = parser.parse_args()

    try:
        query_kind, gene_name, candidate_ids = resolve_candidate_structure_ids(args.query)
        output_dir = resolve_output_dir(args.output_dir, query_kind, gene_name)
        output_dir.mkdir(parents=True, exist_ok=True)

        saved_files = []
        for candidate_id in candidate_ids:
            record = fetch_receptor_structure_record(candidate_id)
            organism = record.get("organism")
            if organism not in TARGET_ORGANISMS or not is_holo_record(record):
                continue

            output_path = output_dir / f"receptor_structure_pdb_{candidate_id}.json"
            write_record_json(record, output_path)
            saved_files.append(output_path)

        for path in saved_files:
            print(path)
        return 0
    except urllib.error.HTTPError as exc:
        print(f"HTTP error while processing {args.query}: {exc}", file=sys.stderr)
        return 1
    except urllib.error.URLError as exc:
        print(f"Network error while processing {args.query}: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:  # noqa: BLE001
        print(f"Failed to fetch related holo structures for {args.query}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
