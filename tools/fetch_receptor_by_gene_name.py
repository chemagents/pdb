import argparse
import sys
import urllib.error
from pathlib import Path

from rcsb_receptor_utils import fetch_receptor_structure_record
from rcsb_receptor_utils import normalize_gene_name
from rcsb_receptor_utils import resolve_default_output_dir
from rcsb_receptor_utils import resolve_gene_scoped_output_dir
from rcsb_receptor_utils import search_structure_ids_by_gene_name
from rcsb_receptor_utils import write_record_json


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Find receptor structures by gene name in RCSB PDB and save normalized receptor_structure JSON files."
        )
    )
    parser.add_argument("gene_name", help="Gene symbol, for example ACHE")
    parser.add_argument(
        "--output-dir",
        default=str(resolve_default_output_dir()),
        help=(
            "Base directory where receptor_structure_pdb_<CODE>.json files will be written. "
            "Files are saved into a gene-specific subdirectory."
        ),
    )
    args = parser.parse_args()

    try:
        gene_name = normalize_gene_name(args.gene_name)
        structure_ids = search_structure_ids_by_gene_name(gene_name)
        if not structure_ids:
            raise RuntimeError(f"No structures found for gene_name {gene_name}")

        output_dir = resolve_gene_scoped_output_dir(args.output_dir, gene_name)
        output_dir.mkdir(parents=True, exist_ok=True)

        saved_files = []
        for structure_id in structure_ids:
            record = fetch_receptor_structure_record(structure_id)
            output_path = output_dir / f"receptor_structure_pdb_{structure_id}.json"
            write_record_json(record, output_path)
            saved_files.append(output_path)

        for saved_file in saved_files:
            print(saved_file)
        return 0
    except urllib.error.HTTPError as exc:
        print(f"HTTP error while fetching {args.gene_name}: {exc}", file=sys.stderr)
        return 1
    except urllib.error.URLError as exc:
        print(f"Network error while fetching {args.gene_name}: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:  # noqa: BLE001
        print(f"Failed to fetch structures for {args.gene_name}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
