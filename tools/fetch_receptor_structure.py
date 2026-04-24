import argparse
import sys
import urllib.error
from pathlib import Path

from rcsb_receptor_utils import fetch_receptor_structure_record
from rcsb_receptor_utils import normalize_structure_id
from rcsb_receptor_utils import resolve_default_output_dir
from rcsb_receptor_utils import resolve_record_output_dir
from rcsb_receptor_utils import write_record_json


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Fetch receptor structure metadata from RCSB PDB and save normalized JSON."
    )
    parser.add_argument("structure_id", help="PDB structure code, for example 4EY7")
    parser.add_argument(
        "--output-dir",
        default=str(resolve_default_output_dir()),
        help="Directory where receptor_structure_pdb_<CODE>.json will be written",
    )
    args = parser.parse_args()

    try:
        structure_id = normalize_structure_id(args.structure_id)
        record = fetch_receptor_structure_record(structure_id)
    except urllib.error.HTTPError as exc:
        print(f"HTTP error while fetching {args.structure_id}: {exc}", file=sys.stderr)
        return 1
    except urllib.error.URLError as exc:
        print(f"Network error while fetching {args.structure_id}: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:  # noqa: BLE001
        print(f"Failed to build record for {args.structure_id}: {exc}", file=sys.stderr)
        return 1

    output_dir = resolve_record_output_dir(args.output_dir, record)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"receptor_structure_pdb_{structure_id}.json"
    write_record_json(record, output_path)
    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
