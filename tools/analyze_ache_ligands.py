import json

from rcsb_receptor_utils import fetch_receptor_structure_record


SAMPLE_IDS = [
    "1B41",
    "1F8U",
    "1N5M",
    "4EY4",
    "4EY5",
    "4EY6",
    "4EY7",
    "4EY8",
    "5DTI",
    "6CQW",
    "6CQX",
    "6CQY",
]


def main() -> None:
    rows = []
    for structure_id in SAMPLE_IDS:
        try:
            record = fetch_receptor_structure_record(structure_id)
            rows.append(
                {
                    "structure_id": structure_id,
                    "title": record.get("structure_title"),
                    "ligand_state": record.get("ligand_state"),
                    "bound_ligands": record.get("bound_ligands"),
                    "active_site_status": record.get("active_site_status"),
                }
            )
        except Exception as exc:  # noqa: BLE001
            rows.append(
                {
                    "structure_id": structure_id,
                    "error": str(exc),
                }
            )

    print(json.dumps(rows, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
