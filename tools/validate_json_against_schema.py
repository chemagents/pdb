import argparse
import json
import re
import sys
from pathlib import Path


SCHEMA_DIR = Path(__file__).resolve().parent.parent / "schemas"
ENTITY_TO_SCHEMA_PATH = {
    "receptor_structure": SCHEMA_DIR / "receptor_structure.schema.en.yaml",
    "binding_site": SCHEMA_DIR / "binding_site.schema.en.yaml",
    "ligand_complex": SCHEMA_DIR / "ligand_complex.schema.en.yaml",
    "ligand": SCHEMA_DIR / "ligand.schema.en.yaml",
}


def parse_scalar(value: str):
    text = value.strip()
    if not text:
        return ""
    if text == "null":
        return None
    if text == "true":
        return True
    if text == "false":
        return False
    if (text.startswith('"') and text.endswith('"')) or (
        text.startswith("'") and text.endswith("'")
    ):
        return text[1:-1]
    if re.fullmatch(r"-?\d+", text):
        return int(text)
    if re.fullmatch(r"-?\d+\.\d+", text):
        return float(text)
    if text.startswith("[") and text.endswith("]"):
        try:
            return json.loads(text)
        except json.JSONDecodeError:
            return text
    return text


def load_schema_fields(schema_path: str | Path) -> list[dict]:
    lines = Path(schema_path).read_text(encoding="utf-8").splitlines()
    fields = []
    current_field = None
    in_fields = False
    in_allowed_values = False
    in_item_schema = False

    for line in lines:
        stripped = line.strip()

        if stripped == "fields:":
            in_fields = True
            continue

        if in_fields and line.startswith("interpretation_rules:"):
            break

        if not in_fields:
            continue

        if line.startswith("  - name:"):
            if current_field:
                fields.append(current_field)
            current_field = {
                "name": parse_scalar(line.split(":", 1)[1]),
                "type": None,
                "required": False,
                "allowed_values": None,
                "constraints": {},
                "item_schema": None,
            }
            in_allowed_values = False
            in_item_schema = False
            continue

        if current_field is None:
            continue

        if line.startswith("    allowed_values:"):
            current_field["allowed_values"] = []
            in_allowed_values = True
            in_item_schema = False
            continue

        if in_allowed_values:
            if line.startswith("      - "):
                current_field["allowed_values"].append(parse_scalar(line.split("-", 1)[1]))
                continue
            in_allowed_values = False

        if line.startswith("    item_schema:"):
            current_field["item_schema"] = {}
            in_item_schema = True
            continue

        if line.startswith("    type:"):
            current_field["type"] = parse_scalar(line.split(":", 1)[1])
            in_item_schema = False
        elif line.startswith("    required:"):
            current_field["required"] = bool(parse_scalar(line.split(":", 1)[1]))
            in_item_schema = False
        elif line.startswith("    constraints:"):
            current_field["constraints"] = {}
            in_item_schema = False
        elif line.startswith("      regex:"):
            current_field.setdefault("constraints", {})["regex"] = parse_scalar(
                line.split(":", 1)[1]
            )
        elif in_item_schema and line.startswith("      "):
            item_text = line.strip()
            if ":" in item_text:
                key, value = item_text.split(":", 1)
                current_field["item_schema"][key.strip()] = parse_scalar(value)
        elif line.startswith("    "):
            in_item_schema = False

    if current_field:
        fields.append(current_field)

    return fields


def validate_scalar_type(value, expected_type: str) -> bool:
    if expected_type == "string":
        return isinstance(value, str)
    if expected_type == "integer":
        return isinstance(value, int) and not isinstance(value, bool)
    if expected_type == "number":
        return (isinstance(value, int) or isinstance(value, float)) and not isinstance(value, bool)
    if expected_type == "boolean":
        return isinstance(value, bool)
    return True


def validate_field(field: dict, value, path_prefix: str) -> list[str]:
    errors = []
    field_name = field["name"]
    field_path = f"{path_prefix}{field_name}"

    if value is None:
        return errors

    field_type = field.get("type")
    if field_type in {"string", "integer", "number", "boolean"}:
        if not validate_scalar_type(value, field_type):
            errors.append(f"{field_path}: expected {field_type}, got {type(value).__name__}")
            return errors
    elif field_type == "array[string]":
        if not isinstance(value, list):
            errors.append(f"{field_path}: expected array[string], got {type(value).__name__}")
            return errors
        for index, item in enumerate(value):
            if not isinstance(item, str):
                errors.append(
                    f"{field_path}[{index}]: expected string, got {type(item).__name__}"
                )
    elif field_type == "array[object]":
        if not isinstance(value, list):
            errors.append(f"{field_path}: expected array[object], got {type(value).__name__}")
            return errors
        item_schema = field.get("item_schema") or {}
        for index, item in enumerate(value):
            if not isinstance(item, dict):
                errors.append(
                    f"{field_path}[{index}]: expected object, got {type(item).__name__}"
                )
                continue
            for required_key, expected_item_type in item_schema.items():
                if required_key not in item:
                    errors.append(f"{field_path}[{index}].{required_key}: missing required key")
                    continue
                if not validate_scalar_type(item[required_key], expected_item_type):
                    errors.append(
                        f"{field_path}[{index}].{required_key}: expected {expected_item_type}, "
                        f"got {type(item[required_key]).__name__}"
                    )

    allowed_values = field.get("allowed_values")
    if allowed_values and value not in allowed_values:
        errors.append(f"{field_path}: value {value!r} not in allowed_values")

    regex_pattern = (field.get("constraints") or {}).get("regex")
    if regex_pattern and isinstance(value, str) and not re.fullmatch(regex_pattern, value):
        errors.append(f"{field_path}: value {value!r} does not match regex {regex_pattern!r}")

    return errors


def validate_record(record: dict, schema_fields: list[dict]) -> list[str]:
    errors = []
    known_field_names = {field["name"] for field in schema_fields}

    for field in schema_fields:
        field_name = field["name"]
        if field.get("required") and field_name not in record:
            errors.append(f"{field_name}: missing required field")
            continue
        if field_name in record:
            errors.extend(validate_field(field, record[field_name], ""))

    for field_name in record:
        if field_name not in known_field_names:
            errors.append(f"{field_name}: unexpected field not declared in schema")

    return errors


def infer_entity_name(json_path: str | Path, data: dict) -> str:
    path = Path(json_path)
    if path.name.startswith("receptor_structure_"):
        return "receptor_structure"
    if "site_id" in data:
        return "binding_site"
    if "complex_id" in data:
        return "ligand_complex"
    if "ligand_entity_id" in data:
        return "ligand"
    raise ValueError(
        "Unable to infer entity type. Pass --entity receptor_structure|binding_site|ligand_complex|ligand."
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate a normalized JSON record against the corresponding YAML schema."
    )
    parser.add_argument("json_path", help="Path to normalized JSON record")
    parser.add_argument(
        "--entity",
        choices=sorted(ENTITY_TO_SCHEMA_PATH.keys()),
        help="Explicit entity type if it cannot be inferred from the JSON record",
    )
    args = parser.parse_args()

    try:
        json_path = Path(args.json_path)
        data = json.loads(json_path.read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            raise RuntimeError("Top-level JSON value must be an object record")

        entity_name = args.entity or infer_entity_name(json_path, data)
        schema_fields = load_schema_fields(ENTITY_TO_SCHEMA_PATH[entity_name])
        errors = validate_record(data, schema_fields)
        if errors:
            for error in errors:
                print(error, file=sys.stderr)
            return 1

        print(f"VALID {entity_name} {json_path}")
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"Validation failed for {args.json_path}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
