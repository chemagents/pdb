import json
import re
import urllib.parse
import urllib.request
from pathlib import Path


RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{structure_id}"
RCSB_POLYMER_ENTITY_URL = (
    "https://data.rcsb.org/rest/v1/core/polymer_entity/{structure_id}/{entity_id}"
)
RCSB_NONPOLYMER_ENTITY_URL = (
    "https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{structure_id}/{entity_id}"
)
RCSB_MMCIF_DOWNLOAD_URL = "https://files.rcsb.org/download/{structure_id}.cif"
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query?json={encoded_query}"

USER_AGENT = "receptor-structure-fetcher/2.1"
REQUEST_TIMEOUT = 60

TARGET_ORGANISMS = {"Homo sapiens", "Mus musculus"}
SOLVENT_CODES = {"HOH", "WAT"}
NON_FUNCTIONAL_COMPONENT_CODES = {
    "HOH",
    "WAT",
    "DOD",
    "GOL",
    "EDO",
    "PEG",
    "PG4",
    "MPD",
    "MES",
    "HEP",
    "HEZ",
    "TRS",
    "ACT",
    "ACY",
    "FMT",
    "DMS",
    "SO4",
    "PO4",
    "NO3",
    "CO3",
    "CL",
    "NA",
    "K",
    "CA",
    "MG",
    "ZN",
    "IOD",
    "BR",
    "NAG",
    "BMA",
    "MAN",
    "FUC",
    "GAL",
    "GLC",
    "NDG",
    "SIA",
}

GLYCAN_COMPONENT_CODES = {"NAG", "BMA", "MAN", "FUC", "GAL", "GLC", "NDG", "SIA"}
PEG_OR_POLYOL_COMPONENT_CODES = {"PEG", "PG4", "GOL", "EDO", "MPD"}
PEPTIDE_KEYWORDS = {"PEPTIDE", "TOXIN", "FASCICULIN"}
TITLE_APO_PATTERNS = (" apo ", "apo state", "apo-form", "apo form")
TITLE_COMPLEX_PATTERNS = ("complex with", "in complex with", "bound to", "inhibited by")
TITLE_PEPTIDE_PATTERNS = ("fasciculin", "toxin", "peptide")
TITLE_REACTIVATOR_PATTERNS = ("hi-6", "obidoxime", "pralidoxime", "oxime")
TITLE_COVALENT_PATTERNS = ("sarin", "vx", "tabun", "soman", "organophosphate", "phosphonyl")
REVERSIBLE_LIGAND_NAMES = ("donepezil", "huperzine", "galantamine", "gallamine")
REACTIVATOR_NAMES = ("hi-6", "obidoxime", "pralidoxime")
ORGANOPHOSPHATE_NAMES = ("vx", "sarin", "soman", "tabun")
WEAK_SMALL_MOLECULE_CODES = {"PE8", "P6G"}
NAME_TOKEN_RE = re.compile(r"[a-z0-9]+")
STOPWORD_TOKENS = {
    "a",
    "acetylcholinesterase",
    "and",
    "bound",
    "complex",
    "crystal",
    "form",
    "human",
    "in",
    "inhibited",
    "mouse",
    "of",
    "recombinant",
    "state",
    "structure",
    "the",
    "with",
}
TITLE_NAME_HINTS = {
    "donepezil": {"donepezil"},
    "huperzine": {"huperzine"},
    "galantamine": {"galantamine"},
    "gallamine": {"gallamine", "triethylethanaminium", "benzene", "tris"},
    "fasciculin": {"fasciculin"},
    "hi-6": {"hi", "6", "oxime"},
}

PDB_ID_RE = re.compile(r"^[0-9][A-Z0-9]{3}$")
GENE_NAME_RE = re.compile(r"^[A-Z0-9][A-Z0-9_-]{1,31}$")


def fetch_json(url: str) -> dict:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": USER_AGENT,
            "Accept": "application/json",
        },
    )
    with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT) as response:
        return json.load(response)


def fetch_text(url: str) -> str:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": USER_AGENT,
        },
    )
    with urllib.request.urlopen(req, timeout=REQUEST_TIMEOUT) as response:
        return response.read().decode("utf-8")


def search_rcsb(payload: dict) -> dict:
    encoded_query = urllib.parse.quote(json.dumps(payload, separators=(",", ":")))
    url = RCSB_SEARCH_URL.format(encoded_query=encoded_query)
    return fetch_json(url)


def first_or_none(items):
    if isinstance(items, list) and items:
        return items[0]
    return None


def normalize_structure_id(structure_id: str) -> str:
    code = structure_id.strip().upper()
    if not code:
        raise ValueError("structure_id is empty")
    return code


def normalize_gene_name(gene_name: str) -> str:
    value = gene_name.strip().upper()
    if not value:
        raise ValueError("gene_name is empty")
    return value


def detect_query_kind(value: str) -> str:
    normalized = value.strip().upper()
    if not normalized:
        raise ValueError("query is empty")
    if PDB_ID_RE.fullmatch(normalized):
        return "pdb_id"
    if GENE_NAME_RE.fullmatch(normalized):
        return "gene_name"
    return "gene_name"


def resolve_gene_scoped_output_dir(output_dir: str | Path, gene_name: str) -> Path:
    base_dir = Path(output_dir)
    normalized_gene_name = normalize_gene_name(gene_name)
    return base_dir / normalized_gene_name


def resolve_default_output_dir() -> Path:
    return Path(__file__).resolve().parent.parent / "output"


def resolve_record_output_dir(base_output_dir: str | Path, record: dict) -> Path:
    gene_name = record.get("gene_name")
    if isinstance(gene_name, str) and gene_name.strip():
        return resolve_gene_scoped_output_dir(base_output_dir, gene_name)
    return Path(base_output_dir)


def build_external_ids(structure_id: str, uniprot_id: str | None) -> list[dict]:
    external_ids = [
        {
            "registry": "PDB",
            "id": structure_id,
            "relation": "primary_entry",
        }
    ]
    if uniprot_id:
        external_ids.append(
            {
                "registry": "UniProt",
                "id": uniprot_id,
                "relation": "sequence_reference",
            }
        )
    return external_ids


def build_related_entities(gene_name: str | None) -> list[dict]:
    if not gene_name:
        return []
    return [
        {
            "entity_type": "receptor_structure",
            "entity_id": f"receptor_structure_gene:{normalize_gene_name(gene_name)}",
            "relation_type": "has_state_record",
            "note": "Grouped with other receptor_structure records for the same gene symbol.",
        }
    ]


def get_primary_polymer_entity_id(entry: dict) -> str:
    polymer_entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get(
        "polymer_entity_ids", []
    )
    if not polymer_entity_ids:
        raise RuntimeError("No polymer entities found")
    return str(polymer_entity_ids[0])


def get_uniprot_id(polymer_entity: dict) -> str | None:
    identifiers = polymer_entity.get(
        "rcsb_polymer_entity_container_identifiers", {}
    ).get("reference_sequence_identifiers", [])
    for item in identifiers:
        if item.get("database_name") == "UniProt":
            return item.get("database_accession")
    return None


def get_gene_name(polymer_entity: dict) -> str | None:
    for source in polymer_entity.get("rcsb_entity_source_organism", []):
        gene_names = source.get("rcsb_gene_name", [])
        for gene_name in gene_names:
            value = gene_name.get("value")
            if value:
                return value
    return None


def get_organism(polymer_entity: dict) -> tuple[str | None, int | None]:
    sources = polymer_entity.get("rcsb_entity_source_organism", [])
    if not sources:
        return None, None
    source = sources[0]
    return source.get("ncbi_scientific_name"), source.get("ncbi_taxonomy_id")


def get_mutations(polymer_entity: dict) -> list[str]:
    mutations = polymer_entity.get("entity_src_gen", [])
    result = []
    for item in mutations:
        value = item.get("pdbx_gene_src_mutation")
        if value and value not in result:
            result.append(value)
    return result


def get_sequence_range(polymer_entity: dict) -> str | None:
    identifiers = polymer_entity.get("rcsb_polymer_entity_container_identifiers", {})
    ranges = identifiers.get("auth_to_entity_poly_seq_mapping", [])
    starts = []
    ends = []
    for item in ranges:
        beg_seq_id = item.get("beg_seq_id")
        end_seq_id = item.get("end_seq_id")
        if isinstance(beg_seq_id, int) and isinstance(end_seq_id, int):
            starts.append(beg_seq_id)
            ends.append(end_seq_id)
    if starts and ends:
        return f"{min(starts)}-{max(ends)}"
    return None


def fetch_nonpolymer_candidates(structure_id: str, entry: dict) -> list[dict]:
    candidates = []
    entry_ids = entry.get("rcsb_entry_container_identifiers", {})
    for entity_id in entry_ids.get("non_polymer_entity_ids", []):
        entity = fetch_json(
            RCSB_NONPOLYMER_ENTITY_URL.format(
                structure_id=structure_id,
                entity_id=entity_id,
            )
        )
        component_id = entity.get("pdbx_entity_nonpoly", {}).get("comp_id")
        if not component_id:
            component_id = entity.get(
                "rcsb_nonpolymer_entity_container_identifiers", {}
            ).get("nonpolymer_comp_id")
        if component_id:
            candidates.append(
                {
                    "entity_id": str(entity_id),
                    "comp_id": str(component_id).upper(),
                    "name": (
                        entity.get("pdbx_entity_nonpoly", {}).get("name")
                        or entity.get("rcsb_nonpolymer_entity", {}).get("pdbx_description")
                        or ""
                    ).upper(),
                    "formula_weight": entity.get("rcsb_nonpolymer_entity", {}).get("formula_weight"),
                    "instance_count": entity.get("rcsb_nonpolymer_entity", {}).get("pdbx_number_of_molecules"),
                    "classification": None,
                    "signals": [],
                    "score": 0,
                }
            )
    return candidates


def classify_candidate(candidate: dict) -> str:
    code = candidate.get("comp_id")
    name = candidate.get("name", "")
    if code in {"HOH", "WAT", "DOD"}:
        return "water"
    if code in GLYCAN_COMPONENT_CODES:
        return "glycan"
    if code in PEG_OR_POLYOL_COMPONENT_CODES or str(code).startswith("PEG"):
        return "peg_or_polyol_fragment"
    if code in NON_FUNCTIONAL_COMPONENT_CODES:
        return "buffer_or_additive"
    if code in WEAK_SMALL_MOLECULE_CODES:
        return "unknown"
    if any(keyword in name for keyword in PEPTIDE_KEYWORDS):
        return "peptide_or_protein_binder"
    return "small_molecule_candidate"


def score_candidate(candidate: dict, title_text: str, apo_support: bool) -> None:
    code = candidate["comp_id"]
    name = candidate.get("name", "").lower()
    title = f" {title_text.lower()} "
    classification = candidate["classification"]

    if classification in {"water", "buffer_or_additive", "glycan", "peg_or_polyol_fragment"}:
        candidate["score"] -= 5
        candidate["signals"].append(f"excluded:{classification}")
        return

    if classification == "small_molecule_candidate":
        candidate["score"] += 4
        candidate["signals"].append("small_molecule_candidate")
    elif classification == "peptide_or_protein_binder":
        candidate["score"] += 3
        candidate["signals"].append("peptide_or_protein_binder")
    elif classification == "unknown":
        candidate["score"] += 1
        candidate["signals"].append("unknown_nonpolymer_candidate")

    if f" {code.lower()} " in title:
        candidate["score"] += 5
        candidate["signals"].append("title_mentions_code")

    normalized_name = re.sub(r"[^a-z0-9]+", " ", name).strip()
    if normalized_name and normalized_name in title:
        candidate["score"] += 5
        candidate["signals"].append("title_mentions_name")

    title_tokens = {
        token
        for token in NAME_TOKEN_RE.findall(title)
        if token not in STOPWORD_TOKENS and len(token) > 2
    }
    name_tokens = {
        token
        for token in NAME_TOKEN_RE.findall(name)
        if token not in STOPWORD_TOKENS and len(token) > 2
    }
    overlap = title_tokens & name_tokens
    if overlap:
        candidate["score"] += min(4, len(overlap) * 2)
        candidate["signals"].append(f"title_name_token_overlap:{','.join(sorted(overlap))}")

    for title_hint, hint_tokens in TITLE_NAME_HINTS.items():
        if title_hint in title and name_tokens & hint_tokens:
            candidate["score"] += 4
            candidate["signals"].append(f"title_hint_match:{title_hint}")

    if any(pattern in title for pattern in TITLE_COMPLEX_PATTERNS):
        candidate["score"] += 2
        candidate["signals"].append("complex_title_context")

    for ligand_name in REVERSIBLE_LIGAND_NAMES + REACTIVATOR_NAMES + ORGANOPHOSPHATE_NAMES:
        if ligand_name in title and ligand_name in name:
            candidate["score"] += 4
            candidate["signals"].append("title_matches_known_ligand_name")

    if apo_support and classification == "small_molecule_candidate":
        candidate["score"] -= 2
        candidate["signals"].append("apo_title_conflict")


def infer_inhibitor_class(title_text: str, meaningful: list[dict], ligand_state: str) -> str | None:
    title = title_text.lower()
    if ligand_state == "apo":
        return "none"
    if any(token in title for token in ORGANOPHOSPHATE_NAMES):
        return "organophosphate"
    if "carbamate" in title:
        return "carbamate"
    if ligand_state == "reactivator_complex" and any(token in title for token in REACTIVATOR_NAMES):
        return "reactivator"
    if any(token in title for token in REVERSIBLE_LIGAND_NAMES):
        return "reversible_inhibitor"
    if any(candidate["comp_id"] in {"HI6"} for candidate in meaningful):
        return "reactivator"
    return "unknown"


def summarize_ligand_context(title_text: str, candidates: list[dict]) -> dict:
    title = f" {title_text.lower()} "
    apo_support = any(pattern in title for pattern in TITLE_APO_PATTERNS)
    peptide_support = any(pattern in title for pattern in TITLE_PEPTIDE_PATTERNS)
    reactivator_support = any(pattern in title for pattern in TITLE_REACTIVATOR_PATTERNS)
    covalent_support = any(pattern in title for pattern in TITLE_COVALENT_PATTERNS)

    for candidate in candidates:
        candidate["classification"] = classify_candidate(candidate)
        score_candidate(candidate, title_text, apo_support)

    meaningful = [candidate for candidate in candidates if candidate["score"] >= 5]
    weak = [candidate for candidate in candidates if 2 <= candidate["score"] < 5]
    small_molecules = [candidate for candidate in meaningful if candidate["classification"] == "small_molecule_candidate"]
    peptides = [candidate for candidate in meaningful if candidate["classification"] == "peptide_or_protein_binder"]

    bound_ligands = []
    for candidate in meaningful:
        code = candidate["comp_id"]
        if code not in bound_ligands:
            bound_ligands.append(code)

    primary_ligand_code = None
    if len(small_molecules) == 1:
        primary_ligand_code = small_molecules[0]["comp_id"]
    elif len(meaningful) == 1:
        primary_ligand_code = meaningful[0]["comp_id"]

    notes = []
    if covalent_support and meaningful:
        ligand_state = "covalent_inhibitor_complex"
        active_site_status = "covalently_modified"
    elif len(meaningful) >= 2 and (covalent_support or reactivator_support):
        ligand_state = "mixed_complex"
        active_site_status = "covalently_modified" if covalent_support else "occupied"
    elif len(small_molecules) == 1:
        code = small_molecules[0]["comp_id"]
        if code in {"VX", "SAR", "TAB", "SOM"}:
            ligand_state = "covalent_inhibitor_complex"
            active_site_status = "covalently_modified"
        elif code in {"HI6"}:
            ligand_state = "reactivator_complex"
            active_site_status = "occupied"
        else:
            ligand_state = "noncovalent_complex"
            active_site_status = "occupied"
    elif peptide_support:
        ligand_state = "holo"
        active_site_status = "blocked"
        bound_ligands = []
        primary_ligand_code = None
        notes.append("Peptide/protein-bound inhibitory context outside current small-molecule prioritization logic.")
    elif not meaningful and not weak:
        ligand_state = "apo"
        active_site_status = "free"
    else:
        ligand_state = "unknown"
        active_site_status = "unknown"
        notes.append("Uncertain ligand interpretation; weak nonpolymer candidates present.")

    inhibitor_class = infer_inhibitor_class(title_text, meaningful, ligand_state)
    if ligand_state == "unknown" and apo_support and not meaningful:
        ligand_state = "apo"
        active_site_status = "free"
        inhibitor_class = "none"
        notes = [note for note in notes if note != "Uncertain ligand interpretation; weak nonpolymer candidates present."]

    excluded = [candidate["comp_id"] for candidate in candidates if candidate["comp_id"] not in bound_ligands]
    if apo_support and ligand_state not in {"apo", "unknown"}:
        notes.append("Title suggests apo-like context but meaningful ligand evidence was retained.")
    if excluded:
        notes.append(f"Excluded non-holo-defining components: {', '.join(excluded)}.")
    if primary_ligand_code is None and len(meaningful) > 1:
        notes.append("Multiple meaningful ligands; no single primary ligand selected.")

    return {
        "bound_ligands": bound_ligands,
        "primary_ligand_code": primary_ligand_code,
        "ligand_state": ligand_state,
        "active_site_status": active_site_status,
        "inhibitor_class": inhibitor_class,
        "notes_suffix": " ".join(notes).strip() or None,
    }


def infer_oligomeric_state(entry: dict) -> tuple[str, str | None, str | None]:
    assembly_count = entry.get("rcsb_entry_info", {}).get("assembly_count")
    if assembly_count == 1:
        return "monomer", "deposited_assembly", "A"
    if assembly_count == 2:
        return "dimer", "deposited_assembly", "A2"
    if assembly_count == 3:
        return "trimer", "deposited_assembly", "A3"
    if assembly_count == 4:
        return "tetramer", "deposited_assembly", "A4"
    if isinstance(assembly_count, int) and assembly_count > 4:
        return "oligomer", "deposited_assembly", None
    return "unknown", None, None


def infer_construct_type(polymer_entity: dict) -> str:
    description = (polymer_entity.get("entity", {}) or {}).get("pdbx_description", "")
    text = description.lower()
    if "catalytic domain" in text:
        return "catalytic_domain"
    if "extracellular" in text:
        return "extracellular_domain"
    if "fragment" in text or "truncated" in text:
        return "truncated_construct"
    if "fusion" in text:
        return "fusion_construct"
    return "unknown"


def fetch_entry_and_polymer_entity(structure_id: str) -> tuple[dict, dict]:
    normalized_structure_id = normalize_structure_id(structure_id)
    entry = fetch_json(RCSB_ENTRY_URL.format(structure_id=normalized_structure_id))
    polymer_entity_id = get_primary_polymer_entity_id(entry)
    polymer_entity = fetch_json(
        RCSB_POLYMER_ENTITY_URL.format(
            structure_id=normalized_structure_id,
            entity_id=polymer_entity_id,
        )
    )
    return entry, polymer_entity


def build_receptor_structure_record(structure_id: str, entry: dict, polymer_entity: dict) -> dict:
    entry_info = entry.get("rcsb_entry_info", {})
    entity = polymer_entity.get("entity", {})
    citation = entry.get("rcsb_primary_citation", {})
    exptl = first_or_none(entry.get("exptl", [])) or {}
    refine = first_or_none(entry.get("refine", [])) or {}
    crystal_grow = first_or_none(entry.get("exptl_crystal_grow", [])) or {}
    diffrn = first_or_none(entry.get("diffrn", [])) or {}

    organism, taxonomy_id = get_organism(polymer_entity)
    uniprot_id = get_uniprot_id(polymer_entity)
    structure_title = entry.get("struct", {}).get("title") or ""
    ligand_context = summarize_ligand_context(
        structure_title,
        fetch_nonpolymer_candidates(structure_id, entry),
    )
    oligomeric_state, assembly_type, stoichiometry = infer_oligomeric_state(entry)

    protein_name = entity.get("pdbx_description")
    if protein_name:
        protein_name = protein_name.title()

    gene_name = get_gene_name(polymer_entity)
    related_entities = build_related_entities(gene_name)

    return {
        "structure_registry": "PDB",
        "structure_id": structure_id,
        "external_ids": build_external_ids(structure_id, uniprot_id),
        "pdb_id": structure_id,
        "protein_name": protein_name,
        "gene_name": gene_name,
        "organism": organism,
        "taxonomy_id": taxonomy_id,
        "uniprot_id": uniprot_id,
        "structure_title": structure_title,
        "experimental_method": exptl.get("method") or entry_info.get("experimental_method"),
        "resolution_angstrom": first_or_none(entry_info.get("resolution_combined", [])),
        "ligand_state": ligand_context["ligand_state"],
        "bound_ligands": ligand_context["bound_ligands"],
        "primary_ligand_code": ligand_context["primary_ligand_code"],
        "inhibitor_class": ligand_context["inhibitor_class"],
        "active_site_status": ligand_context["active_site_status"],
        "aged_state": None,
        "oligomeric_state": oligomeric_state,
        "assembly_type": assembly_type,
        "stoichiometry": stoichiometry,
        "chains": polymer_entity.get("rcsb_polymer_entity_container_identifiers", {}).get(
            "auth_asym_ids", []
        ),
        "construct_type": infer_construct_type(polymer_entity),
        "sequence_range": get_sequence_range(polymer_entity),
        "mutations": get_mutations(polymer_entity),
        "glycosylation": None,
        "post_translational_modifications": [],
        "space_group": entry.get("symmetry", {}).get("space_group_name_hm"),
        "crystallization_conditions": crystal_grow.get("pdbx_details"),
        "temperature_k": diffrn.get("ambient_temp"),
        "r_work": refine.get("ls_rfactor_rwork"),
        "r_free": refine.get("ls_rfactor_rfree"),
        "missing_residues": (entry_info.get("deposited_unmodeled_polymer_monomer_count") or 0) > 0,
        "related_entities": related_entities,
        "usable_for_docking": None,
        "usable_for_md": None,
        "primary_reference": citation.get("pdbx_database_id_doi") or citation.get("pdbx_database_id_pub_med"),
        "source_url": f"https://www.rcsb.org/structure/{structure_id}",
        "notes": (
            "Auto-generated from RCSB Data API."
            + (f" {ligand_context['notes_suffix']}" if ligand_context["notes_suffix"] else "")
        ),
    }


def fetch_receptor_structure_record(structure_id: str) -> dict:
    normalized_structure_id = normalize_structure_id(structure_id)
    entry, polymer_entity = fetch_entry_and_polymer_entity(normalized_structure_id)
    return build_receptor_structure_record(normalized_structure_id, entry, polymer_entity)


def write_record_json(record: dict, output_path: str | Path) -> None:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(record, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def write_json(data: dict | list, output_path: str | Path) -> None:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(data, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def ensure_mmcif_file(structure_id: str, output_path: str | Path) -> Path:
    normalized_structure_id = normalize_structure_id(structure_id)
    path = Path(output_path)
    if path.exists():
        return path

    mmcif_text = fetch_text(
        RCSB_MMCIF_DOWNLOAD_URL.format(structure_id=normalized_structure_id)
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(mmcif_text, encoding="utf-8")
    return path


def build_gene_name_query(gene_name: str) -> dict:
    return {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                "operator": "exact_match",
                "value": normalize_gene_name(gene_name),
            },
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": 1000},
        },
    }


def build_uniprot_query(uniprot_id: str) -> dict:
    return {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": uniprot_id,
            },
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": 1000},
        },
    }


def extract_structure_ids_from_search_response(response: dict) -> list[str]:
    structure_ids = []
    for item in response.get("result_set", []):
        identifier = item.get("identifier")
        if not identifier:
            continue
        structure_id = str(identifier).split("_")[0].upper()
        if structure_id and structure_id not in structure_ids:
            structure_ids.append(structure_id)
    return structure_ids


def search_structure_ids_by_gene_name(gene_name: str) -> list[str]:
    response = search_rcsb(build_gene_name_query(gene_name))
    return extract_structure_ids_from_search_response(response)


def search_related_structure_ids(uniprot_id: str | None, gene_name: str | None) -> list[str]:
    if uniprot_id:
        response = search_rcsb(build_uniprot_query(uniprot_id))
        return extract_structure_ids_from_search_response(response)
    if gene_name:
        return search_structure_ids_by_gene_name(gene_name)
    raise RuntimeError("Unable to search related structures without UniProt or gene name")


def extract_identity_from_structure(structure_id: str) -> tuple[str | None, str | None, str | None]:
    _, polymer_entity = fetch_entry_and_polymer_entity(structure_id)
    return (
        get_uniprot_id(polymer_entity),
        get_gene_name(polymer_entity),
        (polymer_entity.get("entity", {}) or {}).get("pdbx_description"),
    )


def is_holo_record(record: dict) -> bool:
    ligands = record.get("bound_ligands") or []
    return bool(ligands) and record.get("ligand_state") != "apo"
