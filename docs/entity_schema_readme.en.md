# Entity Schema Guide

These schema files exist so that humans and agent systems share a strict, reusable understanding of sensitive real-world modeling objects. They also serve as building blocks for a domain ontology.

## Purpose

Use these schemas when an agent needs to:

- extract structured data from registries such as PDB, UniProt, AlphaFold, PubChem, or related sources
- normalize records across multiple structure and chemistry registries
- preserve consistent semantics for receptor structures, ligand complexes, and ligands
- preserve consistent semantics for receptor structures, binding sites, ligand complexes, and ligands
- avoid silent guessing in ambiguous cases
- build ontology-aligned datasets for docking, molecular simulation, annotation, or knowledge integration

## Files

- `schemas/receptor_structure.schema.en.yaml`
- `schemas/receptor_structure.schema.ru.yaml`
- `schemas/ligand_complex.schema.en.yaml`
- `schemas/ligand_complex.schema.ru.yaml`
- `schemas/ligand.schema.en.yaml`
- `schemas/ligand.schema.ru.yaml`
- `schemas/binding_site.schema.en.yaml`
- `schemas/binding_site.schema.ru.yaml`
- `schemas/domain_ontology.master.en.yaml`
- `schemas/domain_ontology.master.ru.yaml`

## Agent Reading Rules

When an agent reads a schema file, it must treat the schema as a semantic contract.

The agent should:

1. identify the entity type and record granularity before extracting any values
2. preserve field names exactly as defined in the schema
3. obey allowed values where enumerations are provided
4. use `null` when reliable data is unavailable
5. avoid inferring values that are not supported by trusted sources
6. preserve registry provenance using fields such as `structure_registry`, `structure_id`, `external_ids`, and `registry_ids`
7. represent cross-entity links through `related_entities`
8. keep one record per declared granularity unit

## Normalization Rules

- Do not assume PDB is the only source of structural truth.
- Do not merge distinct chemical states if the schema treats them as separately meaningful.
- Do not treat buffer components, salts, solvent molecules, or cryoprotectants as meaningful ligands unless the source or curator explicitly identifies them as such.
- Prefer biologically relevant assembly semantics over purely deposited geometry when the source provides both.
- Keep receptor structure, binding site, ligand complex, and ligand as separate entities even when they are tightly related.

Practical note: current heuristic apo/holo separation and primary ligand code selection may still be incomplete. When a record contains multiple nonpolymer components, the agent should conservatively distinguish the key functional ligand from additives, ions, glycans, and other non-functional components.

For `receptor_structure`, treat `bound_ligands` as the filtered list of meaningful ligand component codes, and treat `primary_ligand_code` as a stricter field: it should only be populated when one functionally central ligand can be selected conservatively. If several meaningful ligands remain or the ranking is unclear, keep `primary_ligand_code` null and explain the ambiguity in `notes`.

For `binding_site`, treat `residue_members` as the conservative residue set for the site interpretation, and treat `active_residue_members` as a stricter subset. Do not automatically promote all ligand-contacting residues to active residues. If the source only supports ligand-contact or software-generated binding-site evidence, keep the site interpretation conservative and avoid catalytic guessing.

## Suggested Workflow For Agents

1. choose the target entity schema
2. read its `description`, `entity_definition`, and `record_granularity`
3. inspect `fields` and create a field-by-field extraction checklist
4. inspect `interpretation_rules` before resolving ambiguous cases
5. produce output that matches `output_template`
6. record unresolved ambiguity in `notes`
7. link related objects using `related_entities`

## Minimal Output Discipline

- keep exact field names
- keep types stable
- keep enums stable
- keep provenance explicit
- keep uncertainty visible

## Ontology Use

These schemas are suitable as practical ontology seeds. They define:

- entity boundaries
- field semantics
- cross-entity relations
- normalization expectations
- registry linkage patterns

They can later be translated into JSON Schema, SHACL, OWL, RDF mappings, database models, or agent memory rules.
