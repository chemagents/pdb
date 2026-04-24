# ORCHESTRATOR.md

## Ontology-first rule

Сначала смысл, потом код.

Цель этого правила: вся логика работы скилла, модулей и подмодулей должна ложиться на онтологию предметной области, а целесообразность каждого действия должна объясняться в терминах онтологии так, чтобы это было понятно человеку без подготовки в программировании.

Из этого следуют четыре требования:

- любая новая логика сначала проверяется на соответствие онтологии, и только потом реализуется в коде;
- если поведение кода расходится с онтологией, приоритет остается за онтологией;
- неоднозначность не скрывается, а явно фиксируется в артефактах и пояснениях;
- качество изменений оценивается не только по работоспособности, но и по смысловой согласованности с предметной областью.

## Status

Этот файл - orchestration document для поэтапной достройки проекта.

Он не является semantic source of truth для сущностей проекта, но является ценным reusable artifact для:

- поэтапной сборки текущего pipeline
- handoff между агентами и сессиями
- повторного использования orchestration patterns в будущих проектах

Когда целевая stage-architecture будет реализована, этот файл может перестать быть ежедневной operational necessity для именно этого репозитория, но не должен автоматически считаться мусором или подлежащим удалению.

## Metaphor

Используй следующую рабочую аналогию:

- `ORCHESTRATOR.md` = ДНК проекта на стадии сборки
- `AGENTS.md` = РНК, то есть operational transcription rules для текущего репозитория
- scripts, schemas, ontology, helpers = белки
- generated and processed files in `output/` = метаболиты

Из этого следуют два правила:

1. `ORCHESTRATOR.md` не должен подменять собой рабочую систему.
2. Любое стабильное содержимое этого файла должно по возможности мигрировать в постоянные структуры проекта, но сам файл может оставаться как reusable construction blueprint.

## Purpose

Текущая цель проекта шире, чем простой extraction одного `receptor_structure`.

Проект должен быть достроен до pipeline, который по `gene_name` или имени белка, например `ACHE`, умеет:

1. найти релевантные holo-структуры
2. консервативно отобрать usable structures
3. выделить подходящие protomer instances
4. привести accepted protomer/ligand instances к общему reference frame
5. собрать consensus ligand cloud
6. разделить альтернативные site clusters
7. построить docking box
8. выпустить QC artifacts и visualization artifacts для human review

## Priority Of Documents

Во всех semantic и repository-local вопросах приоритет имеет:

1. прямой запрос пользователя
2. `AGENTS.md`
3. docs и schemas, на которые ссылается `AGENTS.md`
4. этот `ORCHESTRATOR.md`

Этот файл не имеет права отменять semantic contract из `AGENTS.md`.

Он задает только:

- порядок сборки проекта
- stage model
- artifact model
- growth constraints
- stop conditions for agent-constructor

## Source Intent

При проектировании следующего этапа agent-constructor должен учитывать два разных источника истины:

### Repository truth

`AGENTS.md` описывает, чем проект является сейчас:

- early toolkit
- extraction-first
- schema-driven normalization project
- not yet complete agent platform

### Product / workflow intent

Файл `docs/Автонахождение коробки для докинга.html` описывает, куда проект должен эволюционировать practically:

- docking box не должен строиться по одной случайной структуре
- нужно использовать consensus over multiple holo structures
- meaningful ligands должны отделяться от additives and artifacts
- alignment, clustering, box estimation и QC должны быть explicit stages
- visual overlay accepted ligands in one reference frame является обязательным QC layer

## What This File Is For

Этот файл существует для agent-constructor, который достраивает проект stage-by-stage.

Он должен помогать агенту:

- не превращать проект в giant script
- не строить black-box pipeline
- не расползаться в чрезмерный refactor
- оставлять reviewable artifacts after every stage
- двигаться минимальными, завершенными инкрементами

## What This File Must Not Become

Этот файл не должен становиться:

- копией `AGENTS.md`
- копией design discussion из HTML-диалога
- backlog без executable consequences
- permanent source of truth

Если какая-то идея уже стабилизировалась, она должна быть перенесена:

- в `AGENTS.md`, если это permanent repo rule
- в код, если это executable behavior
- в schemas / ontology, если это semantic contract
- в output manifests, если это runtime traceability rule

## Construction Principles

Проект надо строить как наблюдаемый pipeline из маленьких стадий.

Предпочтительная модель:

1. один входной workflow
2. несколько stage scripts
3. явные stage outputs
4. явные QC outputs
5. возможность rerun одной стадии
6. отсутствие скрытого состояния в голове модели

Каждый полезный следующий шаг должен улучшать одно из трех:

- correctness
- observability
- resumability

## Required Separation Of Responsibilities

### Code is responsible for

- reading mmCIF / PDB / RCSB data
- biological assembly reconstruction
- chain and protomer extraction
- sequence mapping and residue mapping
- structure superposition
- coordinate transformation
- contact calculation
- ligand cloud aggregation
- clustering
- box geometry
- JSON writing

### LLM is responsible for

- designing next minimal stage
- conservative classification in ambiguous cases
- deciding what to reject vs keep when the logic is not purely mechanical
- writing summaries and QC explanations
- checking whether outputs satisfy project intent and repo rules

### LLM must not be the sole source of truth for

- chain mapping
- transforms
- atom coordinates
- distances and contact counts
- final box coordinates

## Stage Model

The target pipeline should emerge in the following order.

The exact file names may evolve, but the stage boundaries should remain explicit.

### Stage 0. Candidate Acquisition

Goal:

- fetch related candidate structures for a gene / target

Minimum outputs:

- candidate index
- per-structure receptor records

Current state:

- partially exists already in `fetch_receptor_structure.py`, `fetch_related_holo_receptor_structures.py`, `fetch_receptor_by_gene_name.py`

### Stage 1. Structure Selection

Goal:

- decide which candidate holo structures are usable for downstream site inference

Input rule:

- this stage must operate on already collected candidate structure description files
- it must not re-fetch candidates from the network if pre-fetched `receptor_structure` JSON files already exist
- candidate acquisition and candidate selection are distinct stages and should remain distinct in code and artifacts

Selection should conservatively consider:

- organism / target identity
- holo relevance
- meaningful ligands vs additives
- structure quality signals
- obvious construct/pathology issues

Minimum outputs:

- accepted / rejected / uncertain decision log
- reasons for each decision

### Stage 2. Protomer Extraction

Goal:

- move from structure-level records to biologically relevant protomer instances

This is the first stage where coordinate-bearing structure files should normally become explicit local artifacts.

At this stage the pipeline should, at minimum for accepted structures:

- materialize local source structure files, preferably `mmCIF`
- link protomer artifacts back to those local coordinate sources
- avoid postponing all coordinate materialization until alignment if protomer work already depends on structure-level geometry context

This stage must make explicit:

- which assembly context is used
- which chain / protomer is considered relevant
- which ligand instances belong to which protomer context

Minimum outputs:

- protomer instance records
- ligand instance records linked to protomer context
- local coordinate-bearing source files for accepted structures, unless they already exist in a stable upstream artifact location

### Stage 3. Canonical Alignment

Goal:

- define one reference protomer and align accepted protomer instances into its frame

If the workspace already contains a human-authored referent artifact for one of the candidate structures, this stage should prefer that structure as the base reference contract over a purely heuristic reference choice.

This stage must emit:

- chosen reference protomer
- mapping / transform metadata
- alignment QC

### Stage 4. Ligand Cloud And Site Clustering

Goal:

- transfer accepted ligand instances to reference frame
- build consensus ligand cloud
- split alternative site clusters

This stage must not silently merge:

- monomer-local sites
- interface sites
- unrelated outliers

### Stage 5. Box Construction

Goal:

- build docking box from accepted consensus cloud

This stage must explain:

- which ligands were used
- what trimming was applied
- what margin was applied
- why this cluster was chosen

### Stage 6. Visualization And Review

Goal:

- create artifacts that let a human quickly validate whether the pipeline is sane

Minimum desired visuals:

- reference protomer + all accepted ligands
- reference protomer + rejected or outlier ligands
- reference protomer + final box
- interactive 3D scene or script that renders the protein, aligned ligand bouquet, and the automatically inferred docking box
- 2D comparison render that overlays the automatically inferred box with a human reference box when a referent config exists for the chosen reference structure

These visualization artifacts are not optional polish.

They are part of the required observable output of the pipeline and must be generated automatically by the visualization stage.

## Artifact Model

Every substantial stage should leave artifacts that can be inspected without rerunning previous stages.

Preferred artifact classes:

- stage index / manifest
- per-structure or per-protomer records
- QC JSON
- human-readable markdown summary
- visualization artifacts

Artifacts should be:

- deterministic
- easy to diff
- sufficient for handoff to another agent or human

## Growth Constraints

Agent-constructor should prefer the smallest next stable step.

Prefer:

- one new stage script over one giant orchestrator rewrite
- one minimal helper module when reuse becomes real
- one output manifest over hidden control flow
- one conservative classification rule over speculative biological inference

Avoid until clearly justified:

- package-wide refactor
- heavy framework introduction
- broad abstraction layers before second reuse
- hidden orchestration state

## Minimal Increment Rule

When deciding what to implement next, prefer a step that satisfies all of the following:

1. creates one new durable capability
2. leaves reviewable artifacts
3. can be tested on one known gene or one known PDB set
4. does not require speculative semantics
5. does not require rewriting previous stages

## Stage Gate Protocol

Before implementing or materially changing any stage, agent-constructor must perform an explicit stage-gate check.

The check is mandatory and must answer all of the following before code is written:

1. What is the immediately previous stage?
2. What are the expected input artifacts of the current stage?
3. Do those input artifacts already exist in `output/`?
4. If they exist, can the stage operate purely from those local artifacts?
5. If the implementation still wants to fetch or recompute upstream data, is that actually part of the current stage boundary?

If the current stage can operate from existing local artifacts, the default behavior must be:

- prefer local artifacts
- avoid re-fetching upstream data
- avoid silently merging current stage responsibilities with earlier stages

This rule exists specifically to prevent accidental boundary violations such as:

- selection stage re-fetching candidate records
- downstream stages recomputing upstream curation without need
- stage scripts hiding missing orchestration logic by calling the network again

## Autonomous Build Requirement

The project must be built so that the next working stage emerges autonomously from the documented orchestration model, not only after user correction.

This means agent-constructor must not rely on the user to point out obvious stage-boundary violations that were already decidable from:

- `AGENTS.md`
- this `ORCHESTRATOR.md`
- existing artifacts in `output/`

When multiple plausible implementations exist, the default autonomous choice should be the one that:

- respects stage boundaries most strictly
- maximizes reuse of already generated artifacts
- minimizes hidden recomputation
- leaves the clearest audit trail

## Human Review Stop Conditions

Stop and ask for review if any of the following becomes materially important:

- meaningful ligand identity is ambiguous
- chain / protomer equivalence is unclear
- assembly choice changes site interpretation
- interface vs monomer-local classification is unresolved
- multiple site clusters compete with no conservative winner
- final box is dominated by one or two outliers
- visual QC contradicts numerical QC
- implementing the next step would require repository-wide reorganization rather than local staged growth

## Migration Rule

Whenever a piece of orchestration guidance becomes stable, migrate it out of this file where appropriate.

Examples:

- stable permanent repository rule -> `AGENTS.md`
- stable executable workflow -> script or helper module
- stable semantic field contract -> schema or ontology
- stable output traceability pattern -> generated manifest structure in code

After migration, this file may:

- keep a shorter summary and references
- keep reusable stage patterns
- keep generic construction heuristics useful for future related projects

## Long-Term Role

If the project matures, `ORCHESTRATOR.md` should evolve from a highly project-specific build note into one of these stable roles:

1. reusable construction blueprint for related agent-built projects
2. historical architecture record explaining why the pipeline was staged this way
3. high-level orchestration guide complementary to `AGENTS.md`

Therefore, do not delete this file by default.

Delete or archive it only if the user explicitly decides that:

- its contents have been fully superseded
- its reuse value is low
- or a better stable replacement document already exists
