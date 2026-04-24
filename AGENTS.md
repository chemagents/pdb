# AGENTS.md

## Project Purpose

Этот репозиторий - ранняя toolkit-заготовка для agentic workflows, которые извлекают и нормализуют данные о `receptor_structure`, `binding_site`, `ligand_complex` и `ligand` из `rcsb.org` и связанных registry.

Агент должен воспринимать проект так:

- в первую очередь это data-extraction project
- во вторую очередь это schema-driven normalization project
- в будущем это может стать agent platform, но сейчас это еще не полноценная система

Главная цель проекта - помогать человеку и coding agent получать строгие, повторяемые, schema-aligned JSON records для основных сущностей проекта.

Главный source of truth для семантики находится в `docs/` и `schemas/`, а не в текущем поведении отдельных script.

## Directory Map

- `README.md`: обзор проекта на русском языке
- `docs/entity_schema_readme.ru.md`: русское руководство по чтению и применению schema
- `docs/entity_schema_readme.en.md`: английское руководство по чтению и применению schema
- `schemas/*.schema.*.yaml`: schema-файлы сущностей для normalized storage
- `schemas/domain_ontology.master.*.yaml`: ontology index и связи между сущностями
- `tools/`: рабочая зона для extractors, normalization scripts и utility code
- `tools/fetch_receptor_structure.py`: получает одну receptor structure из RCSB и пишет normalized JSON
- `tools/fetch_related_holo_receptor_structures.py`: ищет связанные human/mouse holo structures и переиспользует base fetcher
- `tools/select_holo_structures.py`: stage 1 selection script, который читает pre-fetched `receptor_structure` records и пишет selection artifacts
- `tools/validate_json_against_schema.py`: lightweight validator для generated JSON records against schema field names, required fields, enums и basic types
- `output/`: generated JSON outputs на одном уровне с `AGENTS.md`; внутри records распределяются по папкам с `<gene_name>` для отдельных белков

`tools/` - это workspace для extraction и normalization logic.

`output/` содержит generation results и должен рассматриваться как generated artifacts, а не source code. Эта папка находится на одном уровне с `AGENTS.md`, а generated records должны складываться в gene-scoped subdirectories вроде `output/ACHE/`.

## How To Read The Project

Когда агент начинает extraction или normalization, он должен идти в таком порядке:

1. Сначала прочитать `docs/entity_schema_readme.ru.md` или `docs/entity_schema_readme.en.md`.
2. Затем выбрать target schema в `schemas/`.
3. Рассматривать `fields`, enum-like `allowed_values` и `interpretation_rules` как semantic contract.
4. Формировать JSON, совместимый с `output_template` выбранной schema.
5. Представлять cross-entity links через `related_entities`.

Нельзя начинать с догадок о смысле полей только по существующему поведению script.

## Entity Workflow

Ожидаемый entity workflow в этом репозитории:

1. определить, какая сущность извлекается: `receptor_structure`, `binding_site`, `ligand_complex` или `ligand`
2. прочитать соответствующий schema-file в `schemas/`
3. извлекать только то, что подтверждается trusted source data
4. нормализовать значения к точным field names и enum values из schema
5. сохранять результат как JSON, совместимый с `output_template`
6. перед downstream use по возможности валидировать generated JSON against target schema
7. связывать связанные records через `related_entities`

Текущие роли сущностей:

- `receptor_structure`: одна normalized запись о receptor/protein structure
- `binding_site`: одна normalized residue-level интерпретация сайта связывания, привязанная к конкретному receptor_structure context
- `ligand_complex`: одна normalized интерпретация receptor-ligand complex
- `ligand`: одна normalized chemical ligand entity

Отдельно в pipeline уже отражена stage `structure_selection`:

- это не новая ontology entity, а operational artifact layer поверх набора `receptor_structure` records
- она читает уже собранные candidate `receptor_structure` JSON files и не должна повторно делать candidate acquisition, если эти records уже существуют
- ее задача - фиксировать conservative accepted / uncertain / rejected решения для downstream docking-oriented workflow
- она должна оставлять reviewable artifacts: per-structure decision JSON, stage index JSON и human-readable summary

## Rules For Extraction And Normalization

Эти правила обязательны для агентов, работающих в репозитории:

- читать `docs/entity_schema_readme.ru.md` или `.en.md` до изменения extraction semantics
- выбирать target schema до извлечения значений
- сохранять field names точно как в schema
- сохранять enum-like values точно как в schema
- считать `fields`, `allowed_values` и `interpretation_rules` semantic contract
- использовать `null`, если значение отсутствует или не подтверждается trusted source
- не делать silent guessing
- сохранять provenance через registry и identifier fields
- строго сохранять entity boundaries и не смешивать `receptor_structure`, `binding_site`, `ligand_complex` и `ligand`
- сохранять output как JSON, совместимый с соответствующим `output_template`
- связывать сущности через `related_entities`
- относиться консервативно к solvent, salts, buffers и другим non-functional additives, если schema или source явно не требуют иного
- помнить, что текущая реализация `apo`/`holo` остается эвристической и может требовать curator review
- использовать `primary_ligand_code` только когда один ключевой ligand code можно выбрать консервативно и аудируемо
- формировать `bound_ligands` так, чтобы главный ligand по возможности попадал в этот candidate set для downstream восстановления координат protein-ligand binding region

## Current Gaps

В репозитории уже есть рабочие schema-files и functional extractor `fetch_receptor_structure.py`, но это еще не complete agent pipeline.

Текущие gaps:

- нет unified CLI или single entry point для всех workflow
- пока нет extractor или curator workflow для `binding_site`
- пока нет extractor для `ligand_complex`
- пока нет extractor для `ligand`
- normalization logic в текущих scripts неполная
- логика разделения `apo`/`holo` и выделения primary ligand code пока несовершенна
- нет formal integration в более широкий agent pipeline
- нет dedicated automated test suite

Уточнение по current state:

- отдельная stage `structure_selection` уже реализована в `tools/select_holo_structures.py`, но downstream geometry stages после нее пока не реализованы полностью
- lightweight JSON validation against schema expectations теперь доступна отдельным script-ом, но полноценной declarative validation framework пока нет

Агент не должен переоценивать зрелость проекта.

Важно: текущий `receptor_structure` extractor уже использует RCSB nonpolymer entity data, а не только `nonpolymer_bound_components`, чтобы лучше находить candidate ligand codes. Но это еще не полностью решает semantic задачу выделения одного ключевого лиганда и уверенного разделения `apo`/`holo` во всех случаях.

Для текущей семантики `receptor_structure` подразумевается следующее:

- `bound_ligands` хранит только meaningful ligand component codes после conservative filtering
- `bound_ligands` также служит practical candidate set для downstream binding-site coordinate inference; если в структуре есть главный ligand, он по возможности должен присутствовать в этом списке
- `primary_ligand_code` хранит один functionally central ligand code, если его можно выбрать без guessing
- если meaningful ligands несколько и явного winner нет, `primary_ligand_code` должен оставаться null
- peptide/protein-bound inhibitory contexts не должны автоматически схлопываться в `apo`
- ambiguity должна отражаться в `notes`, а не маскироваться уверенным enum

## Recommended Next Tasks

Полезные следующие задачи:

- добавить extractor для `ligand_complex` в `tools/`
- добавить extractor для `ligand` в `tools/`
- добавить extractor или curator workflow для `binding_site`
- усилить JSON validation against schema expectations и встроить ее глубже в generation workflow
- улучшить holo vs apo discrimination
- улучшить выделение primary ligand code для receptor structures
- улучшить normalization coverage в текущем `receptor_structure` extractor
- добавить summary или index outputs для bulk runs
- ввести unified CLI, когда несколько extractors стабилизируются
- добавить lightweight automated tests для известных PDB examples

## External Rules

На момент написания в репозитории не найдено repository-specific Cursor rules или Copilot instructions.

Не найдены:

- `.cursor/rules/`
- `.cursorrules`
- `.github/copilot-instructions.md`

Если такие файлы появятся позже, обнови этот документ и считай их более приоритетными repository rules.

## Environment Assumptions

- текущая platform: Windows / PowerShell
- основной language в проекте: Python
- сейчас нет `pyproject.toml`, `requirements.txt`, `pytest.ini`, `tox.ini` или `setup.cfg`
- formal build system пока не настроен
- dedicated test suite пока нет

Агент не должен придумывать отсутствующий tooling и описывать его как уже существующую инфраструктуру.

## Build, Lint, and Test Commands

Так как репозиторий легковесный, практические команды сейчас - это запуск script и syntax validation.

### Setup

Используй active Python interpreter из окружения.

Примеры:

```bash
python --version
python -m pip --version
```

Проверить один generated JSON against schema:

```bash
python tools/validate_json_against_schema.py output/ACHE/receptor_structure_pdb_4EY7.json
```

### Run Main Scripts

Получить одну normalized receptor structure:

```bash
python tools/fetch_receptor_structure.py 4EY7 --output-dir output
```

Получить связанные human и mouse holo structures для того же receptor:

```bash
python tools/fetch_related_holo_receptor_structures.py 4EY7 --output-dir output
```

Отобрать usable holo structures из уже собранных `receptor_structure` records:

```bash
python tools/select_holo_structures.py ACHE --output-dir output
```

### Build / Syntax Check

Formal build step отсутствует. В качестве базовой проверки syntax используй Python bytecode compilation.

Проверить весь `tools/`:

```bash
python -m compileall tools
```

Проверить один script:

```bash
python -m py_compile tools/fetch_receptor_structure.py
python -m py_compile tools/fetch_related_holo_receptor_structures.py
python -m py_compile tools/select_holo_structures.py
python -m py_compile tools/validate_json_against_schema.py
```

### Lint

Отдельный linter в репозитории не настроен.

Если нужна легкая проверка без изменения project configuration, используй:

```bash
python -m py_compile tools/fetch_receptor_structure.py
```

Если позже появится linter, нужно явно задокументировать точную команду.

### Test

Automated test suite пока нет.

Текущая практическая verification - это запуск script на известных PDB identifier и проверка generated JSON files.

Рекомендуемые manual verification commands:

```bash
python tools/fetch_receptor_structure.py 5U09 --output-dir output
python tools/fetch_receptor_structure.py 4EY7 --output-dir output
python tools/fetch_related_holo_receptor_structures.py 4EY7 --output-dir output
python tools/select_holo_structures.py ACHE --output-dir output
python tools/validate_json_against_schema.py output/ACHE/receptor_structure_pdb_4EY7.json
```

### Single Test Guidance

Formal single-test command пока нет, потому что `pytest` tests в проекте отсутствуют.

До появления test suite ближайший эквивалент single test - это запуск одного script для одного известного structure ID, например:

```bash
python tools/fetch_receptor_structure.py 4EY7 --output-dir output
```

Если в будущем появится `pytest` suite, single test command должен выглядеть так:

```bash
pytest path/to/test_file.py::test_name
```

Не утверждай, что `pytest` доступен, если он реально не добавлен в репозиторий.

## Source of Truth for Data Semantics

При изменении extraction или normalization logic сначала читай:

- `docs/entity_schema_readme.ru.md`
- `docs/entity_schema_readme.en.md`
- `schemas/receptor_structure.schema.en.yaml`
- `schemas/receptor_structure.schema.ru.yaml`
- `schemas/binding_site.schema.en.yaml`
- `schemas/binding_site.schema.ru.yaml`
- `schemas/ligand_complex.schema.en.yaml`
- `schemas/ligand.schema.en.yaml`
- `schemas/domain_ontology.master.en.yaml`

Ключевые semantic rules, уже зафиксированные там:

- сохранять точные field names
- использовать `null`, а не guessing
- держать provenance явным
- сохранять строгие entity boundaries
- не смешивать концептуально разные entities
- не считать solvent, salts или buffers meaningful ligands без явного основания
- не считать текущую эвристическую apo/holo классификацию полностью завершенной

Ключевые operational rules для `structure_selection` stage:

- эта стадия должна работать по pre-fetched `receptor_structure` JSON artifacts, а не заменять собой candidate acquisition
- результаты selection должны оставаться inspectable и rerunnable как отдельные stage outputs
- selection decisions должны фиксировать причины accepted / uncertain / rejected статуса явно, а не только итоговый shortlist

## Code Style

### General

- держать изменения маленькими и локальными
- предпочитать straightforward standard-library solutions
- не добавлять frameworks или heavy dependencies без необходимости
- сохранять текущую script-oriented architecture, если нет явной причины на repository-wide refactor
- держать generated data behavior согласованным со schema, даже если extractor пока неполный

### Imports

- группировать imports в таком порядке: standard library, затем local imports
- использовать one import per line, если grouping не делает код заметно яснее
- предпочитать explicit imports из local modules
- не добавлять unused imports

Пример текущего стиля:

```python
import argparse
import json
import sys
from pathlib import Path

from fetch_receptor_structure import fetch_json
```

### Formatting

- следовать PEP 8
- использовать 4-space indentation
- держать строки читаемыми, без лишней горизонтальной плотности
- использовать trailing commas в multiline literals, если это улучшает diff
- предпочитать простой control flow вместо плотных comprehensions, если логика становится хуже для аудита

### Types

- добавлять type hints для новых functions
- сохранять текущий стиль type hints, например `str | None`, `list[str]`, `tuple[str | None, int | None]`
- если helper возвращает heterogeneous JSON-like data, можно использовать `dict`, если более строгая typing только усложнит код
- не вводить избыточные typing abstractions без реальной пользы

### Naming

- functions: `snake_case`
- variables: `snake_case`
- constants: `UPPER_SNAKE_CASE`
- file names: `snake_case.py`
- использовать имена, отражающие biological/data semantics
- предпочитать термины, уже используемые в schema, например `structure_id`, `gene_name`, `organism`, `bound_ligands`

### Function Design

- держать helpers сфокусированными на одном шаге extraction или normalization
- по возможности предпочитать pure transformation helpers
- переиспользовать существующие utilities вместо дублирования API-fetch logic
- если один script зависит от другого, предпочитать import reusable functions или вызывать существующий script только когда это действительно часть design
- если в workspace присутствует явный human-authored reference artifact для одной из candidate structures, alignment stage должен предпочитать этот reference contract эвристическому выбору

### Error Handling

- обрабатывать `urllib.error.HTTPError` и `urllib.error.URLError` явно для network operations
- падать с понятным stderr message, включающим `structure_id`
- возвращать non-zero exit code из CLI entry points при ошибке
- использовать `RuntimeError` для logical failures в data acquisition или transformation
- не проглатывать exceptions молча
- не допускать тихой порчи `output/`; если script удаляет файлы при filtering, это должно быть намеренно и ясно

### CLI Conventions

- использовать `argparse`
- давать короткий, но специфичный `description`
- принимать `structure_id` как positional argument для fetch scripts
- использовать `--output-dir` для переопределения output location
- по умолчанию направлять output directories в repository `output/` на одном уровне с `AGENTS.md`, а для отдельных белков сохранять records в поддиректории `output/<gene_name>/`

### Network and API Use

- использовать RCSB Data API и Search API через standard library, пока проект формально не перешел на другой client
- устанавливать `User-Agent`
- устанавливать timeout для network calls
- нормализовать structure codes до использования
- считать network failures нормальными и recoverable на уровне process

### JSON Output

- писать UTF-8 с `ensure_ascii=False` и завершающей newline
- по возможности сохранять stable field ordering в порядке schema
- использовать предсказуемые output filenames, например `receptor_structure_pdb_<CODE>.json`
- не добавлять fields, которых нет в schema, если schema специально не расширена

### Visualization Outputs

Для pipeline автоматического определения docking box визуализация является обязательной частью stage outputs.

После построения box pipeline должен по возможности автоматически выпускать:

- интерактивный 3D artifact или Python script для его построения, где видны protein context, ligand bouquet и автоматически найденный docking box
- 2D overlay render для human review
- если для reference structure существует human-authored referent config, отдельный 2D comparison render с наложением автоматически найденного box и human reference box

Эти artifacts должны собираться автоматически внутри visualization stage, а не вручную по просьбе пользователя после основного прогона.

### Schema Alignment

- считать schema-files semantic contract
- точно совпадать с enum-like values
- использовать `null` для unknown values вместо выдуманных placeholder
- сохранять one record per declared granularity unit
- отражать uncertainty в `notes`, когда это нужно

### Output Directory Rules

- `output/` содержит generated artifacts, а не source code
- новые scripts, генерирующие normalized records, должны продолжать заполнять `output/` на одном уровне с `AGENTS.md`
- records для отдельных белков должны сохраняться в gene-scoped subdirectories внутри `output/`, например `output/ACHE/`
- не удалять unrelated files из `output/`
- если script удаляет generated files в рамках filtering, он должен удалять только файлы, созданные этой же логикой
- при использовании unified staged pipeline итоговый агрегированный summary artifact должен сохраняться как `output/<gene_name>/07_exports/site_pipeline_summary.json`

## Agent Change Strategy

При внесении изменений:

- сначала инспектировать `docs/` и `schemas/` перед изменением semantics
- держать extractor behavior консервативным
- предпочитать улучшение correctness вместо расширения speculative metadata
- проверять изменения запуском relevant script на одном известном PDB code
- явно сообщать, если validation была только manual из-за отсутствия tests
- перед реализацией новой stage-script логики сначала проверить stage inputs и existing generated artifacts в `output/`
- если текущая стадия по смыслу должна работать на уже собранных JSON artifacts, предпочитать чтение локальных файлов повторному network fetch
- не смешивать acquisition, selection, protomer extraction, alignment и box construction в одной реализации только потому, что это проще закодировать
- если stage boundary неясна, остановиться, перечитать `ORCHESTRATOR.md` и existing artifacts, а не угадывать поведение

## Base Reference Rule

Для downstream alignment stage reference protomer должен определяться не по хардкоженному списку gene names, а по явным артефактам текущего workflow.

Предпочтительный порядок:

1. explicit local reference artifact, если такой создан в stage outputs
2. human-authored referent config, если он существует для одной из candidate structures текущего target set
3. deterministic stage-derived canonical candidate

Практическое правило:

- если для структуры существует human-authored referent config и normalized `receptor_structure` record показывает monomeric single-chain context, этот structure+chain следует рассматривать как base protomer candidate для alignment
- если normalized `receptor_structure` record показывает monomeric context, но содержит несколько symmetry-equivalent chain IDs, agent должен выбирать базовый chain по детерминированному правилу, а не отбрасывать referent-backed structure полностью
- если referent-backed candidates несколько, agent должен либо выбирать их по детерминированному правилу с понятным provenance, либо останавливать workflow для curator review, если неоднозначность materially affects downstream geometry

## What Not To Do

- не считать, что в проекте уже есть linting или tests, если их нет
- не вводить schema-breaking JSON keys без явной необходимости
- не считать все nonpolymer components biologically meaningful ligands без filtering logic
- не хардкодить organism-specific assumptions сверх уже задокументированного поведения
- не переписывать репозиторий в package layout без явного запроса

## Recommended Next Improvements

Полезные future improvements:

- stronger JSON validation against schema-files and stage artifacts
- lightweight test suite для известных PDB IDs
- shared helpers для RCSB search и normalization
- summary/index files для bulk extraction runs
- better holo vs apo discrimination
